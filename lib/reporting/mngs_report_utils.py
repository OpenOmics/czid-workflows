"""Standard-library reporting shared by the short- and long-read workflows."""

import argparse
import csv
import gzip
import json
import math
from datetime import datetime


# PipelineReportService::CSV_SHORT_READS_COLUMNS in chanzuckerberg/czid-web.
SHORT_READ_COLUMNS = [
    "tax_id", "tax_level", "genus_tax_id", "name", "common_name", "category",
    "is_phage", "agg_score", "max_z_score",
    "nt_z_score", "nt_rpm", "nt_count", "nt_contigs", "nt_contig_r",
    "nt_percent_identity", "nt_alignment_length", "nt_e_value", "nt_bg_mean",
    "nt_bg_stdev", "nt_bg_mean_mass_normalized", "nt_bg_stdev_mass_normalized",
    "nr_z_score", "nr_rpm", "nr_count", "nr_contigs", "nr_contig_r",
    "nr_percent_identity", "nr_alignment_length", "nr_e_value", "nr_bg_mean",
    "nr_bg_stdev", "nr_bg_mean_mass_normalized", "nr_bg_stdev_mass_normalized",
    "species_tax_ids", "known_pathogen",
]

LONG_READ_COLUMNS = [
    "tax_id", "tax_level", "genus_tax_id", "name", "common_name", "category", "is_phage",
    "nt_bpm", "nt_base_count", "nt_count", "nt_contigs", "nt_contig_b",
    "nt_percent_identity", "nt_alignment_length", "nt_e_value",
    "nr_bpm", "nr_base_count", "nr_count", "nr_contigs", "nr_contig_b",
    "nr_percent_identity", "nr_alignment_length", "nr_e_value", "species_tax_ids", "known_pathogen",
]


def add_contig_metrics(rows, contig_summary_json, contig_stats_json, contigs_fasta,
                       technology="short-read", contig_base_counts_json=None):
    """
    Adapt PipelineRun's contig loading and PipelineReportService.merge_contigs.

    The platform resolves species to genus using its taxonomy database. Here the
    pipeline summary supplies both assignments, avoiding a new reference input.
    Read totals come from assembly statistics, not the taxon-summary counts.
    """
    with open(contig_summary_json) as source:
        summary = json.load(source)
    with open(contig_stats_json) as source:
        read_counts = json.load(source)
    base_counts = None
    if technology == "long-read":
        with open(contig_base_counts_json) as source:
            base_counts = json.load(source)
    min_reads = 1 if technology == "long-read" else 4
    metric = "contig_b" if technology == "long-read" else "contig_r"
    with open(contigs_fasta) as source:
        contig_names = {line[1:].rstrip() for line in source if line.startswith(">")}

    assignments = {}
    for entry in summary:
        count_type = entry["count_type"]
        level = int(entry["tax_level"])
        if count_type not in ("NT", "NR") or level not in (1, 2):
            continue
        tax_id = int(entry["taxid"])
        for name in entry["contig_counts"]:
            if name == "*" or name not in contig_names:
                continue
            lineage = assignments.setdefault((name, count_type), {})
            if level in lineage and lineage[level] != tax_id:
                raise ValueError("Conflicting contig assignments for {} in {}".format(name, count_type))
            lineage[level] = tax_id

    for (name, count_type), lineage in assignments.items():
        # CZ ID loads assignments from species records; genus-only entries do
        # not independently assign a contig. The pipeline already filters on
        # unique supporting reads; apply the web loader's technology-specific check.
        reads = read_counts.get(name, 0)
        if 1 not in lineage or reads < min_reads:
            continue
        prefix = count_type.lower()
        for level, tax_id in lineage.items():
            row = rows.get((level, tax_id))
            if row is None or prefix + "_count" not in row:
                continue
            row[prefix + "_contigs"] = row.get(prefix + "_contigs", 0) + 1
            contribution = base_counts.get(name, 0) if base_counts is not None else reads
            row[prefix + "_" + metric] = row.get(prefix + "_" + metric, 0) + contribution


CATEGORIES = {"2": "bacteria", "2157": "archaea", "2759": "eukaryota",
              "10239": "viruses", "12884": "viroids", "": "uncategorized"}
INVALID_CALL_BASE_ID = -100000000


def add_taxon_information(rows, lineage_csv, lineage_version):
    """
    Join only the requested taxa from the versioned CZ ID reference.
    Select records for the explicitly requested taxonomy snapshot.
    """
    datetime.strptime(lineage_version, "%Y-%m-%d")
    wanted = set()
    for row in rows.values():
        for tax_id in (row["tax_id"], row["genus_tax_id"]):
            wanted.add(tax_id if tax_id >= INVALID_CALL_BASE_ID else abs(tax_id) % abs(INVALID_CALL_BASE_ID))
    lineages = {}
    version_available = False
    opener = gzip.open if str(lineage_csv).endswith(".gz") else open
    with opener(lineage_csv, "rt", newline="") as source:
        reader = csv.DictReader(source)
        required = {"taxid", "version_start", "version_end", "species_name", "genus_name",
                    "species_common_name", "genus_common_name", "family_name", "superkingdom_taxid", "is_phage"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Taxonomy CSV is missing columns: {}".format(sorted(required - set(reader.fieldnames or []))))
        for entry in reader:
            if not entry["version_start"] <= lineage_version <= entry["version_end"]:
                continue
            version_available = True
            tax_id = int(entry["taxid"])
            if tax_id in wanted:
                if tax_id in lineages:
                    raise ValueError("Overlapping lineage versions for taxon {}".format(tax_id))
                lineages[tax_id] = entry
    if not version_available:
        raise ValueError("No taxonomy records cover lineage version {}".format(lineage_version))

    for row in rows.values():
        tax_id = row["tax_id"]
        rank = "species" if row["tax_level"] == 1 else "genus"
        entry = lineages.get(tax_id, {})
        category_entry = entry
        if tax_id < INVALID_CALL_BASE_ID:
            parent_id = abs(tax_id) % abs(INVALID_CALL_BASE_ID)
            category_entry = lineages.get(parent_id, {})
            parent_rank = "genus" if rank == "species" else "family"
            if rank == "species":
                parent_id = row["genus_tax_id"]
            parent = lineages.get(parent_id, {}).get(parent_rank + "_name") or "taxon {}".format(parent_id)
            name = "non-{}-specific reads in {} {}".format(rank, parent_rank, parent)
        elif tax_id < 0:
            name = "all taxa with neither family nor genus classification"
            if tax_id == -201:
                name = "all artificial constructs"
            elif tax_id not in (-1, -100, -200, -300, -400, -500, -600, -650, -700):
                name += " {}".format(tax_id)
        else:
            name = entry.get(rank + "_name") or "unnamed {} taxon {}".format(rank, tax_id)
        row["name"] = name
        row["common_name"] = entry.get(rank + "_common_name", "")
        row["category"] = CATEGORIES.get(category_entry.get("superkingdom_taxid", ""), "uncategorized")
        row["is_phage"] = "true" if entry.get("is_phage") == "1" else "false"


def add_species_membership(rows):
    """Group species, including NT-only/NR-only taxa, into genus rows."""
    children = {}
    for (level, tax_id), species in list(rows.items()):
        if level != 1:
            continue
        genus_id = species["genus_tax_id"]
        if (2, genus_id) not in rows:
            # ReportsHelper.cleanup_missing_genus_counts creates a singleton
            # bucket so species with missing genus records remain in the report.
            genus_id = -1900000000 - tax_id
            genus = dict(species, tax_id=genus_id, tax_level=2, genus_tax_id=genus_id)
            genus["name"] = "Ad hoc bucket for {}".format(species.get("name", str(tax_id)).lower())
            rows[(2, genus_id)] = genus
            species["genus_tax_id"] = genus_id
        children.setdefault(genus_id, []).append(tax_id)
    for genus_id, species_ids in children.items():
        # Platform sorts by aggregate score; until scores are available use IDs.
        rows[(2, genus_id)]["species_tax_ids"] = json.dumps(sorted(species_ids))


def normalization_denominator(technology, total_count, before_subsampling, after_subsampling, ercc_counts=None):
    """Adapt PipelineRun.rpm/bpm and MngsReadsStatsLoadService's count inputs."""
    def read_count(path, keys):
        with open(path) as source:
            counts = json.load(source)
        present = [key for key in keys if key in counts]
        if len(present) != 1:
            raise ValueError("{} must contain exactly one of {}".format(path, keys))
        value = counts[present[0]]
        if isinstance(value, bool) or str(value) != str(int(value)) or int(value) < 0:
            raise ValueError("Invalid nonnegative integer count in {}".format(path))
        return int(value)

    if technology == "short-read":
        total = read_count(total_count, ("fastqs",))
        before = read_count(before_subsampling, ("hisat2_host_filtered_out", "hisat2_human_filtered_out"))
        after = read_count(after_subsampling, ("subsampled_out",))
        # Bowtie2's headerless TSV counts individual reads, including both mates.
        ercc = 0
        with open(ercc_counts) as source:
            for line in source:
                if not line.strip():
                    continue
                target, count = line.split()
                count = int(count)
                if not target.startswith("ERCC-") or count < 0:
                    raise ValueError("Invalid Bowtie2 ERCC count: {}".format(line.strip()))
                ercc += count
        total -= ercc
    else:
        total = read_count(total_count, ("original_bases",))
        before = read_count(before_subsampling, ("human_filtered_bases",))
        after = read_count(after_subsampling, ("subsampled_bases",))
    if after > before:
        raise ValueError("Count after subsampling exceeds count before subsampling")
    fraction = float(after) / before if before > 0 else 1.0
    denominator = total * fraction
    if not math.isfinite(denominator) or denominator <= 0:
        raise ValueError("Normalization requires a positive adjusted sample total")
    return denominator


def write_report(counts_json, output_csv, contig_summary_json=None, contig_stats_json=None, contigs_fasta=None,
                 lineage_csv=None, lineage_version=None, technology="short-read", contig_base_counts_json=None,
                 total_count=None, before_subsampling=None, after_subsampling=None, ercc_counts=None):
    """
    Pivot NT/NR species and genus records.

    Preserve the pipeline's count and decimal alignment metrics. E-values in
    the input are log10 values, exported as strings just as in report_csv.
    Family and merged_NT_NR records are outside this report's schema.
    """
    if technology not in ("short-read", "long-read"):
        raise ValueError("Unsupported report technology: {}".format(technology))
    columns = SHORT_READ_COLUMNS if technology == "short-read" else LONG_READ_COLUMNS
    normalization_inputs = (total_count, before_subsampling, after_subsampling)
    if technology == "short-read":
        normalization_inputs += (ercc_counts,)
    elif ercc_counts is not None:
        raise ValueError("ERCC counts are only used by short-read reports")
    denominator = None
    if any(value is not None for value in normalization_inputs):
        if not all(value is not None for value in normalization_inputs):
            raise ValueError("Supply all normalization inputs together")
        denominator = normalization_denominator(technology, total_count, before_subsampling,
                                                after_subsampling, ercc_counts)
    contig_inputs = (contig_summary_json, contig_stats_json, contigs_fasta)
    if technology == "long-read":
        contig_inputs += (contig_base_counts_json,)
    elif contig_base_counts_json is not None:
        raise ValueError("Contig base counts are only used by long-read reports")
    if bool(lineage_csv) != bool(lineage_version):
        raise ValueError("Taxonomy CSV and lineage version must be supplied together")
    if any(value is not None for value in contig_inputs) and not all(value is not None for value in contig_inputs):
        raise ValueError("Supply all contig inputs together (including base counts for long reads)")
    with open(counts_json) as source:
        counts = json.load(source)["pipeline_output"]["taxon_counts_attributes"]

    rows = {}
    seen = set()
    for count in counts:
        count_type = count["count_type"]
        level = int(count["tax_level"])
        if count_type not in ("NT", "NR") or level not in (1, 2):
            continue
        tax_id = int(count["tax_id"])
        genus = int(count["genus_taxid"])
        key = (level, tax_id)
        if (key, count_type) in seen:
            raise ValueError("Duplicate {} record for taxon {} at level {}".format(count_type, tax_id, level))
        seen.add((key, count_type))
        row = rows.setdefault(key, {"tax_id": tax_id, "tax_level": level, "genus_tax_id": genus})
        if row["genus_tax_id"] != genus:
            raise ValueError("Conflicting genus IDs for taxon {}".format(tax_id))
        prefix = count_type.lower()
        for field in ("count", "percent_identity", "alignment_length"):
            row[prefix + "_" + field] = count[field]
        if technology == "long-read":
            row[prefix + "_base_count"] = count["base_count"]
        if denominator is not None:
            field = "base_count" if technology == "long-read" else "count"
            metric = "bpm" if technology == "long-read" else "rpm"
            row[prefix + "_" + metric] = count[field] / denominator * 1000000.0
        e_value = count["e_value"]
        row[prefix + "_e_value"] = "10^{}".format(e_value) if e_value is not None else None

    if contig_summary_json is not None:
        add_contig_metrics(rows, contig_summary_json, contig_stats_json, contigs_fasta,
                           technology, contig_base_counts_json)

    if lineage_csv:
        add_taxon_information(rows, lineage_csv, lineage_version)
    add_species_membership(rows)

    with open(output_csv, "w", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=columns)
        writer.writeheader()
        # Genus followed by its species; score-based platform sorting needs enrichment.
        for row in sorted(rows.values(), key=lambda r: (r["genus_tax_id"], -r["tax_level"], r["tax_id"])):
            writer.writerow(row)


def main(technology):
    parser = argparse.ArgumentParser(description="Generate a CZ ID {} CSV report".format(technology))
    parser.add_argument("counts_json")
    parser.add_argument("output_csv")
    parser.add_argument("--total-count", help="fastqs.count (short reads) or original_bases.count (long reads)")
    parser.add_argument("--before-subsampling", help="Final HISAT2 read count or human_filtered_bases.count")
    parser.add_argument("--after-subsampling", help="subsampled_out.count or subsampled_bases.count")
    if technology == "short-read":
        parser.add_argument("--ercc-counts", help="Headerless bowtie2_ERCC_counts.tsv (empty means zero ERCC reads)")
    parser.add_argument("--contig-summary", help="assembly/combined_contig_summary.json")
    parser.add_argument("--contig-stats", help="assembly/contig_stats.json")
    parser.add_argument("--contigs-fasta", help="assembly/contigs.fasta")
    if technology == "long-read":
        parser.add_argument("--contig-base-counts", help="contig_base_counts.json")
    parser.add_argument("--lineage-csv", help="CZ ID versioned-taxid-lineages.csv.gz")
    parser.add_argument("--lineage-version", help="Taxonomy snapshot date YYYY-MM-DD, matching the alignment index")
    args = parser.parse_args()
    write_report(args.counts_json, args.output_csv, args.contig_summary, args.contig_stats, args.contigs_fasta,
                 args.lineage_csv, args.lineage_version, technology, getattr(args, "contig_base_counts", None),
                 args.total_count, args.before_subsampling, args.after_subsampling, getattr(args, "ercc_counts", None))
