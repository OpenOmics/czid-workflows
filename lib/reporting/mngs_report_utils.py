"""Standard-library reporting shared by the short- and long-read workflows."""

import argparse
import csv
import json
import math
from pathlib import Path
import io
import subprocess
import tempfile


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
    "species_tax_ids", "known_pathogen", "lineage_7rank"
]

LONG_READ_COLUMNS = [
    "tax_id", "tax_level", "genus_tax_id", "name", "common_name", "category", "is_phage",
    "nt_bpm", "nt_base_count", "nt_count", "nt_contigs", "nt_contig_b",
    "nt_percent_identity", "nt_alignment_length", "nt_e_value",
    "nr_bpm", "nr_base_count", "nr_count", "nr_contigs", "nr_contig_b",
    "nr_percent_identity", "nr_alignment_length", "nr_e_value", 
    "species_tax_ids", "known_pathogen", "lineage_7rank"
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
        contribution = base_counts.get(name, 0) if base_counts is not None else reads
        for level, tax_id in lineage.items():
            row = rows.get((level, tax_id))
            if row is None or prefix + "_count" not in row:
                continue
            row[prefix + "_contigs"] = row.get(prefix + "_contigs", 0) + 1
            row[prefix + "_" + metric] = row.get(prefix + "_" + metric, 0) + contribution


INVALID_CALL_BASE_ID = -100000000


LINEAGE_FORMAT = "{domain|acellular root|superkingdom};{phylum};{class};{order};{family};{genus};{species}"
RANK_SEPARATOR = "\x1f"
CATEGORIES = ((12884, "viroids"), (2, "bacteria"), (2157, "archaea"),
              (2759, "eukaryota"), (10239, "viruses"))


def taxonkit_annotations(wanted, taxonomy_files, taxonomy_date, taxonkit, threads):
    """Query report IDs with TaxonKit."""
    def run(arguments, input_text=None):
        result = subprocess.run([taxonkit] + arguments, input=input_text,
                                universal_newlines=True, stdout=subprocess.PIPE, check=True)
        return result.stdout

    version = run(["version"]).strip()
    with tempfile.TemporaryDirectory(prefix="report-taxonomy-") as work:
        for name, path in taxonomy_files.items():
            path = Path(path).resolve()
            (Path(work) / name).symlink_to(path)

        # TaxonKit v0.20.0 follows a single merge hop. Resolve chains explicitly
        # before querying it, preserving the mapping back to original report IDs.
        merged = {}
        with open(taxonomy_files["merged.dmp"]) as source:
            for line in source:
                if not line.strip():
                    continue
                fields = [value.strip() for value in line.split("|")]
                old, new = int(fields[0]), int(fields[1])
                merged[old] = new
        resolved = {}
        for taxid in wanted:
            target, seen = taxid, set()
            while target in merged:
                if target in seen:
                    raise ValueError("Cycle in merged.dmp for tax ID {}".format(taxid))
                seen.add(target)
                target = merged[target]
            resolved[taxid] = target
        targets = set(resolved.values())
        base = ["--data-dir", work, "--threads", str(threads)]
        annotations = {}
        if targets:
            output = run(["lineage"] + base + ["-c", "-n", "-r", "-t"],
                         "\n".join(str(taxid) for taxid in sorted(targets)) + "\n")
            for fields in csv.reader(io.StringIO(output), delimiter="\t", quoting=csv.QUOTE_NONE):
                taxid, status, _lineage, lineage_ids, name, rank = fields
                taxid, status = int(taxid), int(status)
                if status in (-1, 0):
                    annotations[taxid] = {"status": "unknown" if status == -1 else "deleted"}
                    continue
                ancestors = {int(value) for value in lineage_ids.split(";") if value}
                ancestors.add(taxid)
                annotations[taxid] = {
                    "status": "current", "resolved_tax_id": taxid,
                    "name": name, "rank": rank,
                    "category": next((label for ancestor, label in CATEGORIES if ancestor in ancestors), "uncategorized"),
                }
            valid = {taxid for taxid, entry in annotations.items() if entry["status"] == "current"}
            if valid:
                output = run(["reformat2"] + base + ["-I", "1", "-t", "-f", LINEAGE_FORMAT.replace(";", RANK_SEPARATOR)],
                             "\n".join(str(taxid) for taxid in sorted(valid)) + "\n")
                for fields in csv.reader(io.StringIO(output), delimiter="\t", quoting=csv.QUOTE_NONE):
                    taxid, lineage, lineage_ids = fields
                    names, rank_ids = lineage.split(RANK_SEPARATOR), lineage_ids.split(RANK_SEPARATOR)
                    taxid = int(taxid)
                    genus = rank_ids[5]
                    annotations[taxid]["genus_tax_id"] = int(genus) if genus else None
                    annotations[taxid]["lineage_7rank"] = ";".join(
                        name.replace("\\", "\\\\").replace(";", "\\;") for name in names)
        result = {}
        for taxid, target in resolved.items():
            entry = dict(annotations[target])
            if target != taxid:
                entry["status"] = "merged" if entry["status"] == "current" else "unresolved_merge"
                entry["resolved_tax_id"] = target
            result[taxid] = entry
    metadata = {
        "schema_version": 2, "annotation_method": "taxonkit",
        "taxonomy_date": taxonomy_date, "taxonomy_date_source": "user-supplied snapshot label",
        "taxonkit_version": version, "lineage_format": LINEAGE_FORMAT,
        "lineage_escaping": "Backslash escapes literal backslash and semicolon in rank names",
    }
    return metadata, result


def add_taxon_information(rows, taxonomy_files, taxonomy_date, taxonkit="taxonkit", threads=4):
    """Annotate original pipeline IDs directly from a frozen NCBI taxdump."""
    wanted = set()
    for row in rows.values():
        for tax_id in (row["tax_id"], row["genus_tax_id"]):
            if tax_id < INVALID_CALL_BASE_ID:
                tax_id = abs(tax_id) % abs(INVALID_CALL_BASE_ID)
            if tax_id > 0:
                wanted.add(tax_id)
    metadata, annotations = taxonkit_annotations(wanted, taxonomy_files, taxonomy_date, taxonkit, threads)

    audit = {}
    for key, row in rows.items():
        tax_id = row["tax_id"]
        rank = "species" if row["tax_level"] == 1 else "genus"
        entry = annotations.get(tax_id, {})
        valid = entry.get("status") in ("current", "merged")
        row["common_name"] = ""
        row["is_phage"] = ""
        row["known_pathogen"] = ""
        row["lineage_7rank"] = (entry.get("lineage_7rank") or "") if valid else ""
        row["category"] = entry.get("category") or "uncategorized"
        if tax_id < INVALID_CALL_BASE_ID:
            parent_id = row["genus_tax_id"] if rank == "species" else abs(tax_id) % abs(INVALID_CALL_BASE_ID)
            parent_rank = "genus" if rank == "species" else "family"
            if parent_id < INVALID_CALL_BASE_ID:
                # A synthetic genus encodes the family for family-only hits.
                parent_id = abs(parent_id) % abs(INVALID_CALL_BASE_ID)
                parent_rank = "family"
            parent = annotations.get(parent_id, {})
            parent_name = parent.get("name") or "taxon {}".format(parent_id)
            row["name"] = "non-{}-specific reads in {} {}".format(rank, parent_rank, parent_name)
            row["category"] = parent.get("category") or "uncategorized"
        elif tax_id < 0:
            row["name"] = "all artificial constructs" if tax_id == -201 else "unclassified taxon bucket {}".format(tax_id)
        else:
            row["name"] = entry.get("name") or "Unresolved taxon {}".format(tax_id)
        pipeline_genus = annotations.get(row["genus_tax_id"], {})
        resolved_genus = pipeline_genus.get("resolved_tax_id")
        current_genus = entry.get("genus_tax_id")
        audit[key] = {
            "tax_id": tax_id, "tax_level": row["tax_level"],
            "pipeline_genus_tax_id": row["genus_tax_id"],
            "resolved_tax_id": entry.get("resolved_tax_id"),
            "status": "synthetic" if tax_id < 0 else entry.get("status", "unknown"),
            "current_rank": entry.get("rank"), "current_genus_tax_id": current_genus,
            "rank_changed": str(entry.get("rank") != rank).lower() if valid else "",
            "genus_changed": str(current_genus != resolved_genus).lower()
                if valid and current_genus is not None and pipeline_genus.get("status") in ("current", "merged") else "",
        }
    return metadata, audit


def write_taxonomy_audit(rows, output_csv, metadata, audit):
    """Publish annotation provenance without altering the platform-style columns."""
    with open(str(output_csv) + ".taxonomy.json", "w") as target:
        json.dump(metadata, target, indent=2, sort_keys=True)
        target.write("\n")
    fields = ["tax_id", "tax_level", "pipeline_genus_tax_id", "report_genus_tax_id",
              "resolved_tax_id", "status", "current_rank", "current_genus_tax_id",
              "rank_changed", "genus_changed"]
    with open(str(output_csv) + ".taxonomy.csv", "w", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=fields)
        writer.writeheader()
        for key in sorted(rows):
            row = rows[key]
            record = dict(audit.get(key, {"tax_id": row["tax_id"], "tax_level": row["tax_level"],
                                          "status": "synthetic"}))
            record["report_genus_tax_id"] = row["genus_tax_id"]
            writer.writerow(record)


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
            genus["lineage_7rank"] = ""
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
        value = int(counts[present[0]])
        if value < 0:
            raise ValueError("Negative count in {}".format(path))
        return value

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
                _target, count = line.split()
                count = int(count)
                if count < 0:
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
                 taxonomy_files=None, taxonomy_date=None, technology="short-read", contig_base_counts_json=None,
                 total_count=None, before_subsampling=None, after_subsampling=None, ercc_counts=None,
                 taxonkit="taxonkit", taxonkit_threads=4):
    """
    Pivot NT/NR species and genus records.

    Preserve the pipeline's count and decimal alignment metrics. E-values in
    the input are log10 values, exported as strings just as in report_csv.
    Family and merged_NT_NR records are outside this report's schema.
    """
    columns = SHORT_READ_COLUMNS if technology == "short-read" else LONG_READ_COLUMNS
    normalization_inputs = (total_count, before_subsampling, after_subsampling)
    if technology == "short-read":
        normalization_inputs += (ercc_counts,)
    denominator = None
    if any(value is not None for value in normalization_inputs):
        if not all(value is not None for value in normalization_inputs):
            raise ValueError("Supply all normalization inputs together")
        denominator = normalization_denominator(technology, total_count, before_subsampling,
                                                after_subsampling, ercc_counts)
    contig_inputs = (contig_summary_json, contig_stats_json, contigs_fasta)
    if technology == "long-read":
        contig_inputs += (contig_base_counts_json,)
    if bool(taxonomy_files) != bool(taxonomy_date):
        raise ValueError("Taxdump files and snapshot date must be supplied together")
    if taxonomy_files and not all(taxonomy_files.values()):
        raise ValueError("Supply all four taxdump files together")
    if any(value is not None for value in contig_inputs) and not all(value is not None for value in contig_inputs):
        raise ValueError("Supply all contig inputs together (including base counts for long reads)")
    with open(counts_json) as source:
        counts = json.load(source)["pipeline_output"]["taxon_counts_attributes"]

    rows = {}
    normalization_field = "base_count" if technology == "long-read" else "count"
    normalization_metric = "bpm" if technology == "long-read" else "rpm"
    for count in counts:
        count_type = count["count_type"]
        level = int(count["tax_level"])
        if count_type not in ("NT", "NR") or level not in (1, 2):
            continue
        tax_id = int(count["tax_id"])
        genus = int(count["genus_taxid"])
        key = (level, tax_id)
        prefix = count_type.lower()
        row = rows.setdefault(key, {"tax_id": tax_id, "tax_level": level, "genus_tax_id": genus})
        if prefix + "_count" in row:
            raise ValueError("Duplicate {} record for taxon {} at level {}".format(count_type, tax_id, level))
        if row["genus_tax_id"] != genus:
            raise ValueError("Conflicting genus IDs for taxon {}".format(tax_id))
        for field in ("count", "percent_identity", "alignment_length"):
            row[prefix + "_" + field] = count[field]
        if technology == "long-read":
            row[prefix + "_base_count"] = count["base_count"]
        if denominator is not None:
            row[prefix + "_" + normalization_metric] = count[normalization_field] / denominator * 1000000.0
        e_value = count["e_value"]
        row[prefix + "_e_value"] = "10^{}".format(e_value) if e_value is not None else None

    if contig_summary_json is not None:
        add_contig_metrics(rows, contig_summary_json, contig_stats_json, contigs_fasta,
                           technology, contig_base_counts_json)

    if taxonomy_files:
        metadata, audit = add_taxon_information(rows, taxonomy_files, taxonomy_date, taxonkit, taxonkit_threads)
    add_species_membership(rows)
    if taxonomy_files:
        write_taxonomy_audit(rows, output_csv, metadata, audit)

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
    for name in ("nodes", "names", "merged", "delnodes"):
        parser.add_argument("--taxonomy-" + name, help="Frozen NCBI " + name + ".dmp")
    parser.add_argument("--taxonkit", default="taxonkit", help="TaxonKit executable")
    parser.add_argument("--taxonkit-threads", type=int, default=4)
    parser.add_argument("--taxonomy-date", help="Declared annotation snapshot YYYY-MM-DD")
    args = parser.parse_args()
    taxonomy_files = {name + ".dmp": getattr(args, "taxonomy_" + name)
                      for name in ("nodes", "names", "merged", "delnodes")}
    if not any(taxonomy_files.values()):
        taxonomy_files = None
    write_report(
        counts_json=args.counts_json,
        output_csv=args.output_csv,
        contig_summary_json=args.contig_summary,
        contig_stats_json=args.contig_stats,
        contigs_fasta=args.contigs_fasta,
        taxonomy_files=taxonomy_files,
        taxonomy_date=args.taxonomy_date,
        technology=technology,
        contig_base_counts_json=getattr(args, "contig_base_counts", None),
        total_count=args.total_count,
        before_subsampling=args.before_subsampling,
        after_subsampling=args.after_subsampling,
        ercc_counts=getattr(args, "ercc_counts", None),
        taxonkit=args.taxonkit,
        taxonkit_threads=args.taxonkit_threads,
    )
