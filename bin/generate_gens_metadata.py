#!/usr/bin/env python3
"""Generate Gens sample metadata and annotation tracks from Nallo results."""

import argparse
from collections import Counter, defaultdict
import math
from pathlib import Path


CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
UPD_LABELS = (
    "Total SNPs", "Non-informative", "Non-informative (%)",
    "Mismatch father", "Mismatch father (%)", "Mismatch mother",
    "Mismatch mother (%)", "Anti-UPD", "Anti-UPD (%)",
)


def canonical(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom


def reference_lengths(path):
    lengths = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split("\t")
        chrom = canonical(fields[0])
        if chrom in CHROMS[:22]:
            length = int(fields[1])
            if length <= 0 or chrom in lengths:
                raise ValueError(f"Invalid or duplicate reference length for {chrom}")
            lengths[chrom] = length
    if not lengths:
        raise ValueError("Reference FAI has no autosomal lengths")
    return lengths


def read_roh(path, sample, min_quality=85):
    intervals = defaultdict(list)
    tracks = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[0] != "RG":
            continue
        if len(fields) < 8:
            raise ValueError(f"Malformed ROH region: {line}")
        if fields[1] != sample or float(fields[7]) < min_quality:
            continue
        chrom, start, end = fields[2], int(fields[3]), int(fields[4])
        if start < 1 or end < start:
            raise ValueError(f"Invalid ROH interval: {line}")
        # bcftools roh RG positions are one-based, inclusive; BED is half-open.
        bed_start = start - 1
        tracks.append([chrom, bed_start, end, "LOH", ".", ".", ".", ".", "rgb(255,186,60)"])
        if canonical(chrom) in CHROMS[:22]:
            intervals[canonical(chrom)].append((bed_start, end))
    return intervals, tracks


def covered_bases(intervals, lengths):
    total = 0
    for chrom, regions in intervals.items():
        if chrom not in lengths:
            raise ValueError(f"Missing reference length for ROH chromosome {chrom}")
        last_end = 0
        for start, end in sorted(regions):
            if end > lengths[chrom]:
                raise ValueError(f"ROH exceeds reference chromosome {chrom}")
            total += max(0, end - max(start, last_end))
            last_end = max(last_end, end)
    return total


def copy_numbers(path, sex):
    ratios = defaultdict(list)
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith("@") or line.startswith("CONTIG"):
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            raise ValueError(f"Malformed copy ratio row: {line}")
        chrom = canonical(fields[0])
        if chrom in CHROMS:
            ratio = float(fields[3])
            if not math.isfinite(ratio):
                raise ValueError(f"Nonfinite copy ratio on {chrom}")
            ratios[chrom].append(ratio)
    result = []
    for chrom in CHROMS:
        if chrom not in ratios or (chrom == "Y" and sex == "F"):
            continue
        baseline = 1 if sex == "M" and chrom in ("X", "Y") else 2
        result.append((chrom, "Estimated copy number", round(baseline * 2 ** (sum(ratios[chrom]) / len(ratios[chrom])), 2)))
    return result


def upd_summary(path):
    counts = defaultdict(Counter)
    for line in Path(path).read_text().splitlines():
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            raise ValueError(f"Malformed UPD site: {line}")
        chrom = canonical(fields[0])
        if chrom in CHROMS:
            counts[chrom][fields[3]] += 1
    result = []
    for chrom in CHROMS:
        categories = counts[chrom]
        total = sum(categories.values())
        paternal = categories["UPD_PATERNAL_ORIGIN"]
        maternal = categories["UPD_MATERNAL_ORIGIN"]
        anti = categories["ANTI_UPD"]
        # TODO: Review the inherited OLWGS definition. These three categories
        # are informative in the UPD caller, despite this Gens label.
        non_informative = paternal + maternal + anti
        values = (total, non_informative, percent(non_informative, total),
                  maternal, percent(maternal, total), paternal, percent(paternal, total),
                  anti, percent(anti, total))
        result.extend((chrom, label, value) for label, value in zip(UPD_LABELS, values))
    return result


def percent(part, total):
    return round(100 * part / total, 1) if total else 0


def upd_track(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            raise ValueError(f"Malformed UPD region: {line}")
        details = dict(item.split("=", 1) for item in fields[3].split(";"))
        origin = details.get("ORIGIN")
        if origin not in ("MATERNAL", "PATERNAL"):
            raise ValueError(f"Unknown UPD origin: {origin}")
        # The UPD caller emits BED regions, so coordinates are copied unchanged.
        rows.append([fields[0], fields[1], fields[2],
                     f"Uniparental segment ({origin})", ".", ".", ".", ".",
                     "rgb(255,75,75)" if origin == "MATERNAL" else "rgb(75,75,255)"])
    return rows


def write_tsv(path, header, rows):
    with Path(path).open("w") as handle:
        if header:
            handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--sex", choices=("M", "F"), required=True)
    parser.add_argument("--fai", type=Path, required=True)
    parser.add_argument("--roh", type=Path, required=True)
    parser.add_argument("--ratios", type=Path, required=True)
    parser.add_argument("--upd-regions", type=Path)
    parser.add_argument("--upd-sites", type=Path)
    args = parser.parse_args()
    if bool(args.upd_regions) != bool(args.upd_sites):
        parser.error("UPD regions and sites must be supplied together")
    lengths = reference_lengths(args.fai)
    intervals, roh_rows = read_roh(args.roh, args.sample)
    loh_percent = 100 * covered_bases(intervals, lengths) / sum(lengths.values())
    chrom_rows = copy_numbers(args.ratios, args.sex)
    if args.upd_sites:
        chrom_rows.extend(upd_summary(args.upd_sites))
    prefix = args.sample
    write_tsv(f"{prefix}.meta.tsv", ("type", "value"), [("%Autosomal LOH", loh_percent)])
    write_tsv(f"{prefix}.chrom_meta.tsv", ("Chromosome", "type", "value"), chrom_rows)
    write_tsv(f"{prefix}.gens_track.roh.bed", None, roh_rows)
    if args.upd_regions:
        write_tsv(f"{prefix}.gens_track.upd.bed", None, upd_track(args.upd_regions))


if __name__ == "__main__":
    main()
