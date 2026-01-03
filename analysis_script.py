#!/usr/bin/env python3
from collections import defaultdict, Counter
import re
import csv
import math

FASTA_IN = "results/inutil_aln.fasta"

OUT_BY_PATIENT = "mutations_by_patient.tsv"
OUT_RECURRENT = "recurrent_sites.tsv"
OUT_PATIENT_SUMMARY = "patient_mutation_summary.tsv"
OUT_SUMMARY_STATS = "summary_stats.tsv"  # per-patient clinical summary table

header_re = re.compile(
    r"^(Pt(?P<pt>\d+))_day(?P<day>\d+)_CD4/(?P<cd4>\d+)_VL/(?P<vl>\d+)"
)

def read_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            seqs[name] = "".join(chunks)
    return seqs

def parse_header(h: str) -> dict:
    m = header_re.match(h)
    if not m:
        raise ValueError(f"Header not in expected format: {h}")
    return {
        "pt": int(m.group("pt")),
        "day": int(m.group("day")),
        "cd4": int(m.group("cd4")),
        "vl": int(m.group("vl")),
        "id": m.group(1),
    }

def log10_safe(x: int | float) -> float | str:
    """Return log10(x) if x>0 else empty string."""
    try:
        if x is None or x <= 0:
            return ""
        return math.log10(x)
    except Exception:
        return ""

def response_category(dlog_vl, dcd4) -> str:
    """
    Simple descriptive label. Threshold-free by default:
    - Improved: dlog_vl < 0 and dcd4 > 0
    - Worsened: dlog_vl > 0 and dcd4 < 0
    - Mixed: otherwise (includes any missing values)
    """
    if dlog_vl == "" or dcd4 == "":
        return "Mixed/unknown"
    if dlog_vl < 0 and dcd4 > 0:
        return "Improved (VL↓, CD4↑)"
    if dlog_vl > 0 and dcd4 < 0:
        return "Worsened (VL↑, CD4↓)"
    return "Mixed"

seqs = read_fasta(FASTA_IN)

# group sequences by patient
by_pt = defaultdict(list)
meta = {}
for h, s in seqs.items():
    info = parse_header(h)
    meta[h] = info
    by_pt[info["pt"]].append(h)

# choose day0 and the latest available day
pairs = {}
for pt, hs in by_pt.items():
    day0 = [h for h in hs if meta[h]["day"] == 0]
    if not day0:
        continue
    day0 = day0[0]

    later = sorted([h for h in hs if meta[h]["day"] != 0], key=lambda x: meta[x]["day"])
    if not later:
        continue
    day_last = later[-1]
    pairs[pt] = (day0, day_last)

mut_rows = []
site_counter = Counter()
patient_summary = []
summary_stats_rows = []

for pt, (h0, h1) in sorted(pairs.items()):
    s0 = seqs[h0]
    s1 = seqs[h1]
    if len(s0) != len(s1):
        raise ValueError(f"Alignment lengths differ for patient {pt}")

    # mutation calling
    n_changes = 0
    for i, (a, b) in enumerate(zip(s0, s1), start=1):  # 1-based alignment position
        if a in "-X" or b in "-X":
            continue
        if a != b:
            n_changes += 1
            mut_rows.append({
                "pt": f"Pt{pt:02d}",
                "day0_header": h0,
                "dayN_header": h1,
                "dayN": meta[h1]["day"],
                "CD4_day0": meta[h0]["cd4"],
                "CD4_dayN": meta[h1]["cd4"],
                "VL_day0": meta[h0]["vl"],
                "VL_dayN": meta[h1]["vl"],
                "aln_pos": i,
                "ref_AA": a,
                "alt_AA": b,
            })
            site_counter[i] += 1

    # existing patient mutation summary
    patient_summary.append({
        "pt": f"Pt{pt:02d}",
        "dayN": meta[h1]["day"],
        "n_AA_changes": n_changes,
        "CD4_day0": meta[h0]["cd4"],
        "CD4_dayN": meta[h1]["cd4"],
        "VL_day0": meta[h0]["vl"],
        "VL_dayN": meta[h1]["vl"],
        "VL_fold_change": (meta[h1]["vl"] / meta[h0]["vl"]) if meta[h0]["vl"] else "",
    })

    # NEW: per-patient clinical summary table ("Summary stats")
    cd4_0 = meta[h0]["cd4"]
    cd4_n = meta[h1]["cd4"]
    vl_0 = meta[h0]["vl"]
    vl_n = meta[h1]["vl"]

    cd4_fold = (cd4_n / cd4_0) if cd4_0 else ""
    vl_fold = (vl_n / vl_0) if vl_0 else ""

    summary_stats_rows.append({
        "pt": f"Pt{pt:02d}",
        "CD4_day0": cd4_0,
        "CD4_dayN": cd4_n,
        "CD4_fold_change": cd4_fold,
        "VL_day0": vl_0,
        "VL_dayN": vl_n,
        "VL_fold_change": vl_fold,
        "n_AA_changes": n_changes,
    })

# write outputs
with open(OUT_BY_PATIENT, "w", newline="", encoding="utf-8") as f:
    fieldnames = list(mut_rows[0].keys()) if mut_rows else [
        "pt","day0_header","dayN_header","dayN","CD4_day0","CD4_dayN","VL_day0","VL_dayN",
        "aln_pos","ref_AA","alt_AA"
    ]
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    for r in mut_rows:
        w.writerow(r)

with open(OUT_PATIENT_SUMMARY, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(patient_summary[0].keys()) if patient_summary else [
        "pt","dayN","n_AA_changes","CD4_day0","CD4_dayN","VL_day0","VL_dayN","VL_fold_change"
    ], delimiter="\t")
    w.writeheader()
    for r in patient_summary:
        w.writerow(r)

recurrent = [{"aln_pos": pos, "n_patients_mutated_here": n} for pos, n in site_counter.items() if n >= 2]
recurrent.sort(key=lambda x: (-x["n_patients_mutated_here"], x["aln_pos"]))

with open(OUT_RECURRENT, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["aln_pos","n_patients_mutated_here"], delimiter="\t")
    w.writeheader()
    for r in recurrent:
        w.writerow(r)

# NEW: write "Summary stats"
with open(OUT_SUMMARY_STATS, "w", newline="", encoding="utf-8") as f:
    fieldnames = [
        "pt","CD4_day0","CD4_dayN","CD4_fold_change",
        "VL_day0","VL_dayN","VL_fold_change","n_AA_changes"
    ]
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    for r in summary_stats_rows:
        w.writerow(r)

print(f"Wrote: {OUT_BY_PATIENT}, {OUT_PATIENT_SUMMARY}, {OUT_RECURRENT}, {OUT_SUMMARY_STATS}")
print(f"Patients processed: {len(pairs)}")
print(f"Total AA substitutions (all patients): {len(mut_rows)}")
print(f"Recurrent sites (>=2 patients): {len(recurrent)}")
