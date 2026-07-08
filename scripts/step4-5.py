import csv
from collections import defaultdict
from datetime import datetime

# =====================
# CONFIG
# =====================

NODES_FILE = "data/processed/canonical_nodes.csv"
EDGES_FILE = "data/processed/edges_clean.csv"
REPORT_FILE = "data/processed/step_4_5_validation_report.txt"

MAX_SAMPLE = 10

# Relations where self-loops are biologically invalid
INVALID_SELF_LOOP_RELATIONS = {
    "biolink:orthologous_to",
    "biolink:interacts_with",
    "biolink:genetically_interacts_with",
}

# =====================
# LOAD CANONICAL NODES
# =====================

print("Loading canonical nodes...")

canonical_ids = set()

with open(NODES_FILE, "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        canonical_ids.add(row["canonical_id"])

# =====================
# VALIDATION COUNTERS
# =====================

raw_id_leak = 0
raw_id_samples = []

unresolved_ids = 0
unresolved_samples = []

duplicate_edges = 0

invalid_self_loops = 0
self_loop_samples = []

evidence_errors = 0
evidence_samples = []

seen_edges = set()

total_edges = 0

# =====================
# VALIDATE EDGES
# =====================

print("Validating canonical edges...")

with open(EDGES_FILE, "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)

    for row in reader:
        total_edges += 1

        s = row["source_id"]
        t = row["target_id"]
        r = row["relation"]
        w = row["weight"]
        ds = row["dataset_sources"]

        # ---- RAW ID LEAKAGE (TRUE CHECK)
        if s not in canonical_ids:
            raw_id_leak += 1
            if len(raw_id_samples) < MAX_SAMPLE:
                raw_id_samples.append(f"Edge source: {s}")

        if t not in canonical_ids:
            raw_id_leak += 1
            if len(raw_id_samples) < MAX_SAMPLE:
                raw_id_samples.append(f"Edge target: {t}")

        # ---- UNRESOLVED IDS (NULL / EMPTY)
        if not s or not t:
            unresolved_ids += 1
            if len(unresolved_samples) < MAX_SAMPLE:
                unresolved_samples.append(f"{s} -> {t}")

        # ---- EDGE DEDUPLICATION
        edge_key = (s, r, t)
        if edge_key in seen_edges:
            duplicate_edges += 1
        else:
            seen_edges.add(edge_key)

        # ---- EVIDENCE INTEGRITY
        if not ds or ds.strip() == "":
            evidence_errors += 1
            if len(evidence_samples) < MAX_SAMPLE:
                evidence_samples.append(f"{s} --[{r}]-> {t}")

        # ---- INVALID SELF-LOOPS
        if s == t and r.lower() in INVALID_SELF_LOOP_RELATIONS:
            invalid_self_loops += 1
            if len(self_loop_samples) < MAX_SAMPLE:
                self_loop_samples.append(f"{s} --[{r}]-> {t}")

# =====================
# REPORT GENERATION
# =====================

def status(count):
    return "PASS" if count == 0 else "FAIL"

overall_fatal = raw_id_leak > 0
overall_error = invalid_self_loops > 0

with open(REPORT_FILE, "w", encoding="utf-8") as out:
    out.write("=" * 80 + "\n")
    out.write("ENTITY RESOLUTION VALIDATION REPORT\n")
    out.write("=" * 80 + "\n")
    out.write(f"Generated: {datetime.now()}\n")
    out.write(f"Overall Status: {'FATAL' if overall_fatal else 'PASS'}\n\n")

    out.write("-" * 80 + "\n")
    out.write("SUMMARY STATISTICS\n")
    out.write("-" * 80 + "\n")
    out.write(f"Total Edges: {total_edges}\n")
    out.write(f"Total Nodes: {len(canonical_ids)}\n\n")

    out.write("-" * 80 + "\n")
    out.write("DETAILED CHECKS\n")
    out.write("-" * 80 + "\n\n")

    # Raw ID Leakage
    out.write("❌ Raw ID Leakage Check\n")
    out.write(f"   Status: {status(raw_id_leak)}\n")
    out.write(f"   Violations: {raw_id_leak}\n")
    if raw_id_samples:
        out.write("   Sample:\n")
        for s in raw_id_samples:
            out.write(f"     - {s}\n")
    out.write("\n")

    # Unresolved IDs
    out.write("✅ Unresolved ID Check\n")
    out.write(f"   Status: {status(unresolved_ids)}\n\n")

    # Deduplication
    out.write("✅ Edge Deduplication\n")
    out.write(f"   Status: {status(duplicate_edges)}\n\n")

    # Evidence
    out.write("✅ Evidence Integrity\n")
    out.write(f"   Status: {status(evidence_errors)}\n\n")

    # Self-loops
    out.write("❌ Invalid Self-Loops\n")
    out.write(f"   Status: {status(invalid_self_loops)}\n")
    out.write(f"   Violations: {invalid_self_loops}\n")
    if self_loop_samples:
        out.write("   Sample:\n")
        for s in self_loop_samples:
            out.write(f"     - {s}\n")

print("Validation complete.")
print(f"Report written to: {REPORT_FILE}")
