"""Map the plant natural products catalogued in LOTUS onto the bioactivity
records of ChEMBL, and count what is and is not measured.

Method
------
1. From the LOTUS frozen dataset, keep the structure-organism pairs whose
   organism phylum is Streptophyta (land plants). The distinct InChIKeys of
   those rows are the plant natural-product set.
2. Match those structures to ChEMBL on the connectivity layer of the InChIKey
   (its first 14 characters, i.e. the molecular skeleton without
   stereochemistry). This is deliberately permissive: a permissive match finds
   MORE plant compounds in ChEMBL, so the resulting gap is a lower bound.
   The stricter full-key match is reported alongside for comparison.
3. Intersect with the per-target compound sets from 01_fetch_panel.py.

Writes results.json and prints the same figures.

Usage:  python3 02_map_lotus_chembl.py
"""
import collections, csv, gzip, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
csv.field_size_limit(10 ** 7)

NAMES = {
    "CHEMBL220": "Acetylcholinesterase (human)",
    "CHEMBL1914": "Butyrylcholinesterase (human)",
    "CHEMBL1951": "Monoamine oxidase A (human)",
    "CHEMBL215": "5-Lipoxygenase (human)",
    "CHEMBL248": "Neutrophil elastase (human)",
    "CHEMBL230": "Cyclooxygenase-2 (human)",
    "CHEMBL3649": "Xanthine oxidase (bovine)",
    "CHEMBL3318": "Tyrosinase (mushroom)",
    "CHEMBL1900": "Aldose reductase (human)",
    "CHEMBL2439": "Myeloperoxidase (human)",
}
PHYLUM = "Streptophyta"

for required in ("lotus_meta.csv.gz", "chembl37_chemreps.txt.gz", "panel_compounds.json"):
    if not os.path.exists(os.path.join(DATA, required)):
        sys.exit("missing data/%s - run fetch_data.py and 01_fetch_panel.py first" % required)

# ---- step 1: plant structures from LOTUS -----------------------------------
plant = set()
pairs_seen = 0
with gzip.open(os.path.join(DATA, "lotus_meta.csv.gz"), "rt", encoding="utf-8", errors="replace") as f:
    for row in csv.DictReader(f):
        pairs_seen += 1
        if (row.get("organism_taxonomy_03phylum") or "").strip() != PHYLUM:
            continue
        key = (row.get("structure_inchikey") or "").strip()
        if len(key) >= 14:
            plant.add(key)
blocks = {k[:14] for k in plant}
print("LOTUS rows read:                         %9d" % pairs_seen)
print("plant structures (%s):          %9d" % (PHYLUM, len(plant)))
print("distinct connectivity blocks:            %9d" % len(blocks), flush=True)

# ---- step 2: locate them in ChEMBL -----------------------------------------
block_to_ids = collections.defaultdict(set)
exact = 0
rows = 0
with gzip.open(os.path.join(DATA, "chembl37_chemreps.txt.gz"), "rt", encoding="utf-8", errors="replace") as f:
    f.readline()                                  # header
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 4:
            continue
        rows += 1
        key = p[3].strip()
        if key[:14] in blocks:
            block_to_ids[key[:14]].add(p[0])
            if key in plant:
                exact += 1
print("ChEMBL compounds read:                   %9d" % rows)

matched = {k for k in plant if k[:14] in block_to_ids}
all_ids = set()
for b in block_to_ids.values():
    all_ids |= b
print()
print("plant structures in ChEMBL, exact key:   %9d" % exact)
print("plant structures in ChEMBL, connectivity:%9d  (%.1f%%)" % (len(matched), 100.0 * len(matched) / len(plant)))
print("ChEMBL entries they correspond to:       %9d" % len(all_ids), flush=True)

# ---- step 3: intersect with the panel --------------------------------------
panel = json.load(open(os.path.join(DATA, "panel_compounds.json")))
table, active_ids, pairs_with = [], set(), 0
for tid, ids in panel.items():
    target_set = set(ids)
    entries = target_set & all_ids
    structures = sum(1 for k in matched if block_to_ids[k[:14]] & target_set)
    active_ids |= entries
    pairs_with += structures
    table.append((NAMES.get(tid, tid), tid, len(target_set), len(entries), structures))
table.sort(key=lambda r: -r[2])

print()
print("%-32s %10s %9s %11s" % ("target", "compounds", "plant NPs", "structures"))
for name, tid, total, entries, structures in table:
    print("%-32s %10d %9d %11d" % (name, total, entries, structures))

n_targets = len(panel)
active_structures = sum(1 for k in matched if block_to_ids[k[:14]] & active_ids)
total_pairs = len(plant) * n_targets

print()
print("plant structures with >=1 panel activity:%9d  (%.1f%% of %d)"
      % (active_structures, 100.0 * active_structures / len(plant), len(plant)))
print("ChEMBL entries among those:              %9d  (%.1f%% of %d)"
      % (len(active_ids), 100.0 * len(active_ids) / len(all_ids), len(all_ids)))
print()
print("structure-target pairs defined by panel: %9d" % total_pairs)
print("  with a recorded activity:              %9d" % pairs_with)
print("  with none:                             %9d" % (total_pairs - pairs_with))

json.dump({
    "phylum": PHYLUM,
    "lotus_rows": pairs_seen,
    "plant_structures": len(plant),
    "in_chembl_exact": exact,
    "in_chembl_connectivity": len(matched),
    "chembl_entries": len(all_ids),
    "active_structures": active_structures,
    "active_entries": len(active_ids),
    "pairs_total": total_pairs,
    "pairs_with_activity": pairs_with,
    "per_target": [
        {"target": n, "chembl_id": t, "compounds": c, "plant_np_entries": e, "plant_np_structures": s}
        for n, t, c, e, s in table
    ],
}, open(os.path.join(HERE, "results.json"), "w"), indent=1)
json.dump(sorted(active_ids), open(os.path.join(DATA, "active_plant_ids.json"), "w"))
print("\nwritten: results.json, data/active_plant_ids.json")
