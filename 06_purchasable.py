"""How much of the untested plant chemical space can actually be bought?

Commercial availability is taken from MolPort, a compound-sourcing aggregator
that UniChem indexes, via UniChem's PubChem-to-MolPort whole-source mapping.
LOTUS records a PubChem identifier for most of its structures, so the two can
be joined without any vendor website being queried.

This is a lower bound on availability: MolPort does not index every supplier,
specialist natural-product houses in particular, and a structure that LOTUS
does not resolve to a PubChem identifier cannot be counted at all.

Usage:  python3 06_purchasable.py
"""
import collections, gzip, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
N = os.path.join(HERE, "data")
MAP = os.path.join(N, "src22src28.txt.gz")
for req in ("strep_keys.json", "key2cid.json", "chembl37_chemreps.txt.gz", "panel_compounds.json"):
    if not os.path.exists(os.path.join(N, req)):
        sys.exit("missing data/" + req)
if not os.path.exists(MAP):
    sys.exit("missing data/src22src28.txt.gz (UniChem wholeSourceMapping src_id22)")

plant = set(json.load(open(os.path.join(N, "strep_keys.json"))))
key2cid = json.load(open(os.path.join(N, "key2cid.json")))
key2cid = {k: v for k, v in key2cid.items() if k in plant}
print("plant structures (Streptophyta)      : %7d" % len(plant))
print("  with a PubChem identifier in LOTUS : %7d" % len(key2cid), flush=True)

molport = set()
with gzip.open(MAP, "rt", encoding="utf-8", errors="replace") as f:
    f.readline()
    for line in f:
        cid = line.split("\t", 1)[0].strip()
        if cid:
            molport.add(cid)
print("distinct PubChem ids listed by MolPort: %7d" % len(molport), flush=True)

buyable = {k for k, c in key2cid.items() if c in molport}
print("  plant structures listed by MolPort : %7d  (%.1f%% of those with a CID)"
      % (len(buyable), 100.0 * len(buyable) / len(key2cid)))

# which plant structures already carry a recorded activity on the panel
blocks = {k[:14] for k in plant}
blk2ids = collections.defaultdict(set)
with gzip.open(os.path.join(N, "chembl37_chemreps.txt.gz"), "rt", encoding="utf-8", errors="replace") as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4 and p[3].strip()[:14] in blocks:
            blk2ids[p[3].strip()[:14]].add(p[0])
panel = json.load(open(os.path.join(N, "panel_compounds.json")))
active = set()
for ids in panel.values():
    active |= set(ids)
tested = {k for k in plant if blk2ids.get(k[:14], set()) & active}
print("  plant structures with panel activity: %7d" % len(tested))

pool = buyable - tested
print()
print("ADDRESSABLE POOL (buyable, no recorded panel activity): %d" % len(pool))
print("  that is %.1f%% of the %d plant structures LOTUS resolves to a CID,"
      % (100.0 * len(pool) / len(key2cid), len(key2cid)))
print("  and %d times the 160 compounds the validation design requires." % (len(pool) // 160))
json.dump({"plant": len(plant), "with_cid": len(key2cid), "molport_cids": len(molport),
           "buyable": len(buyable), "tested": len(tested), "pool": len(pool)},
          open(os.path.join(HERE, "purchasable.json"), "w"), indent=1)
print("\nwritten: purchasable.json")
