"""Sanity check: draw a random sample of the ChEMBL entries that the mapping
flagged as plant natural products with panel activity, and print their names.

Every compound printed should be a recognisable plant natural product. If the
list contains synthetic drugs, the structure matching is wrong.

Usage:  python3 03_validate.py [sample_size]
"""
import json, os, random, sys, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20
SEED = 7

path = os.path.join(DATA, "active_plant_ids.json")
if not os.path.exists(path):
    sys.exit("missing data/active_plant_ids.json - run 02_map_lotus_chembl.py first")

ids = json.load(open(path))
print("flagged ChEMBL entries: %d" % len(ids))
random.seed(SEED)
sample = random.sample(ids, min(N, len(ids)))

def get(url):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception:
            if attempt == 3:
                raise
            time.sleep(3)

named = 0
for i in range(0, len(sample), 25):
    chunk = sample[i:i + 25]
    url = "https://www.ebi.ac.uk/chembl/api/data/molecule.json?" + urllib.parse.urlencode({
        "molecule_chembl_id__in": ",".join(chunk),
        "only": "molecule_chembl_id,pref_name",
        "limit": 50,
    })
    for m in get(url)["molecules"]:
        name = m.get("pref_name") or "(no preferred name in ChEMBL)"
        if m.get("pref_name"):
            named += 1
        print("  %-14s %s" % (m["molecule_chembl_id"], name))

print("\n%d of %d sampled entries carry a preferred name." % (named, len(sample)))
print("Entries without one are usually literature compounds that ChEMBL stores")
print("by structure only; that is expected and not an error.")
