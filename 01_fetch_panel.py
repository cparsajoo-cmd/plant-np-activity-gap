"""Retrieve, for each target in the panel, the set of ChEMBL compounds that
carry at least one activity value with an associated pChEMBL value.

Writes data/panel_compounds.json  ->  {target_chembl_id: [molecule_chembl_id, ...]}

Usage:  python3 01_fetch_panel.py
"""
import json, os, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
BASE = "https://www.ebi.ac.uk/chembl/api/data"

PANEL = [
    ("Acetylcholinesterase (human)", "CHEMBL220"),
    ("Butyrylcholinesterase (human)", "CHEMBL1914"),
    ("Monoamine oxidase A (human)", "CHEMBL1951"),
    ("5-Lipoxygenase (human)", "CHEMBL215"),
    ("Neutrophil elastase (human)", "CHEMBL248"),
    ("Cyclooxygenase-2 (human)", "CHEMBL230"),
    ("Xanthine oxidase (bovine)", "CHEMBL3649"),
    ("Tyrosinase (mushroom)", "CHEMBL3318"),
    ("Aldose reductase (human)", "CHEMBL1900"),
    ("Myeloperoxidase (human)", "CHEMBL2439"),
]


def get(url):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return json.load(r)
        except Exception:
            if attempt == 4:
                raise
            time.sleep(3 * (attempt + 1))


def compounds(target_id):
    out = set()
    url = BASE + "/activity.json?" + urllib.parse.urlencode({
        "target_chembl_id": target_id,
        "pchembl_value__isnull": "false",
        "only": "molecule_chembl_id",
        "limit": 1000,
    })
    while url:
        page = get(url)
        for a in page["activities"]:
            if a.get("molecule_chembl_id"):
                out.add(a["molecule_chembl_id"])
        nxt = page["page_meta"].get("next")
        url = ("https://www.ebi.ac.uk" + nxt) if nxt else None
    return out


os.makedirs(DATA, exist_ok=True)
result = {}
for name, tid in PANEL:
    s = compounds(tid)
    result[tid] = sorted(s)
    print("%-32s %-12s %6d compounds" % (name, tid, len(s)), flush=True)

union = set().union(*[set(v) for v in result.values()])
print("\nunion across the panel: %d distinct compounds" % len(union))
with open(os.path.join(DATA, "panel_compounds.json"), "w") as f:
    json.dump(result, f)
print("written: data/panel_compounds.json")
