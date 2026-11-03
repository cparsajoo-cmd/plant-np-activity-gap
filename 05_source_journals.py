"""Which literature does ChEMBL's activity data for the panel come from?

For every target in the panel, every activity carrying a pChEMBL value is
traced back to its source document, and every document to its journal. The
counts are complete, not sampled.

Writes journals.json and prints the table.

Usage:  python3 05_source_journals.py
"""
import collections, json, os, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
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
EXPERIMENTAL = {"CHEMBL220", "CHEMBL1914", "CHEMBL2439"}


def get(url):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return json.load(r)
        except Exception:
            if attempt == 4:
                raise
            time.sleep(3 * (attempt + 1))


def activity_documents(target_id):
    """document_chembl_id -> number of pChEMBL activities on this target."""
    out = collections.Counter()
    url = BASE + "/activity.json?" + urllib.parse.urlencode({
        "target_chembl_id": target_id, "pchembl_value__isnull": "false",
        "only": "document_chembl_id", "limit": 1000})
    while url:
        page = get(url)
        for a in page["activities"]:
            if a.get("document_chembl_id"):
                out[a["document_chembl_id"]] += 1
        nxt = page["page_meta"].get("next")
        url = ("https://www.ebi.ac.uk" + nxt) if nxt else None
    return out


def journals(doc_ids):
    """document_chembl_id -> journal (None where ChEMBL records none)."""
    out = {}
    ids = sorted(doc_ids)
    for i in range(0, len(ids), 25):
        chunk = ids[i:i + 25]
        url = BASE + "/document.json?" + urllib.parse.urlencode({
            "document_chembl_id__in": ",".join(chunk), "only": "document_chembl_id,journal",
            "limit": 50})
        for d in get(url)["documents"]:
            out[d["document_chembl_id"]] = d.get("journal")
        if i and i % 500 == 0:
            print("    %d/%d documents" % (i, len(ids)), flush=True)
    return out


per_target, all_docs = {}, collections.Counter()
for name, tid in PANEL:
    docs = activity_documents(tid)
    per_target[tid] = docs
    all_docs.update(docs)
    print("%-32s %5d activities from %4d documents" % (name, sum(docs.values()), len(docs)), flush=True)

print("\npanel total: %d activities from %d distinct documents" % (sum(all_docs.values()), len(all_docs)))
print("resolving journals ...", flush=True)
jmap = journals(all_docs)
missing = [d for d in all_docs if d not in jmap]
print("documents resolved: %d, unresolved: %d" % (len(jmap), len(missing)))


def tally(docs):
    acts, dcount = collections.Counter(), collections.Counter()
    for d, n in docs.items():
        j = jmap.get(d) or "(no journal recorded)"
        acts[j] += n
        dcount[j] += 1
    return acts, dcount


acts, dcount = tally(all_docs)
total = sum(acts.values())
print("\n%-34s %9s %7s %8s" % ("journal", "activities", "share", "papers"))
print("-" * 62)
for j, n in acts.most_common(20):
    print("%-34s %9d %6.1f%% %8d" % (j[:34], n, 100.0 * n / total, dcount[j]))

exp_docs = collections.Counter()
for tid in EXPERIMENTAL:
    exp_docs.update(per_target[tid])
eacts, edcount = tally(exp_docs)
etotal = sum(eacts.values())
print("\nThe three experimental targets only (%d activities, %d documents):" % (etotal, len(exp_docs)))
print("%-34s %9s %7s" % ("journal", "activities", "share"))
print("-" * 52)
for j, n in eacts.most_common(12):
    print("%-34s %9d %6.1f%%" % (j[:34], n, 100.0 * n / etotal))

json.dump({"panel_activities": total, "panel_documents": len(all_docs),
           "by_journal_activities": dict(acts), "by_journal_papers": dict(dcount),
           "experimental_targets": {"activities": etotal, "documents": len(exp_docs),
                                     "by_journal_activities": dict(eacts)}},
          open(os.path.join(HERE, "journals.json"), "w"), indent=1)
print("\nwritten: journals.json")
