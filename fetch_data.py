"""Download the two public datasets the analysis uses, into ./data/.

Nothing here is specific to the project: both files are public releases with
permanent identifiers.

  LOTUS frozen dataset, 13 April 2026   doi:10.5281/zenodo.19360665
  ChEMBL 37 chemical representations    EBI FTP, chembl_37_chemreps.txt.gz

Usage:  python3 fetch_data.py
"""
import os, sys, urllib.request, time

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

FILES = [
    ("lotus_meta.csv.gz",
     "https://zenodo.org/api/records/19360665/files/260413_frozen_metadata.csv.gz/content",
     90298678),
    ("chembl37_chemreps.txt.gz",
     "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_chemreps.txt.gz",
     292540141),
]


def fetch(path, url, expected):
    done = os.path.getsize(path) if os.path.exists(path) else 0
    if done >= expected:
        print("  already complete:", os.path.basename(path))
        return
    for attempt in range(1, 21):
        done = os.path.getsize(path) if os.path.exists(path) else 0
        if done >= expected:
            print("  complete:", os.path.basename(path))
            return
        req = urllib.request.Request(url, headers={"Range": "bytes=%d-" % done} if done else {})
        try:
            with urllib.request.urlopen(req, timeout=300) as r, open(path, "ab" if done else "wb") as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
        except Exception as e:
            print("  attempt %d interrupted (%s), resuming" % (attempt, type(e).__name__))
            time.sleep(3)
    size = os.path.getsize(path) if os.path.exists(path) else 0
    if size < expected:
        sys.exit("could not download %s (%d of %d bytes)" % (path, size, expected))


os.makedirs(DATA, exist_ok=True)
for name, url, size in FILES:
    print("fetching", name, "(%.0f MB)" % (size / 1e6))
    fetch(os.path.join(DATA, name), url, size)
print("done")
