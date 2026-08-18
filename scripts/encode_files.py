#!/usr/bin/env python3
"""List released files for an ENCODE experiment, optionally filtered.

Kept as a separate script rather than a heredoc: `curl ... | python3 - <<EOF`
cannot work, because the heredoc and the pipe both claim stdin and the heredoc
wins, so the piped JSON is never read.
"""
import argparse
import json
import sys
import urllib.request

BASE = "https://www.encodeproject.org"


def fetch(accession: str) -> dict:
    url = f"{BASE}/experiments/{accession}/?format=json"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        body = r.read().decode()
    try:
        return json.loads(body)
    except json.JSONDecodeError:
        sys.exit(f"ENCODE returned non-JSON for {accession} "
                 f"(first 200 chars): {body[:200]!r}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("accession")
    ap.add_argument("--file-format")
    ap.add_argument("--assembly")
    ap.add_argument("--output-type")
    ap.add_argument("--url-only", action="store_true",
                    help="print only the href of the first match")
    a = ap.parse_args()

    rows = []
    for f in fetch(a.accession).get("files", []):
        if f.get("status") != "released":
            continue
        if a.file_format and f.get("file_format") != a.file_format:
            continue
        if a.assembly and f.get("assembly") != a.assembly:
            continue
        if a.output_type and f.get("output_type") != a.output_type:
            continue
        rows.append((f["accession"], f.get("output_type", ""),
                     f.get("assembly", ""),
                     round(f.get("file_size", 0) / 1e9, 2),
                     BASE + f["href"]))

    if not rows:
        sys.exit("No matching files. Re-run without filters to see what exists.")

    if a.url_only:
        print(rows[0][4])
    else:
        print("accession\toutput_type\tassembly\tsize_GB\turl")
        for r in rows:
            print("\t".join(str(x) for x in r))


if __name__ == "__main__":
    main()
