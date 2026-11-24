#!/usr/bin/env python3
"""Deterministic cache-busting for the local CTAS Explorer.

Computes a content-based build tag from every served asset and data file (assets/*,
data/*.json; sorted by path, each prefixed by its relative path), then stamps `?v=<build>`
onto the CSS/JS asset references in every Explorer HTML page. The JS (assets/app.js) reads
that tag from its own <script src> query and appends it to every JSON fetch (see loadJSON),
so a content change -> new tag -> fresh JS and JSON.

Run as the final Explorer build step (after scripts/build_explorer_data.py). Reproducible
and idempotent: re-running with unchanged content yields the same tag. The HTML pages and
data/build_version.txt are excluded from the hash because they carry the tag.
"""
import hashlib, os, re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # staged site root
HASHED = sorted(
    [os.path.join("assets", f) for f in os.listdir(os.path.join(ROOT, "assets"))
     if f.endswith((".js", ".css"))] +
    [os.path.join("data", f) for f in os.listdir(os.path.join(ROOT, "data")) if f.endswith(".json")])

h = hashlib.sha256()
for p in HASHED:
    h.update(p.encode() + b"\0")
    h.update(open(os.path.join(ROOT, p), "rb").read())
build = h.hexdigest()[:10]

HTML = ["index.html", "loci.html", "method.html"]
for p in HTML:
    fp = os.path.join(ROOT, p)
    s = open(fp, encoding="utf-8").read()
    s = re.sub(r'(href=")(assets/styles\.css)(?:\?v=[0-9a-f]+)?(")', r'\1\2?v=' + build + r'\3', s)
    s = re.sub(r'(src=")(assets/app\.js)(?:\?v=[0-9a-f]+)?(")', r'\1\2?v=' + build + r'\3', s)
    open(fp, "w", encoding="utf-8").write(s)

open(os.path.join(ROOT, "data", "build_version.txt"), "w").write(build + "\n")
print("cache build tag:", build)
print("hashed:", ", ".join(HASHED))
print("stamped:", ", ".join(HTML))
