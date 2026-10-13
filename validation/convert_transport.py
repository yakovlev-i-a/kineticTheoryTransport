#!/usr/bin/env python3
"""Convert GRI-mech sutherland transport entries to standardKinetic entries."""
import re, sys

TRANS = "/media/yakovlev-i-a/adata-legend-512/openfoam-libs/standardKinetic/of-11/run/counterFlowFlame2D_wedge/chemkin/trans.dat"
THERMO = "constant/thermo"

# Parse CHEMKIN transport.dat -> name -> [geometry, welldepth, diameter, dipole, polar, relaxation]
ld = {}
with open(TRANS) as f:
    for line in f:
        line = line.split("!")[0].strip()
        if not line:
            continue
        tok = line.split()
        ld[tok[0]] = list(map(float, tok[1:7]))

text = open(THERMO).read()

# Find the end of the species list block ("species ... N (...) ;")
sm = re.search(r'species\s*\n\s*\d+\s*\n\((.*?)\)\s*;', text, re.S)
head = text[:sm.end()]
body = text[sm.end():]

# Scan body for top-level "NAME { ... }" blocks in file order.
trans_pat = re.compile(r'transport\s*\{\s*As\s+[^;]+;\s*Ts\s+[^;]+;\s*\}', re.S)
block_pat = re.compile(r'^([^\s{}]+)\s*\{', re.M)

out_body = body
offset = 0
n_converted = 0
for m in block_pat.finditer(body):
    spec = m.group(1)
    # brace-match this block
    depth = 0
    start = m.start()
    for j in range(start, len(body)):
        c = body[j]
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                end = j + 1
                break
    block = body[start:end]

    if spec in ld and trans_pat.search(block):
        g, w, d, dip, pol, rel = ld[spec]
        new_block, n = trans_pat.subn(
            "    transport\n"
            "    {\n"
            f"        geometry        {int(g)};\n"
            f"        welldepth       {w};\n"
            f"        diameter        {d};\n"
            f"        dipole          {dip};\n"
            f"        polar           {pol};\n"
            f"        relaxation      {rel};\n"
            "    }",
            block,
        )
        out_body = out_body[:start + offset] + new_block + out_body[end + offset:]
        offset += len(new_block) - len(block)
        n_converted += 1
    elif spec in ld:
        print(f"SKIP (no sutherland transport): {spec}", file=sys.stderr)

open(THERMO, "w").write(head + out_body)
print(f"converted {n_converted} species")
