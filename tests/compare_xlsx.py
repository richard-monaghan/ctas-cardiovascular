#!/usr/bin/env python3
"""Cell-by-cell comparison of rebuilt workbooks with the released supplementary tables.
Usage: python3 tests/compare_xlsx.py <rebuilt.xlsx> <released.xlsx>   (exit 0 if every sheet and cell is identical)"""
import sys, openpyxl
a = openpyxl.load_workbook(sys.argv[1], read_only=True); b = openpyxl.load_workbook(sys.argv[2], read_only=True)
nd = nc = 0; notes = []
if a.sheetnames != b.sheetnames: notes.append(f'sheet lists differ: {a.sheetnames} vs {b.sheetnames}')
for s in [x for x in b.sheetnames if x in a.sheetnames]:
    ra = list(a[s].iter_rows(values_only=True)); rb = list(b[s].iter_rows(values_only=True))
    if len(ra) != len(rb): notes.append(f'{s}: {len(ra)} vs {len(rb)} rows')
    for x, y in zip(ra, rb):
        for u, v in zip(x, y):
            nc += 1
            if u != v: nd += 1
print(f'{sys.argv[1].split("/")[-1]}: {nc} cells compared, {nd} differ; {"; ".join(notes) or "structure identical"}')
sys.exit(0 if nd == 0 and not notes else 1)
