#!/usr/bin/env python3
"""Check that non-legacy ROOT analysis macros bind only branches present in the booked schema."""
from pathlib import Path
import re, sys
ROOT=Path(__file__).resolve().parents[1]
run=(ROOT/'src/RunAction.cc').read_text()
# Map ntuple name -> branch names from booking blocks.
schema={}
for m in re.finditer(r'am->CreateNtuple\("([^"]+)".*?am->FinishNtuple\(k\w+\);',run,re.S):
    schema[m.group(1)]={x.group(1) for x in re.finditer(r'am->CreateNtuple[IFDS]Column\("([^"]+)"\)',m.group(0))}

bad=[]
for path in sorted((ROOT/'analysis').glob('*.C')):
    if 'legacy' in path.name: continue
    text=path.read_text()
    # variable -> tree name via f.Get("Tree") or f->Get("Tree")
    var_to_tree={}
    for m in re.finditer(r'auto\*\s+(\w+)\s*=\s*dynamic_cast<TTree\*>\(f(?:->|\.)Get\("([^"]+)"\)\)',text):
        var_to_tree[m.group(1)]=m.group(2)
    for var, branch in re.findall(r'(\w+)->SetBranchAddress\("([^"]+)"',text):
        tree=var_to_tree.get(var)
        if tree and branch not in schema.get(tree,set()): bad.append(f'{path.name}: {tree}.{branch} not in schema')
if bad:
    print('\n'.join(bad)); sys.exit(1)
print('Analysis branch-vs-schema checks: OK')
