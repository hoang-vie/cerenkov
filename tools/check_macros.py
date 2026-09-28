#!/usr/bin/env python3
"""Static check that project-specific UI commands used by macros are registered in messenger sources."""
from pathlib import Path
import re, sys
ROOT=Path(__file__).resolve().parents[1]
sources='\n'.join(p.read_text() for p in list((ROOT/'src').glob('*Messenger.cc')))
registered=set(re.findall(r'new\s+G4UIcmd\w*\("(/[^"]+)"',sources))
directories=set(re.findall(r'new\s+G4UIdirectory\("(/[^"]+)"',sources))
# Standard Geant4 command prefixes are not project-owned and are excluded.
custom_prefixes=('/detector/','/readout/','/trigger/','/output/','/generator/','/CRY/')
errors=[]
for path in sorted((ROOT/'macros').glob('*.mac')):
    for lineno,line in enumerate(path.read_text().splitlines(),1):
        line=line.strip()
        if not line or line.startswith('#'): continue
        cmd=line.split()[0]
        if cmd.startswith(custom_prefixes) and cmd not in registered:
            errors.append(f'{path.name}:{lineno}: unregistered custom command {cmd}')
if errors:
    print('\n'.join(errors)); sys.exit(1)
print(f'Macro custom-command checks: OK ({len(registered)} registered commands; {len(directories)} UI directories)')
