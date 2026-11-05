#!/usr/bin/env python3
"""Run the four bundled finite R6 checks using the Python standard library."""
import sys,subprocess
from pathlib import Path
if not __debug__:raise RuntimeError('Run without -O: assertions are verification conditions.')
HERE=Path(__file__).resolve().parent
for name in ('r6_main_checks.py','r6_additional_checks.py','algebra_checks.py','structural_checks.py'):
    print('\nRUN',name,flush=True)
    subprocess.run([sys.executable,str(HERE/name),'--check'],check=True)
print('ALL R6 FINITE CHECKS PASSED')
