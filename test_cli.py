#!/usr/bin/env python3
"""Regression checks for the exact, low-overhead automatic display cutoff."""
import subprocess
from pathlib import Path

EXE = str(Path(__file__).resolve().parent / 'conway_fast')

def run(*args):
    proc = subprocess.run([EXE, *map(str, args)], text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    return proc.stdout

def has_solution(text):
    return ' X=' in text and ' Y=' in text

assert has_solution(run('61'))
assert not has_solution(run('61', '--no-print'))
assert has_solution(run('61', '--print'))
assert has_solution(run('61', '--print-limit', '10')) # 1766319049 = 10 digits
assert not has_solution(run('61', '--print-limit', '9'))
assert has_solution(run('61', '--print-limit', '9', '--print'))
assert not has_solution(run('2', '--print-limit', '0'))
assert has_solution(run('2', '--print-limit', '1'))
assert not has_solution(run('100460221'))
assert ' square_factor=' not in run('61')
for bad in ('-1', 'x', '12junk'):
    r = subprocess.run([EXE, '61', '--print-limit', bad],
                        stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    assert r.returncode != 0, bad
print('All automatic / forced / suppressed output tests passed.')
# The general-purpose interface must reject obsolete special-purpose options.
rejected = subprocess.run([EXE, '61', '--cattle'], text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
assert rejected.returncode != 0
assert 'unknown argument: --cattle' in rejected.stderr

# Compact output is an exact unit-plus-exponent representation retained
# internally; the large final coordinates are not materialized.
cattle = '410286423278424'
compact = run(cattle, '--no-print')
assert ' representation=compact' in compact
assert ' exponent=2329' in compact
assert ' square_factor=9314' in compact
assert ' X_bits=' not in compact and ' Y_bits=' not in compact
assert ' representation=compact' in run(cattle)  # auto mode skips huge expansion
expanded = run(cattle, '--no-print', '--expand')
assert ' representation=expanded' in expanded
assert ' X_bits=343065 Y_bits=343040' in expanded
assert ' representation=expanded' in run(cattle, '--no-print', '--verify')
# Small solutions continue to print automatically, including lifted inputs.
assert ' X=17 Y=4' in run('18')
assert ' representation=compact' in run('18','--no-print')
assert ' X=99 Y=14' in run('50')
print('Compact/expanded execution tests passed.')
