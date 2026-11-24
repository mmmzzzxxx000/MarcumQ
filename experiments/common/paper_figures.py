# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Recreate all manuscript image files and numerical table exports."""
import json
import shutil
import subprocess
import sys
from marcum._paths import ROOT


def main():
    for task in ('plot-regions','plot-coefficients','plot-bound','plot-methods',
                 'plot-hybrid','plot-random-speedup','paper-audit'):
        subprocess.run([sys.executable,'-B',str(ROOT/'run.py'),task],cwd=str(ROOT),check=True)
    coverage=json.loads((ROOT/'docs/manuscript_coverage.json').read_text(encoding='utf-8'))
    folder=ROOT/'results/manuscript_figures'
    folder.mkdir(parents=True,exist_ok=True)
    for item in coverage['images']:
        source=ROOT/item['generated']
        assert source.is_file(),str(source)
        shutil.copy2(str(source),str(folder/item['manuscript_file']))
    print('All ten manuscript image files:',folder)
    print('Numerical tables:',ROOT/'results/tables')


if __name__=='__main__': main()
