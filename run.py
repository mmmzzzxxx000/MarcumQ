# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Source-checkout entry point. No experiments run on import."""
import argparse, os, runpy, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.dont_write_bytecode=True
sys.path[:0]=[str(ROOT/'src/python'),str(ROOT/'experiments/common'),
              str(ROOT/'plots/section4'),str(ROOT)]
COMMANDS={
 'section51':'experiments/section51/run.py',
 'section52':'experiments/section52/grid/run.py',
 'section53':'experiments/section53/run.py',
 'section61':'experiments/section61/run.py',
 'section62':'experiments/section62/run.py',
 'cdf':'experiments/section61/run.py',
 'outage':'experiments/section62/run.py',
 'outage-reference':'experiments/section62/reference.py',
 'recurrence':'experiments/section3/recurrence.py',
 'cancellation':'experiments/section4/cancellation.py',
 'random':'experiments/section53/run.py',
 'random-timing':'experiments/section52/random/run.py',
 'plot-random-speedup':'plots/section52/speedup.py',
 'paper-audit':'tests/manuscript_results.py',
 'paper-figures':'experiments/common/paper_figures.py',
 'plot-methods':'plots/section51/run.py',
 'plot-hybrid':'plots/section52/grid.py',
 'plot-bound':'plots/section4/bound.py',
 'plot-regions':'plots/section2/regions.py',
 'plot-coefficients':'plots/section3/coefficient_decay.py',
 'test':'tests/acceptance.py',
 'test-hybrid':'tests/hybrid_interfaces.py',
 'build':'experiments/common/build_all.py',
}
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('task',choices=COMMANDS)
    args,rest=parser.parse_known_args()
    path=ROOT/COMMANDS[args.task]
    sys.path.insert(0,str(path.parent))
    sys.argv=[str(path)]+rest
    runpy.run_path(str(path),run_name='__main__')
