# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Regenerate the 99 Section 6.2 references by independent Bessel integration."""
import argparse
from pathlib import Path
import mpmath as mp
from marcum._paths import ROOT
from common import write_csv


def reference(k, threshold, digits, method):
    with mp.workdps(digits):
        k, threshold = mp.mpf(k),mp.mpf(threshold)
        y = (1+k)*threshold
        return +mp.quad(lambda t: mp.exp(-k-t)*mp.besseli(0,2*mp.sqrt(k*t)),
                        [0,y/2,y],method=method)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--limit',type=int,default=99)
    args=parser.parse_args()
    if not 1 <= args.limit <= 99:
        parser.error('limit must be in [1,99]')
    rows=[]
    for k in list(range(5,36))+[29.5,30.5]:
        for threshold in (.1,.01,.001):
            with mp.workdps(110):
                high=reference(k,threshold,100,'tanh-sinh')
                low=reference(k,threshold,60,'gauss-legendre')
                difference=abs(high-low)/high
                assert difference < mp.mpf('1e-45')
                rows.append(dict(K=k,gamma_th=threshold,
                    reference=mp.nstr(high,100),relative_precision_check=mp.nstr(difference,20)))
            if len(rows)==args.limit:
                write_csv(ROOT/'results/section62/references_regenerated.csv',rows)
                print('Validated Bessel-integral references:',len(rows))
                return


if __name__=='__main__': main()
