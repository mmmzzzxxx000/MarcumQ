# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Reproduce the original coefficient-decay figure using Python.

The parameters use mu (not nu). Values are those of the corollary with p=2.
The visual specification preserves the original MATLAB figure.
"""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, LogFormatterMathtext, NullFormatter
import numpy as np

CASES = ((2., 1.), (2., 5.), (10., 1.), (10., 5.), (20., 1.))
COLORS = ((0., .447, .741), (.85, .325, .098), (.929, .694, .125),
          (.494, .184, .556), (.466, .674, .188))
N_TERMS = 40


def coefficients(a, mu):
    """Original recurrence, including the original arithmetic grouping."""
    u = [0.] * N_TERMS
    u[0] = 1.
    u[1] = (a*a-2*(mu+1))/(4*(mu+1))
    for n in range(2, N_TERMS):
        denominator = 4*n*(mu+n)
        first = (a*a-2*(mu+2*n-1))/denominator
        second = -1./denominator
        u[n] = first*u[n-1]+second*u[n-2]
    return np.array(u)


def main(output=None, reference=None):
    if output is None:
        from marcum._paths import ROOT
        output = ROOT/'results/figures/section3'
        reference = ROOT/'data/section3/coefficient_decay.csv'
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        'font.family': 'serif', 'font.serif': ['Times New Roman', 'DejaVu Serif'],
        'mathtext.fontset': 'cm', 'font.size': 14.,
        'axes.linewidth': .5, 'axes.labelsize': 15.4,
        'xtick.labelsize': 14., 'ytick.labelsize': 14.,
        'pdf.fonttype': 42, 'svg.fonttype': 'none',
    })
    # The original EPS bounding box is 675 x 450 points.
    width_inches, height_inches = 675/72., 450/72.
    fig = plt.figure(figsize=(width_inches, height_inches))
    ax = fig.add_axes([.13, .11, .775, .815])
    ax.set_yscale('log')
    ax.set_xlim(0., 40.)
    ax.set_ylim(5.306575654968534e-62, 1e20)
    ax.set_xticks(np.arange(0, 41, 5))
    ax.yaxis.set_major_locator(FixedLocator(10.**np.arange(-60, 21, 20)))
    ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
    ax.yaxis.set_minor_locator(FixedLocator(10.**np.array([k for k in range(-60, 21, 5) if k%20])))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(which='major', direction='in', top=True, right=True, length=4, width=.5)
    ax.tick_params(which='minor', direction='in', top=True, right=True, length=2.5, width=.5)
    ax.set_axisbelow(True)
    ax.grid(True, which='major', color='#262626', linestyle='-', linewidth=.5, alpha=.15)
    ax.grid(True, which='minor', axis='y', color='#262626', linestyle=':', linewidth=.5, alpha=.15)
    rows=[]
    for (a, mu), color in zip(CASES, COLORS):
        u = coefficients(a, mu)
        magnitude = np.abs(u)
        # Zero is retained in the data but has no position on a logarithmic axis.
        shown = np.ma.masked_less_equal(magnitude, 0.)
        ax.plot(np.arange(N_TERMS), shown, color=color, marker='o', markersize=4.5,
                markeredgewidth=0., linewidth=2., clip_on=False,
                label=r'$a=%.1f,\,\mu=%.1f$' % (a, mu))
        rows += [dict(a=a, mu=mu, n=n, u=float(u[n]), abs_u=float(magnitude[n])) for n in range(N_TERMS)]
    ax.set_xlabel(r'Index $n$', labelpad=4)
    ax.set_ylabel(r'$|u_n|$ (Log Scale)', labelpad=2)
    legend = ax.legend(loc='lower left', bbox_to_anchor=(.15515257, .1475),
                       bbox_transform=fig.transFigure, borderaxespad=0.,
                       frameon=True, fancybox=False, framealpha=1., fontsize=14.,
                       borderpad=.2, labelspacing=.15, handlelength=2., handletextpad=.4)
    legend.get_frame().set_linewidth(.5)
    legend.get_frame().set_edgecolor('#666666')
    fig.savefig(str(output/'u_n.pdf'))
    fig.savefig(str(output/'u_n.svg'))
    fig.savefig(str(output/'u_n.png'), dpi=300)
    fig.savefig(str(output/'u_n_preview.png'), dpi=144)
    plt.close(fig)
    with (output/'coefficients.csv').open('w', encoding='utf-8', newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=['a','mu','n','u','abs_u'])
        writer.writeheader();writer.writerows(rows)
    checks=dict(curves=5, points_per_curve=N_TERMS, zero_values=sum(r['abs_u']==0 for r in rows),
                zero_display='Masked only on the log axis; unchanged values retained in CSV.',
                axis_order_parameter='mu', figure_size_points=[675,450])
    if reference is not None:
        with Path(reference).open(encoding='utf-8',newline='') as handle: expected=list(csv.DictReader(handle))
        assert len(expected)==len(rows)
        for actual,old in zip(rows,expected):
            assert (actual['a'],actual['mu'],actual['n'])==(float(old['a']),float(old['mu']),int(old['n']))
            assert np.isclose(actual['abs_u'],float(old['abs_u']),rtol=2e-12,atol=0.)
        checks['all_200_points_match_original_figure']=True
    (output/'validation.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
    print(json.dumps(checks,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--reference',type=Path)
    args=parser.parse_args()
    main(args.output,args.reference)
