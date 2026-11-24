# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Three bound maps with the discussed early branch switch at M=0.998."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.patheffects as effects
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import mpmath as mp
import numpy as np
from scipy.special import gammaincc, gammaln

from marcum._paths import ROOT
OUT = ROOT / 'results/figures/section4'
OUT.mkdir(parents=True, exist_ok=True)
ORDERS = (1, 10, 20)
A_MAX, B_MAX, STEP = 20., 8., .02
SWITCH, LIMIT = .998, 1000.


def log_m(a, b, nu):
    y = b*b/2
    assert (y > 0).all()
    return np.log(2.) - a*a/4 + y + nu*np.log(y) - gammaln(nu+1)


def upper(a, b, nu):
    lm = log_m(a, b, nu)
    m = np.exp(lm)
    first = lm < np.log(SWITCH)
    tail = np.broadcast_to(gammaincc(nu, b*b/2), m.shape)
    values = np.empty_like(m)
    values[first] = (1+m[first])/(-np.expm1(lm[first]))
    values[~first] = 2*(1+m[~first])/tail[~first]
    assert np.isfinite(values).all() and (values >= 1).all()
    return values, lm


def boundary(a, nu):
    # Each branch is increasing in b. The first branch stays below 999;
    # the accepted region therefore has one outer boundary after switching.
    roots = []
    for kind in ('switch', 'second'):
        low, high = np.zeros_like(a), np.full_like(a, B_MAX)
        for _ in range(58):
            mid = (low+high)/2
            lm = log_m(a, mid, nu)
            if kind == 'switch':
                smaller = lm < np.log(SWITCH)
            else:
                tail = gammaincc(nu, mid*mid/2)
                assert (tail > 0).all()
                log_second = np.log(2.) + np.logaddexp(0, lm) - np.log(tail)
                smaller = log_second < np.log(LIMIT)
            low = np.where(smaller, mid, low)
            high = np.where(smaller, high, mid)
        roots.append((low+high)/2)
    return np.maximum(*roots)


def main():
    a = (np.arange(int(A_MAX/STEP))+.5)*STEP
    b = (np.arange(int(B_MAX/STEP))+.5)*STEP
    aa, bb = np.meshgrid(a, b)
    data, checks, summaries = dict(a=a, b=b), [], []
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 10,
        'xtick.labelsize': 8, 'ytick.labelsize': 8, 'axes.linewidth': .65,
        'pdf.fonttype': 42, 'svg.fonttype': 'none',
    })
    width_mm, height_mm = 183., 65.
    fig, axes = plt.subplots(1, 3, figsize=(width_mm/25.4, height_mm/25.4),
                             sharex=True, sharey=True)
    fig.subplots_adjust(left=.065, right=.84, bottom=.16, top=.83,
                        wspace=.15)
    cmap = plt.get_cmap('viridis_r')
    cmap.set_over(cmap(1.))
    for ax, nu in zip(axes, ORDERS):
        values, lm = upper(aa, bb, nu)
        data['B_nu%d' % nu] = values
        accepted = values <= LIMIT
        assert (np.diff(accepted.astype(int), axis=0) <= 0).all()
        edge = boundary(a, nu)
        assert np.array_equal(accepted, b[:, None] <= edge[None, :])
        data['boundary_nu%d' % nu] = edge
        candidates = [(0, 0), (len(b)-1, 0), (0, len(a)-1),
                      (len(b)-1, len(a)-1)]
        for branch in (lm < np.log(SWITCH), lm >= np.log(SWITCH)):
            candidates.append(np.unravel_index(np.argmin(np.where(
                branch, np.abs(lm-np.log(SWITCH)), np.inf)), lm.shape))
        candidates.append(np.unravel_index(np.argmin(np.abs(values-LIMIT)), lm.shape))
        with mp.workdps(80):
            for row, col in candidates:
                av, bv = mp.mpf(float(a[col])), mp.mpf(float(b[row]))
                y = bv*bv/2
                m = 2*mp.exp(-av*av/4+y)*y**nu/mp.gamma(nu+1)
                exact = ((1+m)/(1-m) if m < mp.mpf('0.998') else
                         2*(1+m)*mp.gamma(nu)/mp.gammainc(nu, y, mp.inf))
                error = float(abs(mp.mpf(float(values[row, col]))/exact-1))
                assert error < 1e-9
                checks.append(dict(nu=nu, a=float(av), b=float(bv), relative_difference=error))
        mesh = ax.imshow(values, origin='lower', extent=(0, A_MAX, 0, B_MAX),
                         interpolation='nearest', aspect='auto', cmap=cmap,
                         norm=Normalize(1., 10000., clip=False))
        line_a = np.linspace(0, A_MAX, 2001)
        line_b = boundary(line_a, nu)
        visible = line_b < B_MAX-1e-10
        line, = ax.plot(line_a[visible], line_b[visible], color='#181818',
                        linewidth=.9, linestyle='--')
        line.set_path_effects([effects.withStroke(linewidth=1.8, foreground='white')])
        ax.set(xlim=(0, A_MAX), ylim=(0, B_MAX), xticks=[0, 10, 20],
               yticks=[0, 2, 4, 6, 8], xlabel=r'$a$')
        ax.set_title(r'$\nu=%d$' % nu, pad=7)
        ax.tick_params(length=2.5, pad=3)
        summaries.append(dict(nu=nu, grid_points=int(values.size),
                              at_most_1000=int(accepted.sum()),
                              above_10000=int((values > 10000).sum())))
    axes[0].set_ylabel(r'$b$', labelpad=3)
    cax = fig.add_axes([.87, .16, .017, .67])
    bar = fig.colorbar(mesh, cax=cax, extend='max', ticks=[1, 2000, 4000, 6000, 8000, 10000])
    bar.ax.set_yticklabels(['1', '2000', '4000', '6000', '8000', '10000'])
    bar.ax.tick_params(labelsize=8, length=2.5, pad=3)
    bar.ax.set_title(r'$B_{a,b,\nu}$', fontsize=10, pad=10)
    stem = OUT / 'B_bound_nu1_10_20'
    fig.savefig(str(stem)+'.png', dpi=600)
    fig.savefig(str(stem)+'.pdf', dpi=450)
    fig.savefig(str(stem)+'.svg', dpi=450)
    plt.close(fig)
    np.savez_compressed(str(OUT/'source_data.npz'), **data)
    metadata = dict(branch_switch=SWITCH, first_branch='(1+M)/(1-M), for M<0.998',
                    second_branch='2*(1+M)/gammaincc(nu,b^2/2), for M>=0.998',
                    domain=dict(a=[0, A_MAX], b=[0, B_MAX]), step=STEP,
                    shared_linear_color_scale=[1, 10000], region_boundary=LIMIT,
                    boundary_validation='All 1,200,000 grid classifications agree with the outer boundary.',
                    summaries=summaries, independent_80_digit_checks=checks,
                    manuscript_and_algorithm_files_modified=False)
    (OUT/'metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    print(json.dumps(summaries, indent=2))


if __name__ == '__main__':
    main()
