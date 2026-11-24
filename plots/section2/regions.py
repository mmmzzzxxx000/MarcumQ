# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Algorithm 939 method regions; Python port of the original MATLAB figure.

Uses the original region-selection order, not the proposed hybrid strategy.
The large-order branch is displayed as AF, as in the original legend.
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
import numpy as np

METHODS = ('NI', 'SE', 'AF', 'RR', 'AF')
COLORS = ('#3b528b', '#22a884', '#95d840', '#fde725', '#95d840')
ORDERS = (1, 20)
EXTENT = 40.
RESOLUTION = 1600


def regions(a, b, nu):
    """Original 939 partition in the paper's (a,b) variables."""
    x, y = a*a/2., b*b/2.
    xi = a*b
    width = np.sqrt(4*x + 2*nu)
    transition = (y > x+nu-width) & (y < x+nu+width)
    return np.select([x < 30., (xi > 30.) & (nu*nu < 2*xi),
                      transition & (nu < 135.), transition & (nu >= 135.)],
                     [1, 2, 3, 4], default=0).astype(np.uint8)


def style():
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 10,
        'xtick.labelsize': 8, 'ytick.labelsize': 8, 'axes.linewidth': .65,
        'pdf.fonttype': 42, 'svg.fonttype': 'none',
    })


def draw(ax, nu, a, b):
    categories = regions(a[None, :], b[:, None], float(nu))
    cmap = ListedColormap(COLORS)
    ax.imshow(categories, origin='lower', extent=[0, EXTENT, 0, EXTENT],
              interpolation='nearest', cmap=cmap,
              norm=BoundaryNorm(np.arange(-.5, 5.5), cmap.N), aspect='equal')
    ax.axvline(np.sqrt(60.), color='#303030', linewidth=.8,
               linestyle=(0, (4, 3)))
    curve = np.linspace(0., np.sqrt(EXTENT*EXTENT-2*nu), 1000)
    ax.plot(curve, np.sqrt(curve*curve+2*nu), color='#303030',
            linewidth=.8, linestyle='-.')
    ax.set_xlim(0, EXTENT)
    ax.set_ylim(0, EXTENT)
    ax.set_xticks([0, 10, 20, 30, 40])
    ax.set_yticks([0, 10, 20, 30, 40])
    ax.set_xlabel(r'$a$', labelpad=3)
    ax.set_ylabel(r'$b$', labelpad=3)
    ax.tick_params(length=2.5, pad=3, width=.65)
    ids = sorted(np.unique(categories))
    legend = ax.legend(handles=[Patch(facecolor=COLORS[i], edgecolor='none',
                                    label=METHODS[i]) for i in ids],
                       loc='upper right', fontsize=8, frameon=True,
                       fancybox=False, framealpha=1., borderpad=.35,
                       labelspacing=.25, handlelength=1.25, handletextpad=.5)
    legend.get_frame().set_linewidth(.5)
    legend.get_frame().set_edgecolor('#777777')
    return {METHODS[i]: int(np.count_nonzero(categories == i)) for i in ids}


def main(output=None):
    if output is None:
        from marcum._paths import ROOT
        output = ROOT/'results/figures/section2'
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    style()
    grid = (np.arange(RESOLUTION)+.5)*EXTENT/RESOLUTION
    counts = {}
    for nu in ORDERS:
        width_inches, height_inches = 83/25.4, 80/25.4
        fig, ax = plt.subplots(figsize=(width_inches, height_inches))
        fig.subplots_adjust(left=.16, right=.97, bottom=.17, top=.97)
        counts[str(nu)] = draw(ax, nu, grid, grid)
        stem = output/('ab_nu%d' % nu)
        fig.savefig(str(stem)+'.pdf', dpi=600)
        fig.savefig(str(stem)+'.svg', dpi=600)
        fig.savefig(str(stem)+'.png', dpi=600)
        plt.close(fig)

    width_inches, height_inches = 183/25.4, 100/25.4
    fig, axes = plt.subplots(1, 2, figsize=(width_inches, height_inches))
    fig.subplots_adjust(left=.065, right=.985, bottom=.23, top=.97, wspace=.23)
    for i, (ax, nu) in enumerate(zip(axes, ORDERS)):
        draw(ax, nu, grid, grid)
        ax.text(.5, -.24, r'(%s) $\nu=%d$' % ('ab'[i], nu),
                ha='center', va='top', fontsize=10, transform=ax.transAxes)
    fig.savefig(str(output/'algorithm939_regions_preview.png'), dpi=300)
    plt.close(fig)

    # Independent translation of the original assignment-based MATLAB logic.
    aa, bb = np.meshgrid(np.linspace(.1, 40, 501), np.linspace(.1, 40, 501))
    for nu in (1., 10., 20., 135.):
        x, y, xi = aa**2/2, bb**2/2, aa*bb
        first = x < 30
        second = (xi > 30) & (nu**2 < 2*xi)
        transition = (y > x+nu-np.sqrt(4*x+2*nu)) & (y < x+nu+np.sqrt(4*x+2*nu))
        third, fourth = transition & (nu < 135), transition & (nu >= 135)
        original = np.zeros(aa.shape, dtype=np.uint8)
        original[first] = 1
        original[~first & second] = 2
        original[~first & ~second & third] = 3
        original[~first & ~second & ~third & fourth] = 4
        assert np.array_equal(original, regions(aa, bb, nu))
    assert counts['20'].get('RR', 0) > 0 and counts['1'].get('RR', 0) == 0
    assert not np.any(regions(aa, bb, 10.) == 3)
    (output/'regions.json').write_text(json.dumps(dict(
        orders=list(ORDERS), extent=[0, 40, 0, 40], resolution=RESOLUTION,
        counts=counts, algorithm='Original Algorithm 939 selection, no proposed hybrid guard',
        selection_validation_points=4*501*501, source_language='Python',
        guide_lines=['a=sqrt(60)', 'b=sqrt(a^2+2*nu)']), indent=2), encoding='utf-8')
    print(json.dumps(counts, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    main(args.output)
