# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Deterministic validation plots; every panel is traceable to steps.csv.

No resampling, uncertainty intervals, interpolation, or data-based case selection.
Exact-zero errors have no logarithm and are recorded rather than floored.
"""
import csv
import json
from pathlib import Path
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp

HERE = Path(__file__).resolve().parent
from marcum._paths import ROOT
OUT = ROOT/'results/section3/recurrence'

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                     'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
                     'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': .7, 'legend.frameon': False,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none',
                     'mathtext.fontset': 'dejavusans', 'savefig.facecolor': 'white'})
BLUE = '#2166AC'
ORANGE = '#C76D26'
GRAY = '#72777E'


def load():
    with (OUT/'steps.csv').open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    data = {}
    for row in rows:
        row['n'] = int(row['n'])
        row['ratio'] = float(row['error_over_maximum_local'])
        row['util'] = float(row['error_over_H_bound'])
        row['H'] = float(row['H_n'])
        row['C'] = float(row['C_a_nu'])
        row['post'] = row['after_first_tiny'] == 'True'
        data.setdefault(row['case'], []).append(row)
    return rows, data


def save(fig, name):
    fig.savefig(OUT/(name+'.pdf'))
    fig.savefig(OUT/(name+'.svg'))
    fig.savefig(OUT/(name+'.png'), dpi=600)
    plt.close(fig)


def decorate(ax, letter):
    ax.set_xlabel('Recurrence step n')
    ax.grid(axis='y', color='#E8EBEE', linewidth=.5)
    ax.set_axisbelow(True)
    ax.text(-.16, 1.10, letter, transform=ax.transAxes, fontsize=11, fontweight='bold')


def main():
    rows, data = load()
    width_inches = 183/25.4
    fig, axes = plt.subplots(1, 2, figsize=(width_inches, 85/25.4))
    fig.subplots_adjust(left=.09, right=.98, bottom=.20, top=.83, wspace=.38)
    ax = axes[0]
    representative = data['grid_a1_v1']
    x = [r['n'] for r in representative]
    ax.plot(x, [r['ratio'] for r in representative], color=BLUE, lw=1.6,
            label='Actual error / maximum local error', zorder=4)
    ax.plot(x, [r['H'] for r in representative], color=ORANGE, lw=1.3, label=r'$H_n$')
    ax.axhline(representative[0]['C'], color=GRAY, ls='--', lw=1.2, label=r'$C_{a,\nu}$')
    tiny = next((r['n'] for r in representative if r['post']), None)
    if tiny:
        ax.axvspan(tiny, 256, color='#EEF0F2', zorder=0)
        ax.text(tiny+3, .10, 'Underflow range', color=GRAY, fontsize=8)
    ax.set_title(r'Representative: $a=1$, $\nu=1$', pad=12)
    ax.set_ylabel('Absolute-error amplification')
    ax.set_ylim(-.05, max(representative[0]['C'], max(r['H'] for r in representative))*1.15)
    ax.set_xlim(2, 256)
    ax.legend(loc='center right', bbox_to_anchor=(.99, .68), handlelength=2, borderaxespad=.3)
    decorate(ax, 'a')

    ax = axes[1]
    xs, maximum, median, post_max = [], [], [], []
    for n in range(2, 257):
        before = [float(r['error_over_H_bound']) for r in rows if r['n'] == n and not r['post']]
        after = [float(r['error_over_H_bound']) for r in rows if r['n'] == n and r['post']]
        xs.append(n)
        maximum.append(max(before) if before else np.nan)
        median.append(float(np.median(before)) if before else np.nan)
        post_max.append(max(after) if after else np.nan)
    ax.axhline(1, color=GRAY, ls=':', lw=1.2, label='Theoretical upper limit')
    ax.plot(xs, maximum, color=BLUE, lw=1.5, label='Batch maximum: before underflow')
    ax.plot(xs, median, color=ORANGE, lw=1.2, label='Batch median: before underflow')
    ax.plot(xs, post_max, color=GRAY, lw=1.1, ls='--', label='After underflow: maximum')
    ax.set_title('%d parameter pairs' % len(data), pad=12)
    ax.set_ylabel('Actual error / bound from '+r'$H_n$')
    ax.set_ylim(-.03, 1.1)
    ax.set_xlim(2, 256)
    ax.legend(loc='center right', bbox_to_anchor=(1, .6), fontsize=7.5, handlelength=1.7)
    decorate(ax, 'b')
    save(fig, 'coefficient_error_validation')

    selected = [('grid_a1_v1', 1., 1.), ('grid_a5_v0.25', 5., .25),
                ('grid_a20_v1', 20., 1.), ('grid_a100_v1', 100., 1.),
                ('grid_a5_v-0.5', 5., -.5), ('grid_a20_v-2.3', 20., -2.3)]
    fig, axes = plt.subplots(2, 3, figsize=(width_inches, 132/25.4))
    fig.subplots_adjust(left=.09, right=.98, bottom=.14, top=.87, hspace=.55, wspace=.4)
    zero_count = 0
    for ax, (label, a, v), letter in zip(axes.flat, selected, 'abcdef'):
        chosen = data[label]
        x, ratio_log, h_log = [], [], []
        for row in chosen:
            x.append(row['n'])
            error = mp.mpf(row['actual_absolute_error'])
            delta = mp.mpf(row['maximum_local_error'])
            if error > 0 and delta > 0:
                ratio_log.append(float(mp.log10(error)-mp.log10(delta)))
            else:
                ratio_log.append(np.nan)
                zero_count += 1
            h_log.append(math.log10(row['H']))
        ax.plot(x, ratio_log, color=BLUE, lw=1.2, label='Actual amplification')
        ax.plot(x, h_log, color=ORANGE, lw=1.2, label=r'$H_n$')
        ax.axhline(math.log10(chosen[0]['C']), color=GRAY, ls='--', lw=1, label=r'$C_{a,\nu}$')
        tiny = next((r['n'] for r in chosen if r['post']), None)
        if tiny:
            ax.axvspan(tiny, 256, color='#EEF0F2', zorder=0)
        ax.set_title(r'$a=%g$, $\nu=%g$' % (a, v), pad=9)
        ax.set_ylabel('log10 amplification')
        ax.set_xlim(2, 256)
        decorate(ax, letter)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=3, bbox_to_anchor=(.52, .995))
    save(fig, 'coefficient_error_examples')
    qa = dict(parameter_pairs=len(data), observations=len(rows),
              exact_zero_log_ordinates_omitted=zero_count,
              all_data_retained_in_csv=True, zero_floor_used=False,
              main_panel_a='Predetermined a=1, nu=1; linear scale retains zero errors',
              main_panel_b='All cases; pointwise max/median, not a confidence interval',
              shaded_ranges='At and after first subnormal or rounded-zero coefficient',
              display='Supplement uses log10 of decimal diagnostic strings; all tests were exact rational')
    (OUT/'figure_qa.json').write_text(json.dumps(qa, indent=2), encoding='utf-8')
    print(json.dumps(qa, indent=2))


if __name__ == '__main__':
    main()
