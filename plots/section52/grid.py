# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Hybrid timing and method regions at orders 1 and 20 in all languages."""
from collections import defaultdict
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.patheffects as effects
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec
import numpy as np
from scipy.ndimage import distance_transform_edt, label as connected_regions


from marcum._paths import ROOT
OUT = ROOT / 'data/section52/grid'
FIGURES = ROOT / 'results/figures/section52/grid'
FIGURES.mkdir(parents=True, exist_ok=True)
ORDERS = (1, 20)
LANGUAGES = ('python', 'cpp', 'fortran')
STRATEGIES = ('algorithm_939', 'proposed')
LABELS = ('Algorithm 939', 'Proposed')
METHODS = ('se', 'af', 'ni', 'se_rr', 'rr')
METHOD_LABELS = ('SE', 'AF', 'NI', 'SE-RR', 'RR')


def boundaries(ids, edges):
    neighbors = defaultdict(set)

    def connect(p, q):
        neighbors[p].add(q)
        neighbors[q].add(p)

    for r, c in np.argwhere(ids[:, :-1] != ids[:, 1:]):
        connect((r,c+1), (r+1,c+1))
    for r, c in np.argwhere(ids[:-1] != ids[1:]):
        connect((r+1,c), (r+1,c+1))
    junctions = {point for point, others in neighbors.items() if len(others)!=2}
    paths = []
    for start in sorted(neighbors, key=lambda p: (p not in junctions, p)):
        while neighbors[start]:
            path = [start]
            while neighbors[path[-1]]:
                current = path[-1]
                following = min(neighbors[current])
                neighbors[current].remove(following)
                neighbors[following].remove(current)
                path.append(following)
                if following in junctions:
                    break
            paths.append([(edges[c],edges[r]) for r,c in path])
    return paths


def label_regions(ax, ids, grid):
    for method, label in enumerate(METHOD_LABELS):
        mask = ids==method
        if not mask.any():
            continue
        regions, _ = connected_regions(mask)
        sizes = np.bincount(regions.ravel())
        largest = 1+np.argmax(sizes[1:])
        for index in range(1,len(sizes)):
            if index!=largest and sizes[index]<45:
                continue
            interior = distance_transform_edt(np.pad(regions==index,1))[1:-1,1:-1]
            row, col = np.unravel_index(np.argmax(interior), ids.shape)
            ax.text(grid[col],grid[row],label, ha='center',va='center', color='#141414',
                    fontsize=7, fontweight='bold', clip_on=True,
                    path_effects=[effects.withStroke(linewidth=1.6,foreground='white')])


def plot_order(by_key, nu):
    grid = np.arange(1,41)/2
    edges = np.arange(.25,20.5,.5)
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
                         'font.size':8,'axes.labelsize':8,'xtick.labelsize':7,
                         'ytick.labelsize':7,'axes.linewidth':.60,
                         'pdf.fonttype':42,'svg.fonttype':'none'})
    width_mm, height_mm = 100, 116
    fig = plt.figure(figsize=(width_mm/25.4,height_mm/25.4))
    layout = GridSpec(3,3,width_ratios=[1,1,.105],wspace=.17,hspace=.24,
                      left=.16,right=.86,top=.93,bottom=.10)
    axes = np.array([[fig.add_subplot(layout[i,j]) for j in range(2)] for i in range(3)])
    cmap = plt.get_cmap('viridis_r')
    cmap.set_bad('#ffffff')
    cmap.set_over(cmap(1.))
    for i, (language, row_name, ceiling) in enumerate(zip(LANGUAGES, ('Python','C++','Fortran'), (100,2,2))):
        norm = Normalize(vmin=0,vmax=ceiling,clip=False)
        for k, strategy in enumerate(STRATEGIES):
            ax = axes[i,k]
            cells = [[by_key[language,nu,strategy,a,b] for a in grid] for b in grid]
            ids = np.asarray([[METHODS.index(p['method']) for p in line] for line in cells])
            values = np.asarray([[p['microseconds'] if p['setting'] is not None else np.nan
                                  for p in line] for line in cells])
            mesh = ax.imshow(np.ma.masked_invalid(values),origin='lower',interpolation='nearest',
                             extent=(.25,20.25,.25,20.25),cmap=cmap,norm=norm,aspect='equal')
            # Draw the diagonal exception as a single line over the measured cells.
            selected_af = (grid**2 > 30) & (nu**2 < 2*grid**2)
            selected_af &= (grid**2/2 >= 30) if strategy=='algorithm_939' else (grid > 7)
            diagonal = np.flatnonzero(selected_af)
            boundary_ids = ids.copy()
            assert np.all(ids[diagonal,diagonal]==METHODS.index('se'))
            boundary_ids[diagonal,diagonal] = METHODS.index('af')
            line = LineCollection(boundaries(boundary_ids,edges), colors='#222222', linewidths=.65,
                                  linestyles=(0,(3.6,2.1)), zorder=3)
            line.set_path_effects([effects.withStroke(linewidth=1.18,foreground='#ffffffcc')])
            ax.add_collection(line)
            if diagonal.size:
                extent = edges[[diagonal[0],diagonal[-1]+1]]
                ax.plot(extent,extent,color='#292929',linewidth=.58,linestyle=':',
                        path_effects=[effects.withStroke(linewidth=.98,foreground='#ffffffbd')],
                        zorder=3)
            label_regions(ax,ids,grid)
            if i==0:
                ax.set_title(LABELS[k],fontsize=8,pad=4)
            bad_r,bad_c = np.where(~np.isfinite(values))
            if bad_r.size:
                ax.scatter(grid[bad_c],grid[bad_r],marker='x',s=4,c='#242424',linewidths=.35,zorder=4)
            ax.set(xlim=(.25,20.25),ylim=(.25,20.25),xticks=[.5,10,20],yticks=[.5,10,20])
            ax.tick_params(length=2,pad=1.5,labelbottom=i==2,labelleft=k==0)
            if i==2:
                ax.set_xticklabels(['0.5','10','20'])
                ax.set_xlabel(r'$a$',labelpad=1)
            if k==0:
                ax.set_yticklabels(['0.5','10','20'])
                ax.set_ylabel(r'$b$',labelpad=1)
        cax = fig.add_subplot(layout[i,2])
        ticks = np.linspace(0,ceiling,5)
        bar = fig.colorbar(mesh,cax=cax,ticks=ticks,extend='max')
        bar.ax.set_yticklabels(['%g'%v for v in ticks])
        bar.ax.tick_params(labelsize=7,length=2,pad=2)
        bar.set_label(r'Time ($\mu$s)',fontsize=7,labelpad=3)
        pos = axes[i,0].get_position()
        fig.text(.024,(pos.y0+pos.y1)/2,row_name,va='center',ha='center',
                 rotation=90,rotation_mode='anchor',fontsize=8)
    fig.savefig(str(FIGURES/('hybrid_timing_nu%d.pdf'%nu)),dpi=600)
    fig.savefig(str(FIGURES/('hybrid_timing_nu%d.svg'%nu)),dpi=600)
    fig.savefig(str(FIGURES/('hybrid_timing_nu%d.png'%nu)),dpi=600)
    plt.close(fig)


def main():
    entries = json.loads((OUT/'timings.json').read_text(encoding='utf-8'))
    assert len(entries)==28800
    by_key = {(p['language'],p['nu'],p['strategy'],p['a'],p['b']):p for p in entries}
    assert len(by_key)==len(entries)
    for nu in ORDERS:
        plot_order(by_key, nu)


if __name__=='__main__':
    main()
