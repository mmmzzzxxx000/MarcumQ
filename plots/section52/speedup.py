# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Compare measured Fortran runtimes with a shared linear time axis."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

from marcum._paths import ROOT
HERE=ROOT/'results/figures/section52/random'
HERE.mkdir(parents=True,exist_ok=True)
data=json.loads((ROOT/'data/section52/random/summary.json').read_text(encoding='utf-8'))
overall=data['averages']['overall']
region=data['averages']['serr_region']
top=data['maxima']['serr_region']
groups=[
    ('Overall mean',overall['algorithm_939_mean_us'],overall['proposed_mean_us'],overall['speedup_ratio']),
    ('SE-RR region mean',region['algorithm_939_mean_us'],region['proposed_mean_us'],region['speedup_ratio']),
    ('Maximum observed speedup'+chr(10)+'in SE-RR region',top['algorithm_939_us'],top['proposed_us'],top['speedup']),
]
with (HERE/'speedup_figure_data.csv').open('w',newline='',encoding='ascii') as f:
    writer=csv.writer(f)
    writer.writerow(['comparison','algorithm_939_us','proposed_us','speedup'])
    writer.writerows(groups)
plt.rcParams.update({
    'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':9,'axes.labelsize':9,'xtick.labelsize':8,
    'pdf.fonttype':42,'svg.fonttype':'none',
    'axes.spines.top':False,'axes.spines.right':False,
})
fig,ax=plt.subplots(figsize=(7.15,3.95))
fig.subplots_adjust(left=.295,right=.86,bottom=.17,top=.87)
purple='#440154'
teal='#FDE725'
base_y=np.array([2.7,1.35,0.])
for y,(label,old,new,speedup) in zip(base_y,groups):
    ax.barh(y+.19,old,height=.30,color=purple,zorder=3)
    ax.barh(y-.19,new,height=.30,color=teal,zorder=3)
    ax.text(old+.025,y+.19,f'{old:.3f}',va='center',fontsize=8)
    ax.text(new+.025,y-.19,f'{new:.3f}',va='center',fontsize=8)
    ax.annotate('',xy=(new,y-.51),xytext=(old,y-.51),
                arrowprops={'arrowstyle':'->','color':'#3b3b3b','lw':1.1})
    ax.text((old+new)/2,y-.68,f'{speedup:.2f}'+r'$\times$',
            ha='center',va='center',fontweight='bold',fontsize=10,color=purple)
ax.set_yticks(base_y)
ax.set_yticklabels([row[0] for row in groups])
ax.tick_params(axis='y',length=0,pad=12,labelsize=9)
ax.set_xlim(0,1.7)
ax.set_xticks(np.arange(0,1.61,.4))
ax.set_ylim(-.93,3.15)
ax.set_xlabel('Computation time ('+r'$\mu$'+'s)',labelpad=8)
ax.xaxis.grid(True,color='#e6e6e6',linewidth=.6,zorder=0)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_color('#aaaaaa')
fig.legend(handles=[Patch(color=purple,label='Algorithm 939'),
                    Patch(color=teal,label='Proposed hybrid')],
           loc='upper center',bbox_to_anchor=(.59,.99),ncol=2,frameon=False,
           columnspacing=2,handlelength=1.5)
# Only the two mean rows share the five-batch statistic; the third is a single point.
fig.savefig(HERE/'fortran_random_speedup.png',dpi=600,facecolor='white')
fig.savefig(HERE/'fortran_random_speedup.pdf',facecolor='white')
fig.savefig(HERE/'fortran_random_speedup.svg',facecolor='white')
plt.close(fig)
