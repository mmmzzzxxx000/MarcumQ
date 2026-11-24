# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Audit the recorded data and export the numerical tables in the manuscript."""
import csv
import json
import statistics
from pathlib import Path
from marcum._paths import ROOT
from common import write_csv, write_json


def read_csv(path):
    with path.open(newline='',encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle))


def main():
    output=ROOT/'results/tables'
    output.mkdir(parents=True,exist_ok=True)
    checks=[]
    methods={}
    with (ROOT/'data/section51/timings.jsonl').open(encoding='utf-8') as handle:
        for line in handle:
            row=json.loads(line)
            if row['setting'] is not None:
                key=(row['nu'],row['a'],row['b'])
                methods.setdefault(key,{})[(row['language'],row['method'])]=row['microseconds']
    table=[]
    for nu,count in [(1,658),(10,662),(20,667)]:
        common=[rows for (v,a,b),rows in methods.items() if v==nu and len(rows)==12]
        assert len(common)==count
        for language in ('python','cpp','fortran'):
            row=dict(nu=nu,language=language,points=count)
            for method in ('se','af','ni','se_rr'):
                row[method+'_us']=statistics.mean(r[(language,method)] for r in common)
            for method in ('se','af','ni'):
                row['speedup_vs_'+method]=row[method+'_us']/row['se_rr_us']
            table.append(row)
    write_csv(output/'individual_methods.csv',table)
    checks.append('Section 5.1: three languages, four methods; common sets contain 658, 662, and 667 points.')
    hybrid=json.loads((ROOT/'data/section52/grid/timings.json').read_text(encoding='utf-8'))
    assert len(hybrid)==28800 and all(r['setting'] is not None and r['relative_error']<=1e-10 for r in hybrid)
    selected=[r for r in hybrid if r['language']=='fortran' and r['nu']==1]
    region={(r['a'],r['b']) for r in selected if r['strategy']=='proposed' and r['method']=='se_rr'}
    grid_means={}
    for strategy,whole,part in [('algorithm_939',.608,.737),('proposed',.458,.155)]:
        rows=[r for r in selected if r['strategy']==strategy]
        grid_means[strategy]=dict(overall=statistics.mean(r['microseconds'] for r in rows),
            region=statistics.mean(r['microseconds'] for r in rows if (r['a'],r['b']) in region))
        assert round(grid_means[strategy]['overall'],3)==whole and round(grid_means[strategy]['region'],3)==part
    write_json(output/'hybrid_grid.json',grid_means)
    checks.append('Section 5.2 grid: all hybrid cells accurate; Fortran nu=1 regional and overall means match the text.')
    cdf=json.loads((ROOT/'data/section61/summary.json').read_text(encoding='utf-8'))['pairs']
    assert len(cdf)==12
    assert all(r['points']==1000 and r['interval']=='0_to_5' for r in cdf)
    assert all(r['gammachi_max_error']<=1e-10 and r['hybrid_max_error']<=1e-10 for r in cdf)
    assert round(statistics.mean(r['speedup'] for r in cdf if r['mode']=='scalar'),2)==1.41
    assert round(statistics.mean(r['speedup'] for r in cdf if r['mode']=='batch'),2)==1.82
    write_csv(output/'cdf.csv',cdf)
    checks.append('Section 6.1: six configurations, 1000 thresholds, both modes at relative error 1e-10; mean speedups 1.41 and 1.82.')

    reference=read_csv(ROOT/'data/section62/references.csv')
    timings=read_csv(ROOT/'data/section62/timings.csv')
    expected={(float(k),g) for k in list(range(5,36))+[29.5,30.5] for g in (.1,.01,.001)}
    keys=lambda rows:{(float(r['K']),float(r['gamma_th'])) for r in rows}
    assert len(reference)==len(timings)==len(expected)==99
    assert keys(reference)==keys(timings)==expected
    assert all(float(r['algorithm939_rel_error'])<=1e-10 and float(r['se_rr_rel_error'])<=1e-10 for r in timings)
    old=statistics.mean(float(r['algorithm939_us']) for r in timings)
    new=statistics.mean(float(r['se_rr_us']) for r in timings)
    assert round(old,3)==1.649 and round(new,3)==.174 and round(old/new,2)==9.48
    write_csv(output/'outage.csv',timings)
    write_json(output/'outage_summary.json',dict(points=99,algorithm939_mean_us=old,hybrid_mean_us=new,speedup=old/new))
    checks.append('Section 6.2: all 99 points, target accuracy met, means 1.649 and 0.174 us, speedup 9.48.')

    random=json.loads((ROOT/'data/section53/summary.json').read_text(encoding='utf-8'))['totals']
    assert random['drawn']==100000000 and random['admitted']==94671159
    assert random['passed']==random['admitted'] and random['failed']==random['invalid']==0
    checks.append('Section 5.3: 100 million draws; 94,671,159 admitted and passed.')

    timing=json.loads((ROOT/'data/section52/random/summary.json').read_text(encoding='utf-8'))
    assert timing['metadata']['low']==[.01,.001,1.] and timing['metadata']['high']==[100.,100.,1000.]
    assert timing['metadata']['count']==10000000
    assert timing['counts']['paired']==9399065 and timing['counts']['below_939_cutoff']==600935
    assert timing['counts']['region_paired']==timing['counts']['region_drawn']==7402955
    assert round(timing['averages']['overall']['speedup_percent'],1)==82.2
    assert round(timing['averages']['serr_region']['speedup_percent'],1)==156.6
    assert round(timing['maxima']['serr_region']['speedup'],2)==5.49
    write_json(output/'random_timing.json',timing)
    checks.append('Section 5.2: correct full-domain 10-million sample; 82.2%, 156.6%, maximum 5.49x.')

    cancellation=json.loads((ROOT/'data/section4/cancellation.json').read_text(encoding='utf-8'))
    assert cancellation['fixed_N']==128 and len(cancellation['cases'])==3
    write_csv(output/'cancellation.csv',cancellation['cases'])
    checks.append('Section 4: three finite-sum cancellation examples at N=128.')
    write_json(output/'manuscript_audit.json',dict(passed=True,checks=checks,
        timings='Recorded paper measurements are preserved. New runs are hardware dependent.',
        scope='Manuscript figures and experiments only; response-only component profiling is excluded.'))
    for item in checks: print('PASS:',item)


if __name__=='__main__': main()
