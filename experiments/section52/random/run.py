# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Fortran timing at matched accuracy with independent C++ reference values."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from collections import Counter
import ctypes as ct
import hashlib
import json
import os
from pathlib import Path
import time

os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
from scipy.special import gammaincc, gammaln

SOURCE = Path(__file__).resolve().parent
from marcum._paths import ROOT
HERE = ROOT / 'results/section52/random'
LIBRARIES = ROOT / 'build/hybrid'
SEED = 2026111702
CHUNK = 100_000
EPS = 1e-10
REPEATS = 5
METHODS = ('SE','AF','NI','SE-RR','RR','AF-large-nu')
LOW = [0.01,0.001,1.]
HIGH = [100.,100.,1000.]
DBL = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
INT = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS')
ref = kernel = None


def load():
    global ref, kernel
    ref = ct.CDLL(str(LIBRARIES/'reference.dll'))
    ref.reference_batch.argtypes = [ct.c_int,DBL,ct.c_int,DBL]
    ref.reference_batch.restype = None
    kernel = ct.CDLL(str(LIBRARIES/'measure.dll'))
    kernel.calibrate.argtypes = [ct.c_int]+[DBL]*6+[INT]*4
    kernel.calibrate.restype = None
    kernel.batch_seconds.argtypes = [ct.c_int]*2+[DBL]*3+[INT]*2
    kernel.batch_seconds.restype = ct.c_double
    kernel.screen.argtypes = [ct.c_int]*2+[DBL]*3+[INT]*4+[DBL]
    kernel.screen.restype = None
    kernel.rebuilt_time.argtypes = [ct.c_int]+[ct.c_double]*3+[ct.c_int64]*2
    kernel.rebuilt_time.restype = ct.c_double


def vector(x):
    return np.ascontiguousarray(x)


def routes(a,b,nu):
    x,y = a*a/2,b*b/2
    xi=2*np.sqrt(x*y)
    w=np.sqrt(4*x+2*nu)
    old=np.select([x<30,(xi>30)&(nu*nu<2*xi),(y>=x+nu-w)&(y<=x+nu+w)],
                  [0,1,np.where(nu<135,4,5)],default=2).astype(np.int32)
    lm=np.log(2.)-x/2+y+nu*np.log(y)-gammaln(nu+1)
    region=(b<=7)&(lm<np.log(.998))
    mask=(b<=7)&~region
    region[mask] = 2*(1+np.exp(lm[mask])) <= 1000*gammaincc(nu[mask],y[mask])
    new=np.select([region,(b<=7)|(a<=7),(a*b>30)&(nu*nu<2*a*b),(y>x+nu-w)&(y<x+nu+w)],
                  [3,0,1,np.where(nu<135,4,5)],default=2).astype(np.int32)
    old[(a==b)&(old==1)]=0
    new[(a==b)&(new==1)]=0
    return old,new


def validate(job):
    index,points,location,ref_threads = job
    path=location/f'chunk_{index:03d}.npz'
    if path.exists():
        with np.load(path) as data:
            return index,len(data['points']),int((data['settings']>=0).all(axis=1).sum()),'saved'
    start=time.perf_counter()
    n=len(points)
    reference=np.empty((n,5),dtype=np.float64)
    ref.reference_batch(n,points,ref_threads,reference)
    a,b,nu=map(vector,points.T)
    old,new=routes(a,b,nu)
    reliable=reference[:,4]>0
    below=reference[:,2]<np.log(1e-290)
    eligible=reliable&~below
    hi,lo=map(vector,reference[eligible,:2].T)
    bound=vector((hi+lo)*reference[eligible,3])
    old_n=np.empty(eligible.sum(),np.int32)
    new_n=np.empty_like(old_n)
    if len(old_n):
        kernel.calibrate(len(old_n),a[eligible],b[eligible],nu[eligible],
                         hi,lo,bound,old[eligible],new[eligible],old_n,new_n)
    settings=np.full((n,2),-1,dtype=np.int32)
    settings[eligible,0]=old_n
    settings[eligible,1]=new_n
    temp=path.with_suffix('.tmp.npz')
    np.savez(temp,points=points,reference=reference,routes=np.column_stack((old,new)),
             settings=settings)
    temp.replace(path)
    return index,n,int((settings>=0).all(axis=1).sum()),round(time.perf_counter()-start,2)


def write_json(path,data):
    temp=path.with_suffix('.tmp')
    temp.write_text(json.dumps(data,indent=2),encoding='utf-8')
    temp.replace(path)


def prepare(folder,count,workers):
    location=folder/'chunks'
    location.mkdir(exist_ok=True)
    start=time.perf_counter()
    rng=np.random.default_rng(SEED)
    jobs=((offset//CHUNK,
           np.exp(rng.uniform(np.log(LOW),np.log(HIGH),(min(CHUNK,count-offset),3))),
           location,max(1,(os.cpu_count() or 16)//workers))
          for offset in range(0,count,CHUNK))
    completed=paired=0
    with ProcessPoolExecutor(max_workers=workers,initializer=load) as pool:
        for index,n,passed,elapsed in pool.map(validate,jobs,chunksize=1):
            completed+=n
            paired+=passed
            if index<8 or (index+1)%5==0 or completed==count:
                print(f'Calibrated {completed:,}/{count:,}; paired {paired:,}; last chunk {elapsed}s; '
                      f'wall {time.perf_counter()-start:.1f}s',flush=True)
    return time.perf_counter()-start


def inspect(folder, count):
    totals=Counter()
    routing=np.zeros((len(METHODS),len(METHODS)),np.int64)
    for index in range((count+CHUNK-1)//CHUNK):
        with np.load(folder/'chunks'/f'chunk_{index:03d}.npz') as data:
            n=len(data['points'])
            old,new=map(vector,data['routes'].T)
            on,nn=map(vector,data['settings'].T)
            references=data['reference']
        good_ref=references[:,4]>0
        below=references[:,2]<np.log(1e-290)
        eligible=good_ref&~below
        valid=(on>=0)&(nn>=0)
        region=new==3
        totals.update(drawn=n,paired=int(valid.sum()),region_drawn=int(region.sum()),
                      region_paired=int((region&valid).sum()),
                      reference_unresolved=int((~good_ref).sum()),
                      below_939_cutoff=int((good_ref&below).sum()),
                      above_floor_939_not_qualified=int((eligible&(on<0)).sum()),
                      above_floor_proposed_not_qualified=int((eligible&(nn<0)).sum()),
                      above_floor_not_paired=int((eligible&~valid).sum()),
                      changed_paired=int((valid&(old!=new)).sum()))
        np.add.at(routing,(old[valid],new[valid]),1)
    assert totals['drawn']==count
    assert totals['paired']+totals['reference_unresolved']+totals['below_939_cutoff']+totals['above_floor_not_paired']==count
    return dict(totals),routing.tolist()


def workloads(folder,index):
    with np.load(folder/'chunks'/f'chunk_{index:03d}.npz') as source:
        points,routes,settings=[source[k] for k in ('points','routes','settings')]
    valid=(settings>=0).all(axis=1)
    return points,routes,settings,valid,routes[:,1]==3


def average_times(folder,count,counts):
    existing=folder/'averages.json'
    if existing.exists():
        return json.loads(existing.read_text(encoding='utf-8'))
    samples=np.zeros((REPEATS,2,2))
    for rep in range(REPEATS):
        start=time.perf_counter()
        for index in range((count+CHUNK-1)//CHUNK):
            points,routes,settings,valid,region=workloads(folder,index)
            side_order=(0,1) if (rep+index)%2==0 else (1,0)
            for group,mask in enumerate((valid,valid&region)):
                xyz=[vector(v) for v in points[mask].T]
                for side in side_order:
                    samples[rep,group,side]+=kernel.batch_seconds(len(xyz[0]),2,*xyz,
                                                               routes[mask,side],settings[mask,side])
        print(f'Average timing pass {rep+1}/{REPEATS}: {time.perf_counter()-start:.1f}s',flush=True)
    result={}
    for group,(key,n) in enumerate((('overall',counts['paired']),('serr_region',counts['region_paired']))):
        old,new=np.median(samples[:,group,:],axis=0)*1e6/n
        ratio=old/new
        result[key]=dict(points=n,algorithm_939_mean_us=old,proposed_mean_us=new,
                         speedup_ratio=ratio,speedup_percent=(ratio-1)*100,
                         totals_seconds=samples[:,group,:].tolist())
    write_json(existing,result)
    return result


def screen_chunk(folder,index):
    path=folder/'screen'/f'candidates_{index:03d}.json'
    if path.exists():
        return json.loads(path.read_text(encoding='utf-8'))
    points,routes,settings,valid,region=workloads(folder,index)
    selected=valid&(routes[:,0]!=routes[:,1])
    xyz=list(map(vector,points[selected].T))
    knobs=list(map(vector,routes[selected].T))+list(map(vector,settings[selected].T))
    ratios=np.empty(len(xyz[0]))
    if len(ratios):
        kernel.screen(len(ratios),8,*xyz,*knobs,ratios)
    by_ratio=np.argsort(ratios)
    indices=set(by_ratio[-16:].tolist())
    region_filtered=region[selected]
    indices.update(by_ratio[region_filtered[by_ratio]][-16:].tolist())
    rows=[dict(screen=float(ratios[i]),point=[float(v[i]) for v in xyz],
               old=int(knobs[0][i]),new=int(knobs[1][i]),
               old_setting=int(knobs[2][i]),new_setting=int(knobs[3][i]),
               region=bool(region_filtered[i])) for i in indices]
    write_json(path,rows)
    return rows


def maximum_times(folder,count):
    start=time.perf_counter()
    (folder/'screen').mkdir(exist_ok=True)
    candidates=[]
    for index in range((count+CHUNK-1)//CHUNK):
        candidates.extend(screen_chunk(folder,index))
        if (index+1)%5==0:
            print(f'Maximum screening {(index+1)*CHUNK:,}/{count:,}; {time.perf_counter()-start:.1f}s',flush=True)
    ranked=sorted(candidates,key=lambda row:row['screen'],reverse=True)
    finalists=ranked[:128]+[row for row in ranked if row['region']][:128]
    unique={}
    for row in finalists:
        unique[tuple(row['point'])]=row
    candidates=list(unique.values())
    print(f'Retesting {len(candidates)} top candidates after screening every changed point.',flush=True)
    reevaluated=[]
    for row in candidates:
        samples=np.empty((9,2))
        for rep in range(len(samples)):
            for side in ((0,1) if rep%2==0 else (1,0)):
                name=('old','new')[side]
                samples[rep,side]=1e6*kernel.rebuilt_time(row[name],*row['point'],row[name+'_setting'],512)
        old,new=np.median(samples,axis=0)
        reevaluated.append(dict(row,algorithm_939_method=METHODS[row['old']],
                                proposed_method=METHODS[row['new']],algorithm_939_us=old,
                                proposed_us=new,speedup=old/new,samples_us=samples.tolist()))
    reevaluated.sort(key=lambda row:row['speedup'],reverse=True)
    write_json(folder/'remeasured_candidates.json',reevaluated)
    return dict(overall=reevaluated[0],serr_region=next(r for r in reevaluated if r['region']),
                remeasured=len(reevaluated),wall_seconds=time.perf_counter()-start)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--count',type=int,default=10_000_000)
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--build-only',action='store_true')
    args=parser.parse_args()
    if args.count < 1 or args.workers < 1:
        parser.error('count and workers must be positive')
    from build_hybrid import build
    build()
    if args.build_only:
        return
    folder=HERE/f'sample_{args.count}'
    folder.mkdir(parents=True,exist_ok=True)
    metadata=dict(seed=SEED,count=args.count,low=LOW,high=HIGH,distribution='independent log-uniform',
                  reference='C++ x87 extended-precision defining Bessel integral; probability bound for Q near one',
                  comparison='Fortran evaluation at the first tested setting meeting 1e-10 relative error; routing excluded',
                  target_relative_error=EPS,methods=METHODS,chunk_size=CHUNK,
                  binaries={name:hashlib.sha256((LIBRARIES/name).read_bytes()).hexdigest()
                            for name in ('reference.dll','measure.dll')})
    meta_path=folder/'protocol.json'
    if meta_path.exists():
        assert json.loads(meta_path.read_text(encoding='utf-8'))==json.loads(json.dumps(metadata))
    else:
        write_json(meta_path,metadata)
    validation_seconds=prepare(folder,args.count,args.workers)
    counts,routing=inspect(folder,args.count)
    write_json(folder/'qualification.json',dict(counts=counts,routing_counts=routing,
                                               validation_seconds=validation_seconds))
    print(f'Qualification: {counts}',flush=True)
    load()
    averages=average_times(folder,args.count,counts)
    print(f'Average speedups: {[(key, row["speedup_ratio"]) for key,row in averages.items()]}',flush=True)
    maxima=maximum_times(folder,args.count)
    result=dict(metadata=metadata,counts=counts,routing_counts=routing,averages=averages,maxima=maxima,
                validation_seconds=validation_seconds)
    write_json(folder/'summary.json',result)
    print(f'Maximum confirmed on retest, SE-RR region: {maxima["serr_region"]["speedup"]:.6g}x',flush=True)


if __name__=='__main__':
    main()
