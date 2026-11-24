# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Check relocated sources, native ports, reference accuracy and current routes."""
import ast, importlib.util, json, math
from pathlib import Path
import mpmath as mp
import numpy as np
from marcum import kernel, se_rr
from marcum._paths import ROOT
from marcum.application1 import Runtime
from source939 import verify_sources
from bound import upper

def module(name,path):
    import sys
    sys.path.insert(0,str(path.parent))
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj

def main():
    count=0
    verify_sources()
    for path in ROOT.rglob('*.py'):
        if 'build' not in path.parts: ast.parse(path.read_text(encoding='utf-8-sig'))
    cases=[(1.,1.,2.),(2.,2.,1.),(8.,5.,1.),(10.,7.,10.),(2.,1.,100.),(2.,1.,101.),(5.,3.,1000.)]
    for a,b,nu in cases:
        for n in (8,20,40):
            values=[kernel('se_rr',a,b,nu,l)[0](n) for l in ('python','cpp','fortran')]
            assert all(math.isfinite(v) for v in values),(a,b,nu,n,values)
            assert max(values)-min(values)<2e-12,(a,b,nu,n,values)
            count+=1
    # Held-out high-precision integral values from the same validated reference columns.
    for nu in (1,10,20):
        for a,b in ((1.,1.),(8.,2.),(15.,3.)):
            path=ROOT/'data/section51/references'/('nu%d'%nu)/('a_%03d.json'%(round(a*4)-1))
            ref=json.loads(path.read_text(encoding='utf-8'))
            with mp.workdps(90):
                truth=mp.mpf(ref['q'][ref['b'].index(b)])
                for language in ('python','cpp','fortran'):
                    evaluate,_=kernel('se_rr',a,b,float(nu),language)
                    assert any(abs(mp.mpf(evaluate(n))-truth)<=mp.mpf('1e-10')*abs(truth) for n in range(1,129)),(a,b,nu,language)
                    count+=1
    for language in ('python','cpp','fortran'):
        b=np.linspace(0,3,31);r=Runtime(language)
        scalar,_=r.prepare('se_rr','scalar',2.,b,1.,np.full(len(b),32))
        batch,_=r.prepare('se_rr','batch',2.,b,1.,np.full(len(b),32))
        assert np.max(np.abs(scalar()-batch()))<2e-14
        count+=1
    a,b=np.meshgrid(np.arange(1,41)/2,np.arange(1,41)/2)
    hybrid=module('hybrid_driver',ROOT/'experiments/section52/grid/run.py')
    for nu in (1,10,20):
        old,new=hybrid.plans(None,a,b,nu)
        bounds,_=upper(a,b,nu)
        assert np.array_equal(new==3,(b<=7)&(bounds<=1000))
        assert not np.any((new==1)&(a==b))
        assert not np.any((old==1)&(a==b))
        count+=1
    from hybrid_runtime import load
    native,_=load()
    cpp,_=hybrid.build_cpp_rr()
    for aa,bb,vv in ((8.,10.,20.),(10.,12.,20.)):
        for setting in (5,7):
            values=[hybrid.python_rr_kernel(aa,bb,vv)[0](setting),
                    cpp.rr_port_q(aa,bb,vv,setting),native.rebuilt_q(4,aa,bb,vv,setting)]
            assert all(math.isfinite(v) for v in values) and max(values)-min(values)<1e-10,values
            count+=1
    # Four-method/native consistency at finite admissible settings.
    for method,point,setting in [('se',(2.,2.,1.),20),('af',(12.,10.,1.),10),('ni',(10.,2.,1.),7)]:
        values=[kernel(method,*point,l)[0](setting) for l in ('python','cpp','fortran')]
        assert max(values)-min(values)<1e-10,(method,values)
        count+=1
    # Compile/import the current CDF driver and exercise both native components.
    import sys
    sys.path.insert(0,str(ROOT/'experiments/section61'))
    cdf=module('cdf_driver',ROOT/'experiments/section61/run.py')
    folder=ROOT/'results/acceptance/cdf';folder.mkdir(parents=True,exist_ok=True)
    for method in ('GammaCHI','Hybrid'):
        lib,_=cdf.build(method,folder)
        bs=np.linspace(0,5,17);case={'a':2.,'nu':1.}
        branches=(bs*bs/2>3).astype(np.int32) if method=='GammaCHI' else cdf.hybrid_mask(2.,1.,bs).astype(np.int32)
        orders=np.full(len(bs),40,dtype=np.int64)
        for mode in ('scalar','batch'):
            evaluate,_=cdf.prepare(lib,mode,case,bs,branches,orders)
            v,flag=evaluate();assert np.isfinite(v).all() and (flag==0).all()
            count+=1
    print('PASS:',count,'numerical/route/component checks; source archive and Python syntax verified.')
if __name__=='__main__':main()
