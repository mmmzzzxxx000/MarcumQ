# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Check all public language wrappers against shared routes and Bessel integrals."""
import ctypes
import json
from pathlib import Path
import numpy as np
import mpmath as mp
from marcum._paths import ROOT
from marcum.hybrid import METHODS, _library, _cpp_library, select_method, hybrid_q
from integral_reference import q_integral


def main():
    fortran, cpp = _library(), _cpp_library()
    signature = [ctypes.c_double]*3+[ctypes.c_int]
    selectors = [fortran.marcum_fortran_select_method, cpp.marcum_cpp_select_method]
    for function in selectors:
        function.argtypes=signature
        function.restype=ctypes.c_int
    # Include grid boundaries, noninteger orders, and large-order transitions.
    rng=np.random.RandomState(20261120)
    random=np.exp(rng.uniform(np.log([.01,.001,.1]),np.log([100.,100.,1000.]),size=(5000,3)))
    grid=np.array([(a,b,nu) for nu in (1.,10.,20.,135.,200.,1000.)
                   for a in (1.,3.,7.,8.,10.,20.,100.)
                   for b in (1.,3.,7.,8.,10.,20.,100.)])
    points=np.concatenate([random,grid])
    for a,b,nu in points:
        for proposed in (False,True):
            expected=METHODS.index(select_method(a,b,nu,proposed))
            assert all(fn(a,b,nu,int(proposed))==expected for fn in selectors),(a,b,nu,proposed,expected)

    cases=[('SE-RR',8.,2.,1.),('SE',1.,7.,1.),('AF',12.,10.,1.),
           ('NI',20.,8.,200.),('RR',8.,10.,20.),('LAE',10.,22.,200.)]
    rows=[]
    for method,a,b,nu in cases:
        assert select_method(a,b,nu)==method
        low,_=q_integral(a,b,nu,dps=40,partition=(2,8))
        high,_=q_integral(a,b,nu,dps=55,partition=(1,4,12))
        with mp.workdps(60):
            assert abs(low-high)<=mp.mpf('1e-30')*high
            for language in ('python','cpp','fortran'):
                settings=(20,40,80,128,256) if method in ('SE','SE-RR') else (2,4,6,8,10,12,16)
                result=None
                for setting in settings:
                    value=hybrid_q(a,b,nu,setting,language=language)
                    error=abs(mp.mpf(value)-high)/high
                    if mp.isfinite(error) and error<=mp.mpf('1e-10'):
                        result=dict(method=method,language=language,a=a,b=b,nu=nu,setting=setting,
                                    q=value,relative_error=float(error),reference=mp.nstr(high,55))
                        break
                assert result is not None,(method,language,value,str(error))
                rows.append(result)
    out=ROOT/'results/acceptance'
    out.mkdir(parents=True,exist_ok=True)
    (out/'hybrid_interfaces.json').write_text(json.dumps(dict(route_checks=2*len(points),
        evaluations=rows,shared_lae='Python/C++ call the official Fortran LAE component.'),indent=2),encoding='utf-8')
    print('PASS:',2*len(points),'route triples per native selector and',len(rows),'six-method language evaluations.')


if __name__=='__main__': main()
