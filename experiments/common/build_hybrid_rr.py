# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Build the unchanged 939 RR routines with fixed-count seed adapters."""
import json, shutil, subprocess
from source939 import ROOT, UPSTREAM, verify_sources
from build_rebuilt import fortran_text, digest
def build():
    verified=verify_sources()
    folder=ROOT/'build/hybrid_timing_fortran_rr'
    folder.mkdir(parents=True,exist_ok=True)
    src=UPSTREAM/'939/Fortran90/Src'
    original=(src/'marcumq.f90').read_text(encoding='utf-8')
    text, changes=fortran_text(original)
    anchor='PUBLIC :: marcum,qser,pser,pqasyxy,MarcumPQtrap,bench_setting'
    assert text.count(anchor)==1
    text=text.replace(anchor,anchor+',prec,qrec')
    for name in ('prec','qrec'):
        first=original.index('    SUBROUTINE '+name+'(')
        last=original.index('    END SUBROUTINE '+name,first)+len('    END SUBROUTINE '+name)
        assert original[first:last] in text
    generated=folder/'marcumq_generated.f90';generated.write_text(text,encoding='utf-8')
    flags=['-O3','-fno-fast-math','-ffp-contract=off','-fPIC','-frecursive','-ffree-line-length-none','-I',str(src),'-J',str(folder)]
    compiler=shutil.which('gfortran')
    if not compiler: raise RuntimeError('Install gfortran and add it to PATH.')
    objects=[];commands=[]
    for source in [src/'cons.f90',src/'incgam.f90',generated,ROOT/'src/fortran/serr.f90',ROOT/'src/fortran/rebuild_rr_adapter.f90']:
        obj=folder/(source.stem+'.o')
        command=[compiler]+flags+['-c',str(source),'-o',str(obj)]
        subprocess.run(command,cwd=folder,check=True);commands.append(command);objects.append(str(obj))
    library=folder/'hybrid_timing.dll'
    command=[compiler,'-shared']+objects+['-static-libgfortran','-static-libgcc','-o',str(library)]
    subprocess.run(command,cwd=folder,check=True);commands.append(command)
    (folder/'manifest.json').write_text(json.dumps(dict(library_sha256=digest(library),upstream=verified,changes=changes,commands=commands),indent=2),encoding='utf-8')
if __name__=='__main__':build()
