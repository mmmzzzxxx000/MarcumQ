# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Build isolated kernels; archive source remains byte-for-byte unchanged.

Generated Fortran edits only expose private components and external counts.
setting=-1 preserves each original numerical stopping condition.
Nonnegative setting chooses a fixed truncation/refinement count externally;
the recurrence, summation, integrand and special-function source stay intact.
All substitutions and hashes are saved in each build manifest.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from source939 import ROOT,UPSTREAM,verify_sources


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def replace_once(text,old,new,changes):
    if text.count(old)!=1:
        raise RuntimeError('Source patch anchor mismatch: '+repr(old))
    changes.append(dict(original=old,replacement=new))
    return text.replace(old,new,1)


def fortran_text(text):
    edits=[]
    patch=lambda old,new:None
    pairs=[
      ('PUBLIC  :: marcum','PUBLIC :: marcum,qser,pser,pqasyxy,MarcumPQtrap,bench_setting\n  INTEGER :: bench_setting=-1'),
      ('DO WHILE ((q0/q>delta).AND.(n<1000))','DO WHILE (((bench_setting<0).AND.(q0/q>delta).AND.(n<1000)).OR. &\n &((bench_setting>=0).AND.(n<bench_setting)))'),
      ('DO WHILE ((k<10000).AND.(m==0))','DO WHILE (((bench_setting<0).AND.(k<10000).AND.(m==0)).OR. &\n &((bench_setting>=0).AND.(k<=bench_setting)))'),
      ('nnmax=startingpser(mu,x,y);','IF (bench_setting>=0) THEN\n      nnmax=bench_setting\n    ELSE\n      nnmax=startingpser(mu,x,y)\n    ENDIF'),
      ('k=startingpser(mu,x,y)+1;','k=nnmax+1;'),
      ('an=sqxi; n=0; n0=100;','an=sqxi; n=0; n0=100;\n      IF (bench_setting>=0) n0=bench_setting'),
      ('DO WHILE ((abs(bn(n))>delta).AND.(n<n0))','DO WHILE ((((bench_setting<0).AND.(abs(bn(n))>delta)).OR.(bench_setting>=0)).AND.(n<n0))'),
      ('DO WHILE(((v>e).AND.(nc<10)).OR.(nc<=2))','DO WHILE (((bench_setting<0).AND.(((v>e).AND.(nc<10)).OR.(nc<=2))).OR. &\n &((bench_setting>=0).AND.(nc<bench_setting)))'),
    ]
    for old,new in pairs:text=replace_once(text,old,new,edits)
    return text,edits


def cpp_text(text):
    edits=[]
    pairs=[
      ('bool source; };','bool source; int64_t fixed=-1; };'),
      ('int starting_p(const Context&c){','int starting_p(const Context&c){\n if(c.fixed>=0)return int(c.fixed);'),
      ('while(q!=0&&q0/q>c.delta&&n<limit)','while(c.fixed>=0?n<c.fixed:(q!=0&&q0/q>c.delta&&n<limit))'),
      ('for(;k<std::min<int64_t>(10000,c.limit);)','for(;k<(c.fixed>=0?c.fixed+1:std::min<int64_t>(10000,c.limit));)'),
      ('if((sum==0&&k>150)||(sum>0&&t/sum<c.delta&&k>10))','if(c.fixed<0&&((sum==0&&k>150)||(sum>0&&t/sum<1e-16&&k>10)))'),
      ('while(std::abs(bn.back())>c.delta&&int(bn.size())-1<limit)','while(c.fixed>=0?int(bn.size())-1<c.fixed:(std::abs(bn.back())>c.delta&&int(bn.size())-1<limit))'),
      ('while((v>eps&&nc<c.refinements)||nc<=2)','while(c.fixed>=0?nc<c.fixed:((v>eps&&nc<c.refinements)||nc<=2))'),
    ]
    for old,new in pairs:text=replace_once(text,old,new,edits)
    # Avoid timing legacy result diagnostics; this adapter needs only p and q.
    start=text.index('void finish(');end=text.index('std::pair<double,double> gamma_pair',start)
    old=text[start:end]
    new='void finish(components_result &r,double p,double q,int status,int target){r.p=p;r.q=q;r.value=target==0?p:q;r.status=status;}\n'
    text=replace_once(text,old,new,edits)
    text+='\n#include "rebuild_adapter.inc"\n'
    return text,edits


def build():
    verified=verify_sources()
    src=UPSTREAM/'939/Fortran90/Src'
    for lang,compiler_name in [('cpp','g++'),('fortran','gfortran')]:
        compiler=shutil.which(compiler_name)
        if not compiler:raise RuntimeError('Compiler missing: '+compiler_name)
        folder=ROOT/'build/rebuilt'/lang;folder.mkdir(parents=True,exist_ok=True)
        commands=[]
        def run(args):
            cmd=[str(a) for a in args];commands.append(cmd)
            r=subprocess.run(cmd,cwd=str(folder),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
            if r.returncode:raise RuntimeError(r.stdout)
        suffix='.dll' if os.name=='nt' else '.so'
        library=folder/('rebuilt'+suffix)
        flags=['-O3','-fno-fast-math','-ffp-contract=off','-fPIC']
        if lang=='cpp':
            original=ROOT/'src/cpp/algorithm939.cpp'
            text,edits=cpp_text(original.read_text(encoding='utf-8'))
            generated=folder/'components_generated.cpp';generated.write_text(text,encoding='utf-8')
            flags+=['-std=c++17','-I',ROOT/'src/cpp','-I',ROOT/'src/cpp']
            run([compiler]+flags+['-shared',generated,ROOT/'src/cpp/serr.cpp','-static-libgcc','-static-libstdc++','-o',library])
            provenance={str(original.relative_to(ROOT)):digest(original)}
        else:
            generated=folder/'marcumq_generated.f90'
            text,edits=fortran_text((src/'marcumq.f90').read_text(encoding='utf-8'))
            generated.write_text(text,encoding='utf-8')
            flags+=['-frecursive','-ffree-line-length-none','-I',src,'-J',folder]
            objects=[]
            for source in (src/'cons.f90',src/'incgam.f90',generated,ROOT/'src/fortran/serr.f90',ROOT/'src/fortran/rebuild_adapter.f90'):
                obj=folder/(source.stem+'.o');run([compiler]+flags+['-c',source,'-o',obj]);objects.append(obj)
            run([compiler,'-shared']+objects+['-static-libgfortran','-static-libgcc','-o',library])
            provenance=verified
        serr_source=ROOT/('src/cpp/serr.cpp' if lang=='cpp' else 'src/fortran/serr.f90')
        provenance[str(serr_source.relative_to(ROOT))]=digest(serr_source)
        record=dict(language=lang,commands=commands,source_hashes=provenance,
                    changes=edits,generated_sha256=digest(generated),library_sha256=digest(library),
                    note='Original stopping at setting=-1; external fixed counts otherwise; unchanged upstream archive.')
        (folder/'manifest.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
        print(lang+': '+str(library),flush=True)


if __name__=='__main__':build()
