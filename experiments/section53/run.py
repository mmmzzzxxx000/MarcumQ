# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Rebuild and run the current Fortran random-point experiment."""
import argparse, os, shutil, subprocess
from pathlib import Path
from marcum._paths import ROOT
from source939 import UPSTREAM,verify_sources
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--count',type=int,default=100000000)
parser.add_argument('--threads',type=int,default=4)
args=parser.parse_args()
if args.count<1 or args.threads<1: parser.error('count and threads must be positive')
verify_sources()
folder=ROOT/'results/section53';folder.mkdir(parents=True,exist_ok=True)
src=UPSTREAM/'939/Fortran90/Src'
compiler=shutil.which('gfortran')
if not compiler: raise RuntimeError('gfortran is required')
exe=folder/'validate.exe'
command=[compiler,'-O3','-fno-fast-math','-ffp-contract=off','-frecursive','-ffree-line-length-none','-fopenmp','-static','-I',str(src),'-J',str(folder)]
command += [str(src/name) for name in ('cons.f90','incgam.f90','marcumq.f90')]
command += [str(ROOT/'src/fortran/serr.f90'),str(Path(__file__).with_name('validate.f90')),'-lgomp','-ldl','-o',str(exe)]
subprocess.run(command,cwd=folder,check=True)
env=dict(os.environ,OMP_NUM_THREADS=str(args.threads))
with (folder/'run.txt').open('w',encoding='utf-8') as handle:
    proc=subprocess.Popen([str(exe),str(args.count)],cwd=folder,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    for line in proc.stdout: print(line,end='',flush=True);handle.write(line);handle.flush()
    if proc.wait(): raise RuntimeError('Random experiment failed; inspect results/section53/run.txt')
