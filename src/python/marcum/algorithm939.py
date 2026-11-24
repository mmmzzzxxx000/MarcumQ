# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Raw binary64 SE/AF/NI ports for the rebuilt, same-accuracy experiment.

Source: Algorithm 939 marcumq.f90, qser/pser, pqasyxy, MarcumPQtrap.
The default setting=-1 retains source stopping. A nonnegative setting
exposes ONLY truncation: last series index, or trapezoidal refinement level.
No reference, new error test, recovery, clipping or SERR kernel occurs here.
Python uses SciPy incomplete-gamma ratios and erfcx as library primitives.
"""
import math
import sys
from scipy import special
from .kernels.series import _starting_p
from .kernels.asymptotic import _ignega
from .kernels.quadrature import zeta, xminsin


def se(a,b,nu,setting=-1):
    x,y=a*a/2,b*b/2
    if y>x+nu:
        q0=q=float(special.gammaincc(nu,y))
        lh=nu*math.log(y)-y-float(special.gammaln(nu+1))
        if lh>math.log(sys.float_info.min*10) and x<100:
            h=math.exp(lh);n=0
            while (n<setting if setting>=0 else q0/q>1e-17 and n<1000):
                q0=x*(q0+h)/(n+1.)
                h=x*y*h/((n+1.)*(nu+n+1))
                q+=q0;n+=1
            return math.exp(-x)*q
        total=0.;k=0
        while k<(setting+1 if setting>=0 else 10000):
            t=math.exp(-x+k*math.log(x)-float(special.gammaln(k+1)))*float(special.gammaincc(nu+k,y))
            total+=t;k+=1
            if setting<0 and ((total==0 and k>150) or (total>0 and t/total<1e-16 and k>10)):
                break
        return total
    N=setting if setting>=0 else _starting_p(nu,x,y,1e-15,1000000)
    n=N+1
    lh=-x-y+n*math.log(x)+(n+nu)*math.log(y)-float(special.gammaln(nu+n+1))-float(special.gammaln(n+1))
    if lh<math.log(sys.float_info.min*10):
        total=0.
        for k in range(n,0,-1):
            total+=math.exp(k*math.log(x)-float(special.gammaln(k+1)))*float(special.gammainc(nu+k,y))
        p=(total+float(special.gammainc(nu,y)))*math.exp(-x)
    else:
        h=math.exp(lh)
        p1=float(special.gammainc(nu+n,y))*math.exp(-x+n*math.log(x)-float(special.gammaln(n+1)))
        p=0.
        while n>0:
            h=h*n*(nu+n)/(x*y)
            p1=n*p1/x+h
            p+=p1;n-=1
    return 1-p


def af(a,b,nu,setting=-1):
    x,y=a*a/2,b*b/2
    xi=a*b;rho=b/a
    z=(y-x)*(y-x)/(x+y+xi)
    rhomu=math.exp(nu*math.log(rho))
    s=1. if y>=x else -1.
    psi0=.5*rhomu*float(special.erfcx(math.sqrt(z)))/math.sqrt(rho)
    v=2*nu-1;rhom=v*(rho-1);rhop=2*(rho+1)
    c=s*rhomu/math.sqrt(8*math.pi);an=math.sqrt(xi);bn=[1.]
    while (len(bn)-1<setting if setting>=0 else abs(bn[-1])>1e-17 and len(bn)-1<100):
        n=len(bn)
        an=(4*nu*nu-(2*n-1)**2)*an/(8*n*xi)
        bn.append(an*(rhom-n*rhop)/(rho*(v+2*n)))
    N=len(bn)-1
    if not N:raise ValueError('pqasyxy requires a nonzero last index')
    nr=min(int(z)+1,N);phi=[0.]*(N+1)
    phi[nr]=math.exp((nr-.5)*math.log(z))*_ignega(nr,z,1e-17)
    for n in range(nr+1,N+1):phi[n]=(-z*phi[n-1]+1)/(n-.5)
    for n in range(nr-1,0,-1):phi[n]=(1-(n+.5)*phi[n+1])/z
    pq=psi0
    for n in range(1,N+1):
        c=-c;pq+=c*bn[n]*phi[n]
    pq*=math.exp(-z)
    return pq if s==1 else 1-pq


def ni(a,b,nu,setting=-1):
    xs=a*a/(2*nu);ys=b*b/(2*nu)
    xis2=4*xs*ys;wxis=math.sqrt(1+xis2)
    def trapsum(h,d):
        b0=3.
        def f(theta):
            nonlocal b0
            if theta>b0:return 0.
            if abs(theta)<1e-10:
                r=(1+wxis)/(2*ys)
                return r/(1-r)*math.exp(-nu*wxis*theta*theta*.5)
            t2=theta*theta;st=math.sin(theta);ct=math.cos(theta)
            ts=theta/st;s2=st*st;wx=math.sqrt(ts*ts+xis2);xm=xminsin(theta)
            p=xm*t2*ts/6
            term1=(p*(ts+1)-t2-s2*xis2)/(ct*wx+wxis)
            p=p*(1+(ts+1)/(wx+wxis))/(1+wxis)
            exponent=nu*(term1-math.log1p(p))
            if exponent<=math.log(1e-16):
                b0=min(theta,b0);factor=0.
            else:factor=math.exp(exponent)
            r=(ts+wx)/(2*ys)
            p=(2*theta*math.sin(theta/2)**2-xm*t2*theta/6)/(2*ys*s2)
            dr=p*(1+ts/wx)
            return factor*(dr*st+(ct-r)*r)/(r*(r-2*ct)+1)
        total=f(0)/2 if d==0 else 0.
        theta=h if d==0 else d
        end=3-h/2 if d==0 else 3.
        while theta<end and theta<b0:
            total+=f(theta);theta+=h
        return total*h
    h=3/8;pq=trapsum(h,0);v=1.;nc=0
    while (nc<setting if setting>=0 else (v>1e-13 and nc<10) or nc<=2):
        nc+=1;q=trapsum(h,h/2)
        v=abs(pq/q-1) if abs(q)>0 else 0.
        h/=2;pq=(pq+q)/2
    zz=zeta(xs,ys)
    pq=pq*math.exp(-nu*.5*zz*zz)/math.pi
    return pq if zz<0 else 1+pq
