# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Algorithm 939's prec, qrec and fc in Python; reuse its quadrature port."""
import math
import sys
import time

from marcum import algorithm939 as upstream


from marcum.algorithm939 import *
from marcum.algorithm939 import special, zeta, xminsin
def seed(x,y,nu,setting,want_p):
    xs=x/nu;ys=y/nu
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
    return (1-pq if zz<0 else -pq) if want_p else (pq if zz<0 else 1+pq)


def fc(order, z):
    tiny = sys.float_info.min * 10
    ratio, c, d, change = tiny, tiny, 0.0, 0.0
    m = 0
    b = 2 * order / z
    while abs(change-1) > 1e-15:
        d = b+d
        if abs(d) < tiny:
            d = tiny
        c = b+1/c
        if abs(c) < tiny:
            c = tiny
        d = 1/d
        change = c*d
        ratio *= change
        m += 1
        b = 2*(order+m)/z
    return ratio


def prec(x, y, order, level):
    limit = y-x+1+math.sqrt(2*(x+y)+1)
    steps = int(limit)+2-int(order)
    shifted = order+steps
    c = math.sqrt(y/x)*fc(shifted, 2*math.sqrt(x*y))
    p1 = seed(x, y, shifted+1, level, True)
    p0 = seed(x, y, shifted, level, True)
    for n in range(steps):
        p = ((1+c)*p0-p1)/c
        p1, p0 = p0, p
        c = y/(shifted-n-1+x*c)
    return 1-p0


def rr(a, b, order, setting):
    x, y = a*a/2, b*b/2
    if y <= x+order:
        return prec(x, y, order, setting)
    limit = y-x+(1-math.sqrt(2*(x+y)+1))
    if limit < 5:
        if x < 200:
            return upstream.se(a, b, order, setting)
        return prec(x, y, order, setting)
    steps = int(order)-(int(limit)-1)
    shifted = order-steps
    c = [math.sqrt(y/x)*fc(order, 2*math.sqrt(x*y))]
    for n in range(1, steps+1):
        c.append(y/(order-n+x*c[-1]))
    q0 = seed(x, y, shifted-1, setting, False)
    q1 = seed(x, y, shifted, setting, False)
    for n in range(steps, 0, -1):
        q0, q1 = q1, (1+c[n])*q1-c[n]*q0
    return q1


def kernel(a, b, order):
    def evaluate(setting):
        return rr(a, b, order, setting)

    def timer(setting, repeats):
        start = time.perf_counter()
        for _ in range(repeats):
            rr(a, b, order, setting)
        return (time.perf_counter()-start)/repeats

    return evaluate, timer
