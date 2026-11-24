# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Neutral Table-1 oracle: defining Bessel integral, not either tested series.

This is arbitrary-precision numerical quadrature, not a rigorous interval
certificate. Agreement under precision and partition changes is recorded.
"""
import mpmath as mp


def q_integral(a,b,nu,*,dps,partition):
    with mp.workdps(dps):
        a,b,nu=(mp.mpf(float(v)) for v in (a,b,nu))
        if not (a>0 and b>=0 and nu>0):
            raise ValueError('Defining integral requires a>0, b>=0, nu>0')
        factor=a**(1-nu)
        def f(t):
            if t==0:return mp.mpf(0)
            return factor*t**nu*mp.exp(-(t*t+a*a)/2)*mp.besseli(nu-1,a*t)
        points=[b]+[b+mp.mpf(v) for v in partition]+[mp.inf]
        q,error=mp.quad(f,points,method='tanh-sinh',error=True)
        return +q,+error


def table1_reference(a,b,nu):
    low,err_low=q_integral(a,b,nu,dps=90,partition=(2,8))
    high,err_high=q_integral(a,b,nu,dps=160,partition=(1,4,12))
    with mp.workdps(170):
        difference=abs(low-high)/abs(high)
        if difference>mp.mpf('1e-70'):
            raise ArithmeticError('Integral references disagree under precision/partition change')
        return dict(q=mp.nstr(high,150),reference_method='defining Bessel integral, mpmath tanh-sinh',
                    dps=[90,160],partitions=[[2,8],[1,4,12]],
                    precision_check_relative_difference=mp.nstr(difference,20),
                    quadrature_estimated_errors=[mp.nstr(err_low,15),mp.nstr(err_high,15)],
                    interval_certificate=False)


def _scaled_integral(a,b,nu,*,dps,variant):
    """Same defining integral, scaled for tiny tails and partitioned at peak.

    Scaling and extra breakpoints belong only to the high-precision oracle.
    They never enter a timed method. Multiplying by exp(log_scale) restores
    the original integral exactly in real arithmetic.
    """
    with mp.workdps(dps):
        a,b,nu=(mp.mpf(float(v)) for v in (a,b,nu))
        if not (a>0 and b>0 and nu>0):raise ValueError('Positive arguments required')
        centre=mp.sqrt(a*a+2*nu)
        anchor=max(b,centre)
        def log_density(t):
            return (1-nu)*mp.log(a)+nu*mp.log(t)-(t*t+a*a)/2+mp.log(mp.besseli(nu-1,a*t))
        log_scale=log_density(anchor)
        offsets=(2,8) if variant==0 else (1,4,12)
        around=(-8,0,8) if variant==0 else (-4,0,4,12)
        finite=sorted(set([b]+[b+v for v in offsets]+[centre+v for v in around if centre+v>b]))
        normalized,error=mp.quad(lambda t:mp.exp(log_density(t)-log_scale),finite+[mp.inf],
                                 method='tanh-sinh',error=True)
        scale=mp.exp(log_scale)
        return +(normalized*scale),+(error*scale),[mp.nstr(v,20) for v in finite]


def region_reference(a,b,nu):
    low,el,pl=_scaled_integral(a,b,nu,dps=90,variant=0)
    high,eh,ph=_scaled_integral(a,b,nu,dps=160,variant=1)
    with mp.workdps(170):
        if not (high>0):raise ArithmeticError('Integral returned nonpositive Q')
        difference=abs(low-high)/abs(high)
        if difference>mp.mpf('1e-65'):
            raise ArithmeticError('Scaled integral references disagree')
        return dict(q=mp.nstr(high,150),reference_method='defining Bessel integral, scaled mpmath tanh-sinh',
                    dps=[90,160],finite_breakpoints=[pl,ph],
                    precision_check_relative_difference=mp.nstr(difference,20),
                    quadrature_estimated_errors=[mp.nstr(el,15),mp.nstr(eh,15)],
                    interval_certificate=False)
