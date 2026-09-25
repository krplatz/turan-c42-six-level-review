#!/usr/bin/env python3
"""Independent numerical check using hypergeometric functions and quadrature.

Requires mpmath. This is a sanity check, NOT a rigorous interval certificate.
The pair integrals use a decaying integral on [0,infinity), not the series
used by verify_certificate.py.
"""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 65
root = Path(__file__).resolve().parent
params = json.loads((root/'parameters.json').read_text())
cert = json.loads((root/'verification.json').read_text())

def rational(s):
    if '/' in str(s):
        a,b=str(s).split('/')
        return mp.mpf(a)/mp.mpf(b)
    return mp.mpf(s)

def point(v):
    return mp.mpc(rational(v['real']), rational(v['imag']))

def I(alpha,x):
    return mp.power(x,alpha)/alpha*mp.hyp2f1(1,alpha,alpha+1,x)

def Q(alpha,t,u):
    L=1-t-u
    if L<=0: return mp.mpc(0)
    def f(q):
        v=L*mp.exp(-q)
        return mp.exp(-alpha*q)/(1-v)*mp.log((1-t-v)*(1-u-v)/(t*u))
    return mp.power(L,alpha)*mp.quad(f,[0,1,4,16,64,mp.inf])

def inside(x,b):
    den=mp.mpf(b['denominator'])
    return mp.mpf(b['lo_numerator'])/den <= x <= mp.mpf(b['hi_numerator'])/den

def check_complex(x,b,name):
    if not inside(mp.re(x),b['real']) or not inside(mp.im(x),b['imag']):
        raise ValueError(name+' is outside the rigorous enclosure')

s=point(params['s']); alpha=1-s; C=rational(params['C'])
tau=rational(params['tau']); levels=[point(p) for p in params['levels']]
thresholds=[rational(t) for t in params['thresholds']]
jumps=[]; prev=s
for z in levels:
    jumps.append(z-prev); prev=z
K=I(alpha,tau); D=mp.re(I(mp.re(alpha),tau))
Y=1+s*K
check_complex(K,cert['K'],'K')
if not inside(D,cert['D']): raise ValueError('D outside interval')
for j,t in enumerate(thresholds):
    v=I(alpha,1-t)-K
    check_complex(v,cert['linear'][j],f'linear[{j}]')
    Y-=jumps[j]*v
count=0
for j,t in enumerate(thresholds):
    for k in range(j,len(thresholds)):
        q=Q(alpha,t,thresholds[k])
        check_complex(q,cert['Q'][j][k],f'Q[{j},{k}]')
        Y+=(mp.mpf('0.5') if j==k else 1)*jumps[j]*jumps[k]*q
        count+=1
check_complex(Y,cert['Y'],'Y')
ratio=abs(Y)/D
print('Independent high-precision check: PASS')
print(f'{count} pair-kernel entries checked against exact enclosures.')
print('I: hypergeometric evaluation; Q: transformed adaptive quadrature.')
for name,value in [('Re(Y)',mp.re(Y)),('Im(Y)',mp.im(Y)),('D',D),('|Y|/D',ratio),('C-|Y|/D',C-ratio),('C^2 D^2-|Y|^2',C*C*D*D-abs(Y)**2)]:
    print(name+' = '+mp.nstr(value,55))
print('Supplementary numerical check only; rigor comes from verify_certificate.py and proof.md.')
(root/'independent_check.json').write_text(json.dumps({'status':'PASS','mpmath_digits':mp.mp.dps,'pair_entries':count,'ratio':mp.nstr(ratio,60),'slack':mp.nstr(C-ratio,60),'method':'hypergeometric I and transformed Q quadrature; numerical only'},indent=2)+'\n')
