#!/usr/bin/env python3
"""Exact algebra checks and negative controls for the C_42 work package.

Finite tests do not replace the general proof. They test the exact
coefficient identity, outward rounding, and rejection of invalid data.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import random
import tempfile
from verify_certificate import IV, SCALE, verify, require

root=Path(__file__).resolve().parent
p=json.loads((root/'parameters.json').read_text())

def c(re=0,im=0): return (F(re),F(im))
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def neg(x): return (-x[0],-x[1])
def sub(x,y): return add(x,neg(y))
def mul(x,y): return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def scale(x,q): return (x[0]*q,x[1]*q)
def total(xs):
    z=c()
    for x in xs: z=add(z,x)
    return z

def point(v): return c(v['real'],v['imag'])
s=point(p['s']); alpha=sub(c(1),s)
levels=[point(v) for v in p['levels']]
thresholds=[F(t) for t in p['thresholds']]

def algebra_check(n):
    A=n//3; J=range(A+1,n-A); final=range(n-A,n+1)
    S=[c()]*(n+1)
    for k in range(1,A+1): S[k]=s
    for m in J:
        for j,t in enumerate(thresholds):
            if m > (t.numerator*n)//t.denominator:
                S[m]=levels[j]
    beta=[c(1)]
    for l in range(1,n+1):
        beta.append(scale(mul(beta[-1],add(alpha,c(l-1))),F(1,l)))
    g={m:scale(sub(S[m],s),F(1,m)) for m in J}
    linear=total(mul(g[m],beta[n-m]) for m in J)
    quadratic=total(mul(mul(g[m],g[k]),beta[n-m-k])
                    for m in J for k in J if m+k<=n)
    P={m:scale(beta[n-m],F(n,m)) for m in final}
    Y=add(scale(add(sub(beta[n],linear),scale(quadratic,F(1,2))),n),
          mul(s,total(P.values())))
    # An algebra-only closing choice. This intentionally does not impose a
    # radius bound at these small n. P_n=1, so S_n=Y closes the equation.
    for m in final: S[m]=c()
    S[n]=Y
    b=[c(1)]
    for k in range(1,n+1):
        rhs=sub(total(b),total(mul(S[m],b[k-m]) for m in range(1,k+1)))
        b.append(scale(rhs,F(1,k)))
    require(b[n]==c(),f'exact nth coefficient failed at n={n}')
    # Recover the power sums of the monic degree-n polynomial
    # (Z-1)(Z^(n-1)+b1 Z^(n-2)+...+b_(n-1)) using Newton identities.
    coeff=[c(1)]+[sub(b[k],b[k-1]) for k in range(1,n+1)]
    moments=[c(n)]
    for k in range(1,n+1):
        value=neg(add(total(mul(coeff[j],moments[k-j]) for j in range(1,k)),
                      scale(coeff[k],k)))
        moments.append(value)
        require(value==S[k],f'exact Newton identity failed at n={n}, k={k}')

for n in range(3,37): algebra_check(n)
print('Exact coefficient closure and Newton identities: PASS for n=3,...,36')

rng=random.Random(20260925)
for _ in range(1000):
    a,b=F(rng.randrange(-10**8,10**8),rng.randrange(1,1000)),F(rng.randrange(-10**8,10**8),rng.randrange(1,1000))
    x,y=IV.exact(a),IV.exact(b)
    checks=[(x+y,a+b),(x-y,a-b),(x*y,a*b)]
    if b: checks.append((x/b,a/b))
    for interval,truth in checks:
        require(F(interval.lo,SCALE)<=truth<=F(interval.hi,SCALE),'outward-rounding check failed')
print('Outward-rounded rational arithmetic: PASS on 1,000 deterministic test pairs')

with tempfile.TemporaryDirectory() as td:
    path=Path(td)/'bad.json'
    bad=dict(p); bad['C']='0.68'; path.write_text(json.dumps(bad))
    try: verify(path)
    except ValueError as e:
        require('profile exceeds' in str(e),'unexpected negative-control failure')
    else: raise ValueError('invalid radius was accepted')
    bad=dict(p); bad['s']={'real':'0','imag':'0'}
    bad['levels']=[{'real':'0','imag':'0'} for _ in p['levels']]
    path.write_text(json.dumps(bad))
    try: verify(path)
    except ValueError as e:
        require('crucial' in str(e),'unexpected zero-profile failure')
    else: raise ValueError('invalid limiting inequality was accepted')
print('Negative controls: invalid radius and invalid limiting inequality both rejected')
print('All audit checks: PASS')
