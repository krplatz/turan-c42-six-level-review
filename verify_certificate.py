#!/usr/bin/env python3
"""Exact fixed-point interval certificate for the Turan C_42 profile.

No third-party packages or floating-point pass/fail comparisons are used.
This checks the finite analytic criterion in proof.md, not its asymptotic
or algebraic proof. It is not a Lean kernel certificate.

Usage: python3 verify_certificate.py [parameters.json] [verification.json]
All decimal input strings denote exact rational numbers.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
import sys

DIGITS = 65
SCALE = 10 ** DIGITS


def require(condition: bool, message: str) -> None:
    """Unlike assert, these checks remain active under python -O."""
    if not condition:
        raise ValueError(message)


def ceildiv(a: int, b: int) -> int:
    require(b > 0, 'ceildiv requires a positive denominator')
    return -((-a) // b)


@dataclass(frozen=True)
class IV:
    lo: int
    hi: int

    def __post_init__(self) -> None:
        require(self.lo <= self.hi, 'reversed interval')

    @staticmethod
    def exact(x: F | int | str) -> IV:
        x = F(x)
        return IV((x.numerator * SCALE) // x.denominator,
                  ceildiv(x.numerator * SCALE, x.denominator))

    @staticmethod
    def bounds(lo: F, hi: F) -> IV:
        return IV(IV.exact(lo).lo, IV.exact(hi).hi)

    def __add__(self, other: IV | F | int) -> IV:
        b = as_iv(other)
        return IV(self.lo + b.lo, self.hi + b.hi)

    __radd__ = __add__

    def __neg__(self) -> IV:
        return IV(-self.hi, -self.lo)

    def __sub__(self, other: IV | F | int) -> IV:
        return self + (-as_iv(other))

    def __rsub__(self, other: IV | F | int) -> IV:
        return as_iv(other) + (-self)

    def __mul__(self, other: IV | F | int) -> IV:
        b = as_iv(other)
        products = (self.lo*b.lo, self.lo*b.hi,
                    self.hi*b.lo, self.hi*b.hi)
        return IV(min(products)//SCALE, ceildiv(max(products), SCALE))

    __rmul__ = __mul__

    def __truediv__(self, other: F | int) -> IV:
        x = F(other)
        require(x != 0, 'division by zero')
        return self * (1 / x)

    def abs_upper(self) -> F:
        return F(max(abs(self.lo), abs(self.hi)), SCALE)

    def widen(self, error: F) -> IV:
        require(error >= 0, 'negative error radius')
        e = IV.exact(error).hi
        return IV(self.lo-e, self.hi+e)

    def payload(self) -> dict[str, str]:
        return {'lo_numerator': str(self.lo), 'hi_numerator': str(self.hi),
                'denominator': str(SCALE),
                'lo_decimal': exact_decimal(self.lo),
                'hi_decimal': exact_decimal(self.hi)}


def as_iv(x: IV | F | int) -> IV:
    return x if isinstance(x, IV) else IV.exact(x)


def exact_decimal(n: int) -> str:
    sign = '-' if n < 0 else ''
    a, b = divmod(abs(n), SCALE)
    return f'{sign}{a}.{b:0{DIGITS}d}'


@dataclass(frozen=True)
class CV:
    re: IV
    im: IV

    @staticmethod
    def exact(re: F | int, im: F | int = 0) -> CV:
        return CV(IV.exact(re), IV.exact(im))

    def __add__(self, b: CV) -> CV:
        return CV(self.re+b.re, self.im+b.im)

    def __neg__(self) -> CV:
        return CV(-self.re, -self.im)

    def __sub__(self, b: CV) -> CV:
        return self + (-b)

    def __mul__(self, b: CV) -> CV:
        return CV(self.re*b.re-self.im*b.im, self.re*b.im+self.im*b.re)

    def scale(self, x: IV | F | int) -> CV:
        return CV(self.re*x, self.im*x)

    def widen(self, error: F) -> CV:
        # A disk of this radius is contained in this axis-aligned rectangle.
        return CV(self.re.widen(error), self.im.widen(error))

    def norm_squared_upper(self) -> F:
        return self.re.abs_upper()**2 + self.im.abs_upper()**2

    def payload(self) -> dict:
        return {'real': self.re.payload(), 'imag': self.im.payload()}


@lru_cache(maxsize=None)
def log_unit(q: F) -> IV:
    """log(q) for exact 1 <= q <= 2 via its atanh series."""
    require(1 <= q <= 2, 'log_unit domain')
    u = (q-1)/(q+1)
    N = 100
    u_iv = IV.exact(u)
    u2 = IV.exact(u*u)
    power = u_iv
    result = IV.exact(0)
    for j in range(N):
        result = result + power * F(2, 2*j+1)
        power = power * u2
    tail = 2*u**(2*N+1)/((2*N+1)*(1-u*u))
    return result + IV.bounds(F(0), tail)


@lru_cache(maxsize=None)
def log_exact(x: F) -> IV:
    require(x > 0, 'logarithm domain')
    q, k = x, 0
    while q < 1:
        q *= 2
        k -= 1
    while q >= 2:
        q /= 2
        k += 1
    return log_unit(q) + log_unit(F(2))*k


def complex_exp(z: CV) -> CV:
    """Exponential Taylor enclosure followed by repeated squaring."""
    reduced = z
    k = 0
    while reduced.re.abs_upper()+reduced.im.abs_upper() > F(1, 2):
        reduced = reduced.scale(F(1, 2))
        k += 1
    rho = reduced.re.abs_upper()+reduced.im.abs_upper()
    N = 80
    term = CV.exact(1)
    result = term
    for j in range(1, N+1):
        term = (term * reduced).scale(F(1, j))
        result = result + term
    # |exp(z)-sum_{j=0}^N z^j/j!| <= exp(rho)*rho^(N+1)/(N+1)!
    # rho <= 1/2 and exp(rho) < 2.
    tail = 2*rho**(N+1)/factorial(N+1)
    result = result.widen(tail)
    for _ in range(k):
        result = result * result
    return result


def x_power(alpha: tuple[F, F], x: F) -> CV:
    require(0 < x < 1, 'power domain')
    logarithm = log_exact(x)
    return complex_exp(CV(logarithm*alpha[0], logarithm*alpha[1]))


def reciprocal_alpha(alpha: tuple[F, F], m: int) -> CV:
    a, b = alpha[0]+m, alpha[1]
    d = a*a+b*b
    require(d > 0, 'zero complex denominator')
    return CV.exact(a/d, -b/d)


def integral_I(alpha: tuple[F, F], x: F, N: int) -> CV:
    """Enclose integral_0^x u^(alpha-1)/(1-u) du."""
    a = alpha[0]
    require(a > 0 and 0 < x < 1 and N > 0, 'I domain')
    power = IV.exact(1)
    answer = CV.exact(0)
    for m in range(N):
        answer = answer + reciprocal_alpha(alpha, m).scale(power)
        power = power*x
    answer = x_power(alpha, x)*answer
    # Using x^a <= 1 and |alpha+m| >= a+m.
    tail = x**N/((a+N)*(1-x))
    return answer.widen(tail)


def integral_Q(alpha: tuple[F, F], t: F, u: F, N: int) -> CV:
    """Enclose the pair kernel from the generalized quadratic construction."""
    L = 1-t-u
    if L <= 0:
        return CV.exact(0)
    a = alpha[0]
    require(a > 0 and 0 < t < 1 and 0 < u < 1, 'Q domain')
    B, E = 1-t, 1-u
    ratio_t, ratio_u = L/B, L/E
    require(0 < ratio_t < 1 and 0 < ratio_u < 1, 'Q series convergence')
    c0 = log_exact(B*E/(t*u))
    require(c0.lo > 0, 'Q leading coefficient must be positive')
    scaled = c0
    answer = reciprocal_alpha(alpha, 0).scale(scaled)
    pt, pu = IV.exact(1), IV.exact(1)
    for m in range(1, N):
        pt, pu = pt*ratio_t, pu*ratio_u
        scaled = scaled*L-(pt+pu)/m
        answer = answer + reciprocal_alpha(alpha, m).scale(scaled)
    answer = x_power(alpha, L)*answer
    # |c_m| <= c0 + B^(-m)/(1-B) + E^(-m)/(1-E).
    # Dropping L^a <= 1 gives this rational tail bound.
    c0_upper = F(c0.hi, SCALE)
    tail = (c0_upper*L**N/(1-L)
            + ratio_t**N/(t*(1-ratio_t))
            + ratio_u**N/(u*(1-ratio_u)))/(a+N)
    return answer.widen(tail)


def parse_complex(data: dict) -> tuple[F, F]:
    return F(data['real']), F(data['imag'])


def rational_string(x: F) -> str:
    return f'{x.numerator}/{x.denominator}'


def verify(path: Path) -> dict:
    raw = path.read_bytes()
    data = json.loads(raw)
    require(data['decimal_scale_digits'] == DIGITS, 'unsupported scale')
    N = data['series_terms']
    require(isinstance(N, int) and 100 <= N <= 1000, 'invalid series length')
    C, tau = F(data['C']), F(data['tau'])
    old = F(data['previous_proposed_bound'])
    s = parse_complex(data['s'])
    levels = [parse_complex(p) for p in data['levels']]
    thresholds = [F(t) for t in data['thresholds']]
    require(F(1, 3) <= tau < F(1, 2), 'quadratic construction needs tau >= 1/3')
    require(0 < C < old, 'claimed bound is not a positive improvement')
    require(len(levels) == len(thresholds) and len(levels) > 0, 'profile length')
    require(thresholds[0] == tau, 'first threshold must be tau')
    require(all(t < u for t, u in zip(thresholds, thresholds[1:])), 'threshold order')
    require(thresholds[-1] < 1-tau, 'empty last middle interval')
    alpha = (1-s[0], -s[1])
    require(alpha[0] > 0, 'endpoint integrability needs Re(alpha)>0')
    radii = [C*C-x*x-y*y for x, y in [s]+levels]
    require(all(r > 0 for r in radii), 'profile exceeds claimed radius')

    jumps: list[CV] = []
    previous = s
    for current in levels:
        jumps.append(CV.exact(current[0]-previous[0], current[1]-previous[1]))
        previous = current
    K = integral_I(alpha, tau, N)
    D = integral_I((alpha[0], F(0)), tau, N).re
    require(D.lo > 0, 'positive denominator not certified')
    linear = [integral_I(alpha, 1-t, N)-K for t in thresholds]
    Q: list[list[CV]] = [[CV.exact(0) for _ in levels] for _ in levels]
    Y = CV.exact(1) + CV.exact(*s)*K
    for j, w in enumerate(jumps):
        Y = Y-w*linear[j]
    for j, t in enumerate(thresholds):
        for k in range(j, len(thresholds)):
            value = integral_Q(alpha, t, thresholds[k], N)
            Q[j][k] = Q[k][j] = value
            weight = F(1, 2) if j == k else F(1)
            Y = Y + (jumps[j]*jumps[k]*value).scale(weight)

    Y2_upper = Y.norm_squared_upper()
    C2D2_lower = C*C*F(D.lo, SCALE)**2
    margin = C2D2_lower-Y2_upper
    require(margin > 0, 'the crucial |Y| < C D inequality was not certified')
    # A simpler, human-readable rational lower bound, checked exactly.
    readable_margin = F(1, 10_000_000)
    while not margin > readable_margin:
        readable_margin /= 10
    return {'status':'PASS', 'criterion':'max(|s|, |eta_j|, |Y|/D) < C',
            'C':data['C'], 'previous_proposed_bound':data['previous_proposed_bound'],
            'parameters_sha256':sha256(raw).hexdigest(),
            'series_terms':N, 'scale_digits':DIGITS,
            'alpha':{'real':rational_string(alpha[0]),'imag':rational_string(alpha[1])},
            'profile_radius_squared_slacks':[rational_string(r) for r in radii],
            'K':K.payload(),'D':D.payload(),
            'linear':[x.payload() for x in linear],
            'Q':[[x.payload() for x in row] for row in Q], 'Y':Y.payload(),
            'Y_squared_upper':rational_string(Y2_upper),
            'C_squared_D_squared_lower':rational_string(C2D2_lower),
            'squared_inequality_margin_lower':rational_string(margin),
            'simple_margin_lower':rational_string(readable_margin),
            'scope':'Finite limiting criterion only. The analytic construction is proved separately in proof.md. Not Lean-verified.'}


def main() -> None:
    root = Path(__file__).resolve().parent
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else root/'parameters.json'
    output_path = Path(sys.argv[2]) if len(sys.argv) > 2 else root/'verification.json'
    result = verify(input_path)
    output_path.write_text(json.dumps(result, indent=2)+'\n')
    print('Turan C_42 six-level interval certificate: PASS')
    print('All decimal inputs interpreted as exact rational numbers.')
    print('Radius checks: all 7 strict inequalities PASS')
    print('|Y|^2 < C^2 D^2: PASS')
    print('Squared-inequality margin > '+result['simple_margin_lower'])
    print('Certified criterion for C = '+result['C'])
    print('Re(Y) in ['+result['Y']['real']['lo_decimal']+', '+result['Y']['real']['hi_decimal']+']')
    print('Im(Y) in ['+result['Y']['imag']['lo_decimal']+', '+result['Y']['imag']['hi_decimal']+']')
    print('D in ['+result['D']['lo_decimal']+', '+result['D']['hi_decimal']+']')
    print('Parameters SHA-256: '+result['parameters_sha256'])
    print('This is an exact arithmetic computation, not a Lean certificate.')


if __name__ == '__main__':
    main()
