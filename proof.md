# A six-level candidate upper bound for Turán's pure power-sum constant

**Research draft, 25 September 2026**

## Status and claim

This note gives an analytic construction and an exact-arithmetic certificate for the candidate bound

\[
\boxed{C_{42}\le \frac{68897}{100000}=0.688970.}
\]

The finite inequality used by the construction has passed an outward-rounded rational interval computation. A second implementation, using hypergeometric functions and adaptive quadrature rather than the certificate's series, agrees. The complete asymptotic and algebraic argument is written below. Neither that argument nor the arithmetic program has been formalized in Lean or independently reviewed by a domain expert. No explicit finite threshold is claimed.

The terminology **candidate bound** describes this review status. There is no conjectural mathematical assumption in the argument below.

### Provenance

The starting point was Sebastian Griego's two-block construction in `sebastian-griego/turan-c42-certificate`, version `v1.0.0`, commit `c8ddce14d9a5e898406d5dc6b8d08bb8a39507c7`. In particular, the generating-function construction, the final block used to cancel a coefficient, and the Newton-identity conversion are inherited ideas, not discoveries claimed here. The modification is a six-level middle profile, an explicit finite-step criterion, and the rational parameter certificate below. The asymptotic bridge is included rather than assumed from the earlier computational certificate.

The current table in Tao's optimization-constants repository lists the Harcos upper bound `0.69368` and Griego's proposed bound `0.6906538`, the latter marked with an asterisk. The proposed bound in this note is smaller by exactly `0.0016838`. Searches made during this session did not locate a previously published bound as small as the present one. This is not a guarantee of priority.

## 1. The precise problem

For an integer \(n\ge1\), let

\[
R_n=\min_{\substack{z_1,\ldots,z_n\in\mathbb C\\\max_j|z_j|=1}}
\max_{1\le k\le n}\left|\sum_{j=1}^n z_j^k\right|,
\qquad C_{42}=\limsup_{n\to\infty}R_n.
\]

The condition is **maximum modulus equal to one**, not that every number has modulus one. There is no division of the power sums by \(n\). To obtain an upper bound for \(C_{42}\), it suffices to construct admissible configurations with the stated bound for every sufficiently large integer \(n\), not merely for a subsequence.

This is a bound-improvement problem associated with the classical pure power-sum problem. It does not ask for a new proof of the positive-lower-bound existence question in Erdős Problem 519, which was resolved previously.

## 2. A finite-step sufficient criterion

Choose a real number \(\tau\) with

\[
\frac13\le\tau<\frac12,
\]

thresholds

\[
\tau=t_1<t_2<\cdots<t_M<1-\tau,
\]

and complex numbers \(s,\eta_1,\ldots,\eta_M\). Write

\[
\alpha=1-s,\qquad a=\operatorname{Re}\alpha>0,
\]

and define the jumps

\[
w_1=\eta_1-s,\qquad w_j=\eta_j-\eta_{j-1}\quad(j\ge2).
\]

For positive real arguments, complex powers use the real logarithm. Set

\[
I_\alpha(x)=\int_0^x\frac{v^{\alpha-1}}{1-v}\,dv,
\quad K=I_\alpha(\tau),\quad D=I_a(\tau)>0,
\]

and

\[
L_j=I_\alpha(1-t_j)-I_\alpha(\tau).
\]

For \(t,u>0\), define \(Q_\alpha(t,u)=0\) when \(t+u\ge1\). Otherwise put

\[
Q_\alpha(t,u)=
\int_0^{1-t-u}\frac{v^{\alpha-1}}{1-v}
\log\!\left(\frac{(1-t-v)(1-u-v)}{tu}\right)\,dv.
\tag{2.1}
\]

Equivalently,

\[
Q_\alpha(t,u)=
\int_{\substack{x\ge t,\ y\ge u\\x+y\le1}}
\frac{(1-x-y)^{\alpha-1}}{xy}\,dx\,dy.
\tag{2.2}
\]

To verify the equivalence, fix \(v=1-x-y\), integrate
\(1/[x(1-v-x)]\) from \(t\) to \(1-v-u\), and use partial fractions. The resulting inner integral is exactly the logarithmic factor in (2.1), divided by \(1-v\). All these integrals are absolutely convergent because \(a>0\).

Define

\[
Y=1-\sum_{j=1}^M w_jL_j
+\frac12\sum_{j=1}^M\sum_{k=1}^M w_jw_kQ_\alpha(t_j,t_k)+sK.
\tag{2.3}
\]

**Proposition.** If a positive real number \(C\) satisfies

\[
|s|<C,\qquad |\eta_j|<C\quad(1\le j\le M),\qquad |Y|<CD,
\tag{2.4}
\]

then \(C_{42}\le C\).

Sections 3 through 6 prove this proposition. Sections 7 and 8 certify a specific rational input satisfying its hypotheses.

## 3. Prescribing the first and middle power sums

Fix a sufficiently large integer \(n\), and put \(A=\lfloor\tau n\rfloor\). Partition the indices \(1,\ldots,n\) into

\[
B=\{1,\ldots,A\},\quad
J=\{A+1,\ldots,n-A-1\},\quad
F=\{n-A,\ldots,n\}.
\]

For \(m\in B\), prescribe \(S_m=s\). For \(m\in J\), prescribe

\[
S_m=s+\sum_{j=1}^M w_j\mathbf 1_{m>\lfloor t_j n\rfloor}.
\tag{3.1}
\]

The first threshold equals \(\tau\), so these middle values are among the \(\eta_j\). Distinct thresholds produce nonempty blocks for all sufficiently large \(n\); empty blocks at smaller \(n\) would not affect the argument.

The final \(S_m\), for \(m\in F\), will be selected to cancel one polynomial coefficient. Define

\[
\beta_0=1,\qquad
\beta_l=\frac{(\alpha)_l}{l!}\quad(l\ge1),
\]

so that \((1-z)^{-\alpha}=\sum_{l\ge0}\beta_lz^l\). Every \(\beta_l\) is nonzero: each factor \(\alpha+j\) has positive real part.

On the middle block let \(g_m=(S_m-s)/m\). For a final index let

\[
P_m=\frac{n}{m}\beta_{n-m}\ne0.
\]

Set

\[
Y_n=n\beta_n
-n\sum_{m\in J}g_m\beta_{n-m}
+\frac n2\sum_{\substack{m,k\in J\\m+k\le n}}
 g_mg_k\beta_{n-m-k}
+s\sum_{m\in F}P_m.
\tag{3.2}
\]

We will prove

\[
\frac{|Y_n|}{\sum_{m\in F}|P_m|}\longrightarrow\frac{|Y|}{D}.
\tag{3.3}
\]

Because the limiting value is strictly below \(C\), the left side is at most \(C\) for every sufficiently large \(n\). For such \(n\), make the explicit choice

\[
S_m=\frac{Y_n}{\sum_{r\in F}|P_r|}\frac{\overline{P_m}}{|P_m|}
\qquad(m\in F).
\tag{3.4}
\]

This also makes sense when \(Y_n=0\). It gives both

\[
|S_m|\le C,\qquad \sum_{m\in F}S_mP_m=Y_n.
\tag{3.5}
\]

There is no probabilistic or optimization assumption in this final-block choice.

## 4. Cancelling the nth coefficient exactly

Let \(b_0,\ldots,b_n\) be the first \(n+1\) coefficients of the formal power series

\[
B(z)=(1-z)^{-1}
\exp\!\left(-\sum_{m=1}^n\frac{S_mz^m}{m}\right).
\tag{4.1}
\]

Through degree \(n\), this is

\[
B(z)=(1-z)^{-\alpha}
\exp\!\left(-\sum_{m>A}\frac{(S_m-s)z^m}{m}\right).
\tag{4.2}
\]

Terms of degree greater than \(n\) can be ignored. Crucially,

\[
3(A+1)>n,
\]

including when \(\tau=1/3\). Consequently, three correction indices cannot contribute to the coefficient of \(z^n\). Also, a final index plus any correction index is at least

\[
(n-A)+(A+1)=n+1.
\]

Thus the final corrections enter the nth coefficient only linearly, and the middle corrections enter only linearly and quadratically. Expanding (4.2) gives the exact finite identity

\[
b_n=\beta_n-\sum_{m\in J}g_m\beta_{n-m}
+\frac12\sum_{\substack{m,k\in J\\m+k\le n}}
 g_mg_k\beta_{n-m-k}
-\sum_{m\in F}\frac{(S_m-s)\beta_{n-m}}{m}.
\]

By (3.2),

\[
nb_n=Y_n-\sum_{m\in F}S_mP_m=0.
\tag{4.3}
\]

This is why several middle levels are permitted without an uncontrolled infinite expansion. More levels alter the linear and quadratic coefficients, but do not introduce cubic terms.

## 5. The asymptotic bridge, including the singular endpoint

### 5.1 Coefficient asymptotics

There is a nonzero constant \(c_\alpha\) such that

\[
\beta_l=c_\alpha l^{\alpha-1}(1+O(l^{-1})),\qquad l\to\infty.
\tag{5.1}
\]

In particular, for a constant \(B_\alpha\) and all \(l\ge0\),

\[
|\beta_l|\le B_\alpha(l+1)^{a-1},\qquad
\sum_{l=0}^r|\beta_l|=O((r+1)^a).
\tag{5.2}
\]

For completeness, (5.1) follows directly from the product

\[
\beta_l=\prod_{j=1}^l\left(1+\frac{\alpha-1}{j}\right).
\]

All the factors have positive real part. Put \(d=\alpha-1\). The expansion
\(\log(1+d/j)=d/j+O(j^{-2})\) shows that
\(\sum_{j\le l}\log(1+d/j)-dH_l\) converges with an \(O(l^{-1})\) tail. Since \(H_l-\log l\) has a finite limit and an \(O(l^{-1})\) error, exponentiation proves (5.1), with \(c_\alpha\ne0\). Enlarging a constant for finitely many small indices gives (5.2). This argument does not require a numerical Gamma-function evaluation.

### 5.2 One-dimensional sums

For a fixed \(t\in[\tau,1-\tau)\), ordinary Riemann sums using (5.1) give

\[
c_\alpha^{-1}n^{1-\alpha}
\sum_{\substack{m>\lfloor tn\rfloor\\m\le n-A-1}}
\frac{\beta_{n-m}}m
\longrightarrow I_\alpha(1-t)-I_\alpha(\tau).
\tag{5.3}
\]

These middle sums stay away from the singular endpoint \(n-m=0\). For the final block, substitute \(l=n-m\). Then

\[
c_\alpha^{-1}n^{1-\alpha}
\sum_{m\in F}\frac{\beta_{n-m}}m\longrightarrow K,
\tag{5.4}
\]

and

\[
|c_\alpha|^{-1}n^{1-a}
\sum_{m\in F}\frac{|\beta_{n-m}|}m\longrightarrow D.
\tag{5.5}
\]

To justify the endpoint in these two limits, first restrict to \(l\ge\varepsilon n\), where (5.1) is uniform and ordinary Riemann sums apply. For \(0\le l<\varepsilon n\), (5.2) bounds the absolute normalized sum by

\[
O\!\left(n^{-a}(\varepsilon n+1)^a\right)
=O\!\left((\varepsilon+1/n)^a\right).
\]

The omitted limiting integral is likewise \(O(\varepsilon^a)\). Taking \(n\to\infty\) and then \(\varepsilon\downarrow0\) proves the claims. This also treats a floor-induced change in an endpoint.

### 5.3 Pair sums

For fixed thresholds \(t,u\), consider

\[
T_n(t,u)=\sum_{\substack{m>\lfloor tn\rfloor,\ k>\lfloor un\rfloor\\m+k\le n}}
\frac{\beta_{n-m-k}}{mk}.
\]

There are no terms when \(t+u\ge1\), including the equality case. Otherwise let
\(T=\lfloor tn\rfloor+1\), \(U=\lfloor un\rfloor+1\), and group by \(l=n-m-k\). Writing \(N_l=n-l\), the inner weight is

\[
h_{n,l}=\sum_{m=T}^{N_l-U}\frac{1}{m(N_l-m)}
=\frac1{N_l}\left(
\sum_{m=T}^{N_l-U}\frac1m+
\sum_{k=U}^{N_l-T}\frac1k\right).
\tag{5.6}
\]

If \(l/n\to v\in[0,1-t-u]\), then

\[
nh_{n,l}\longrightarrow
\frac1{1-v}\log\!\left(\frac{(1-t-v)(1-u-v)}{tu}\right).
\tag{5.7}
\]

The harmonic-sum approximations are uniform away from the integrable coefficient endpoint \(l=0\); the lower bounds on \(m/n\) and \(k/n\) are fixed and positive. The limiting weight extends continuously by zero at \(v=1-t-u\), so the empty last ranges and floor errors have the same limit.

For an explicit domination argument, each nonempty inner sum contains at most \(n\) terms, and \(m>tn\), \(k>un\). Hence

\[
0\le h_{n,l}\le\frac1{tu\,n}.
\]

Together with (5.2), the contribution of \(0\le l<\varepsilon n\), after absolute normalization by \(n^{1-a}\), is again
\(O((\varepsilon+1/n)^a)\). Truncating, passing to Riemann sums using (5.1) and (5.7), and restoring the endpoint therefore proves

\[
c_\alpha^{-1}n^{1-\alpha}T_n(t,u)\longrightarrow Q_\alpha(t,u).
\tag{5.8}
\]

In the middle-block quadratic sum, the upper restrictions \(m,k\le n-A-1\) are automatic: the other index is at least \(A+1\). Expanding the jumps in (3.1) thus produces exactly the pair sums in (5.8).

### 5.4 Concluding the ratio limit

There are finitely many thresholds. Combining (5.1), (5.3), (5.4), and (5.8) in (3.2) gives

\[
c_\alpha^{-1}n^{-\alpha}Y_n\longrightarrow
1-\sum_jw_jL_j+\frac12\sum_{j,k}w_jw_kQ_\alpha(t_j,t_k)+sK=Y.
\]

Equation (5.5) gives

\[
|c_\alpha|^{-1}n^{-a}\sum_{m\in F}|P_m|\longrightarrow D>0.
\]

Taking absolute values proves (3.3). The estimates apply as \(n\) runs over all positive integers. No subsequence or unproved uniformity assumption is being used.

## 6. Turning the prescribed moments into admissible complex numbers

Define the monic polynomial

\[
p_n(Z)=Z^{n-1}+b_1Z^{n-2}+\cdots+b_{n-1}.
\]

Let its roots, with multiplicity, be \(y_2,\ldots,y_n\), and put \(y_1=1\). The reversed polynomial is

\[
\prod_{j=1}^n(1-y_jz)=(1-z)\sum_{l=0}^{n-1}b_lz^l.
\]

Since \(b_n=0\), the right-hand side agrees through degree \(n\) with
\((1-z)B(z)\), and therefore with

\[
\exp\!\left(-\sum_{m=1}^n\frac{S_mz^m}{m}\right).
\]

Taking formal logarithms, or equivalently applying Newton's identities, yields

\[
\sum_{j=1}^ny_j^k=S_k\qquad(1\le k\le n).
\]

Let \(\Lambda=\max_j|y_j|\). Because \(y_1=1\), \(\Lambda\ge1\). Define \(z_j=y_j/\Lambda\). Then

\[
\max_j|z_j|=1,
\qquad
\left|\sum_jz_j^k\right|=\Lambda^{-k}|S_k|\le C
\quad(1\le k\le n).
\]

Thus \(R_n\le C\) for every sufficiently large \(n\), proving the proposition.

## 7. Rational parameters for the claimed bound

Take

\[
C=0.688970,\quad \tau=\frac13,\quad
s=0.35649157-0.58957038i,
\]

so

\[
\alpha=0.64350843+0.58957038i.
\]

Every terminating decimal in this section is an **exact rational number**, not an approximate specification.

| j | Threshold t_j | Real part of eta_j | Imaginary part of eta_j |
|---|---:|---:|---:|
| 1 | 1/3 | 0.50302478 | -0.47079246 |
| 2 | 0.37946944 | 0.53723152 | -0.43134878 |
| 3 | 0.43625048 | 0.57359831 | -0.38166012 |
| 4 | 0.49726115 | 0.60923543 | -0.32173228 |
| 5 | 0.55869412 | 0.64178485 | -0.25058267 |
| 6 | 0.61709940 | 0.66874564 | -0.16570680 |

The smallest of the seven exact squared-radius slacks is

\[
\min\left(C^2-|s|^2,\ C^2-|\eta_1|^2,\ldots,C^2-|\eta_6|^2\right)
=0.0000001848110012>0.
\]

The exact interval computation implies the following deliberately coarsened enclosures:

\[
\operatorname{Re}Y\in[0.423384919210350080838839338429,
0.423384919210350080838839338431],
\]

\[
\operatorname{Im}Y\in[-0.447258960858649747806289714465,
-0.447258960858649747806289714463],
\]

\[
D\in[0.893899193572628881864334777218,
0.893899193572628881864334777219].
\]

More importantly, it proves the strict rational inequality

\[
C^2D^2-|Y|^2>\frac1{10^7}>0.
\tag{7.1}
\]

The full machine-readable interval output is `verification.json`. The simpler intervals displayed above are not substituted for that output in the program; they are just readable outer enclosures.

For orientation only, the separate high-precision evaluation gives

\[
\frac{|Y|}{D}\approx0.688969860824420794067832588083,
\]

leaving a radius gap of approximately \(1.391755792\times10^{-7}\). This decimal approximation is not needed by the exact checker. Conditions (2.4) follow from the exact checks, so the proposition gives the asserted bound.

## 8. What the exact interval computation proves

### 8.1 Integral series and rigorous tails

For \(0<x<1\),

\[
I_\alpha(x)=x^\alpha\sum_{m=0}^{N-1}\frac{x^m}{\alpha+m}+E_I,
\qquad
|E_I|\le\frac{x^N}{(a+N)(1-x)}.
\tag{8.1}
\]

Here the harmless factor \(x^a\le1\) has been dropped. The estimate uses
\(|\alpha+m|\ge a+m\).

For a pair kernel with \(L=1-t-u>0\), set \(B=1-t\), \(E=1-u\), and

\[
c_0=\log\frac{BE}{tu},\qquad
c_m=c_0-\sum_{h=1}^m\left(\frac1{hB^h}+\frac1{hE^h}\right).
\]

Expanding the logarithms and \(1/(1-v)\) gives

\[
Q_\alpha(t,u)=L^\alpha\sum_{m=0}^{N-1}\frac{c_mL^m}{\alpha+m}+E_Q.
\]

Since \(L<B,E\), this series converges absolutely. The estimate

\[
|c_m|\le c_0+\frac{B^{-m}}{1-B}+\frac{E^{-m}}{1-E}
\]

gives

\[
|E_Q|\le\frac1{a+N}\left[
\frac{c_0L^N}{1-L}
+\frac{(L/B)^N}{t(1-L/B)}
+\frac{(L/E)^N}{u(1-L/E)}\right].
\tag{8.2}
\]

Again \(L^a\le1\) has been dropped. For numerical stability the checker forms
\(q_m=c_mL^m\) by

\[
q_0=c_0,\qquad q_m=Lq_{m-1}-\frac{(L/B)^m+(L/E)^m}{m},
\]

instead of forming the potentially larger unscaled \(c_m\).

The checker uses \(N=160\) in (8.1) and (8.2).

### 8.2 Elementary functions and rounding

All intervals have integer endpoints over the common denominator \(10^{65}\). Addition is exact at this scale. Products and rational scaling are rounded outward using integer floor and ceiling division. Complex arithmetic is performed using rectangles containing the real and imaginary parts.

For an exact positive rational \(x\), write \(x=2^kq\), \(1\le q<2\). The logarithm is enclosed using

\[
\log q=2\sum_{j=0}^{99}\frac{v^{2j+1}}{2j+1}+E,
\quad v=\frac{q-1}{q+1},\quad
0\le E\le\frac{2v^{201}}{201(1-v^2)},
\]

and the same formula encloses \(\log2\). The signed integer multiple \(k\log2\) is combined with outward rounding.

For a complex exponential, repeatedly halve the input until the sum of upper bounds for the absolute real and imaginary parts is at most \(1/2\). The degree-80 Taylor polynomial then has error at most

\[
\frac{2\rho^{81}}{81!},\qquad \rho\le\frac12.
\]

This follows from the exponential Taylor remainder and \(e^{1/2}<2\). Repeated outward-rounded squaring restores the original argument. The checker therefore evaluates \(x^\alpha=\exp(\alpha\log x)\) without uncertified library logarithms, trigonometric functions, or exponentials.

All omitted complex-series errors are enclosed in rectangles containing the corresponding error disks. To check (7.1), the program uses the sum of the maximal squared absolute endpoint values for the two coordinates of \(Y\), and the positive lower endpoint of \(D\). It compares the resulting rationals exactly. No square root is needed.

### 8.3 Executed checks and their limitations

The following were actually run:

- The standard-library-only exact interval verifier: PASS.
- The same verifier under `python -O`: identical JSON output. Required checks do not depend on Python `assert`.
- A 65-digit second implementation using hypergeometric evaluation for \(I\) and transformed quadrature for \(Q\): all 21 independent symmetric pair entries and the final quantities agree with the rigorous enclosures.
- Exact rational complex-arithmetic checks of \(b_n=0\) and Newton identities for every \(3\le n\le36\): PASS. These use an algebra-only closing choice, not a small-n radius claim.
- Deterministic rational tests of outward rounding and two negative controls: PASS. An invalid radius and an invalid limiting inequality were both rejected.

The second numerical implementation is not an independent human review. Finite tests do not prove the general proposition; Sections 3 through 6 supply that argument.

## 9. Finite-n behavior and scope of the result

Floating-point evaluations of the actual discrete final-block ratio gave:

| n | Computed \(|Y_n|/\sum_{m\in F}|P_m|\) |
|---:|---:|
| 300 | 0.708114989052 |
| 3,000 | 0.693397938525 |
| 30,000 | 0.689985798609 |
| 300,000 | 0.689201875871 |
| 1,000,000 | 0.689076917474 |

These values are consistent with convergence, but **none is a finite-n witness for the claimed radius**. This table is included to avoid conflating the limiting certificate with a uniform bound at moderate dimensions. No rigorous explicit integer \(N\) such that the construction works for all \(n\ge N\) has been calculated.

The strict limiting inequality and the proved convergence imply existence of such an \(N\), which is sufficient for a limsup bound. The result does not determine \(C_{42}\), prove an optimality claim for six-level profiles, or resolve an unrelated unsolved power-sum variant.

## 10. Reproduction and attribution

Run the exact certificate with:

```bash
python3 verify_certificate.py
python3 -O verify_certificate.py parameters.json verification_optimized.json
python3 audit_checks.py
```

Optional numerical checks and search require `mpmath`, `numpy`, and `scipy`:

```bash
python3 independent_check.py
python3 finite_n_check.py
python3 search_multilevel.py
python3 -c 'from search_multilevel import fit; fit(6)'
```

Search is heuristic and is not guaranteed to reproduce the same final digits across numerical-library versions. Verification of the stored rational point does not depend on reproducing the optimizer.

Parameter file SHA-256:

```text
ca248d8996539ee84d417e4a8ed8a95de33202f4bf34245bd28d6ef82d0e7cbf
```

This work was generated in a ChatGPT research session at Kerrigan Plata's request. GPT-6 Astra Pro selected the target, derived and implemented the profile extension, ran the searches and checks, and drafted this note. The prior construction is credited above. No external expert, journal, or proof assistant has endorsed the result.

### Sources consulted

1. Tao et al., *Optimization Constants in Mathematics*, constant 42a, accessed 25 September 2026: `https://teorth.github.io/optimizationproblems/constants/42a.html`.
2. Sebastian Griego, *An improved asymptotic certificate for Turan's pure power sum constant C_42*, version `v1.0.0`: `https://github.com/sebastian-griego/turan-c42-certificate/tree/v1.0.0`. The mathematical note and README were read; its verifier was not used as this certificate's arithmetic implementation.
3. Thomas Bloom, *Erdős Problem 519*: `https://www.erdosproblems.com/519`.
4. The linked classical reference is A. Biró, *An upper estimate in Turán's pure power sum problem*, Indagationes Mathematicae 11 (2000), 499–508. The historical Harcos bound is attributed here via the primary-maintained problem tables; this session did not independently reconstruct its numerical certificate.
