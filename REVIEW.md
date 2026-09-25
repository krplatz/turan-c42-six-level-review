# Focused review request

Claim: limsup R_n <= 68897/100000, with R_n defined in Section 1 of proof.md.

Please prioritize the proposition in Sections 2-6:

1. Does the multilevel jump profile give precisely the stated linear and quadratic kernels, including floors and block endpoints?
2. Do the coefficient asymptotics and endpoint domination justify the limits for all sufficiently large integers n?
3. Does coefficient cancellation followed by Newton identities and normalization yield admissible configurations with the asserted moment bound?
4. Are the interval series, tail estimates, and outward-rounding implementation sufficient to certify the numerical hypotheses?
5. Is this improvement already known, or covered by a stronger existing result?

The certified squared-inequality margin exceeds 1/10000000. The numerical radius gap is approximately 1.391755792e-7. No explicit eventual threshold N is supplied. The original report states that the final-block ratio at n=1,000,000 remains above the claimed limiting bound; those finite computations are not presented as witnesses.

A targeted objection or reference is useful even if a full review is not possible. All current checks were performed in AI-assisted sessions. No human expert endorsement is claimed.
