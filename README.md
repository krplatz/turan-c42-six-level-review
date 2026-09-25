# Proposed six-level bound for Turan's pure power-sum constant

**Proposed bound: C_42 <= 0.688970. Seeking independent mathematical review.**

This is a research draft, not an independently established record. The construction extends Sebastian Griego's [two-block method](https://github.com/sebastian-griego/turan-c42-certificate/tree/v1.0.0) using six middle levels. The generating-function, final-block cancellation, and Newton-identity machinery are inherited and credited.

## Read and review

- [Mathematical note](proof.md): precise statement, proposition, proof, and parameters.
- [Reviewer guide](REVIEW.md): focused questions and verification scope.
- [Exact rational parameters](parameters.json) and [certificate](verification.json).

The priority question remains open. No explicit finite threshold N or Lean verification is claimed. This is a limsup bound; a successful limiting certificate alone does not verify the analytical reduction.

## Reproduce

Python 3.10 or newer; the exact verifier and audit require no third-party packages:

```sh
python verify_certificate.py parameters.json verification_recomputed.json
python -O verify_certificate.py parameters.json verification_optimized.json
python audit_checks.py
```

Optional alternative numerical check (requires mpmath):

```sh
python independent_check.py
```

The verifier certifies the limiting inequalities with outward-rounded rational intervals and explicit series tails. The alternative numerical implementation is a cross-check, not a rigorous certificate.

## Attribution and review status

Prepared at Kerrigan Plata's request. GPT-6 Astra Pro generated the proposed extension, proof draft, parameter search, and verification code, according to the original work package. A subsequent Codex session inspected the argument and code and reproduced the numerical checks; this was AI-assisted review, not independent human expert review. The exact verifier passed normally and under optimized Python, regenerated certificates matched the supplied certificate, algebra and negative-control checks passed, and the alternative numerical check agreed.

The original mathematical note is preserved. Its references to optional search and finite-n files describe the larger original work package; this focused review repository includes the proof and verification artifacts only.

Current comparison source: [Optimization Constants, 42a](https://teorth.github.io/optimizationproblems/constants/42a.html). No endorsement by its maintainers is implied.
