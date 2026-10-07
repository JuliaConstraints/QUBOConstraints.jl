# QUBOConstraints

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaConstraints.github.io/QUBOConstraints.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaConstraints.github.io/QUBOConstraints.jl/dev)
[![Build Status](https://github.com/JuliaConstraints/QUBOConstraints.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaConstraints/QUBOConstraints.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/JuliaConstraints/QUBOConstraints.jl/branch/main/graph/badge.svg?token=zTso333WD7)](https://codecov.io/gh/JuliaConstraints/QUBOConstraints.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/Q/QUBOConstraints.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/report.html)

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).

## Research roadmap

The planned post-ICN research program, including XCSP3-core coverage, encoding families,
auxiliary-variable bounds, exact synthesis by bounded search, and a deep PerfChecker pass at every
phase, is recorded in
[`RESEARCH_ROADMAP.md`](RESEARCH_ROADMAP.md). It is a research plan, not an implementation-status
claim.

The frozen ICN combinatorics and their QUBO transposition are tracked in
[`docs/research/ICN_TO_QUBO_TRANSPOSITION.md`](docs/research/ICN_TO_QUBO_TRANSPOSITION.md).
The current execution state starts with the reproducible
[`docs/research/Q0_BASELINE.md`](docs/research/Q0_BASELINE.md).
The component, codebook, exactness, oracle, and auxiliary-search contracts established in Q1
are specified in [`docs/research/Q1_CONTRACTS.md`](docs/research/Q1_CONTRACTS.md).

The first executable Q1 component and exhaustive-oracle API is documented in
[`docs/src/components.md`](docs/src/components.md), with a reproducible exact Boolean example.
Measured qualification and remaining limitations are recorded in
[`docs/research/Q1_RESULTS.md`](docs/research/Q1_RESULTS.md).
Q2 follows globally correct unconstrained QUBO formulations. Historical counterexamples are
retained for scientific correction, not compatibility; see
[`docs/research/Q2_CORRECTNESS.md`](docs/research/Q2_CORRECTNESS.md).
Its qualified finite replication/synthesis results are in
[`docs/research/Q2_RESULTS.md`](docs/research/Q2_RESULTS.md).
Structured codebooks, including unary-cardinality, bounded coefficients, Gray and Fibonacci,
are described in [`docs/research/Q3_ENCODINGS.md`](docs/research/Q3_ENCODINGS.md).
Q4 local relations and finite truth constructions are qualified in
[`docs/research/Q4_RESULTS.md`](docs/research/Q4_RESULTS.md).
Q5 arithmetic/counting, semantic witnesses and certified auxiliary bisection are specified in
[`docs/research/Q5_CORRECTNESS.md`](docs/research/Q5_CORRECTNESS.md).
Its finite qualification and measured limits are recorded in
[`docs/research/Q5_RESULTS.md`](docs/research/Q5_RESULTS.md).
The historical semantic review before Q6 is recorded in
[`docs/research/Q6_SEMANTIC_REVIEW.md`](docs/research/Q6_SEMANTIC_REVIEW.md).
Q6 structural constructions and explicitly finite wrappers are locally qualified in
[`docs/research/Q6_RESULTS.md`](docs/research/Q6_RESULTS.md), with remaining scope and
performance limits. Q7 uses the explicitly approved whole-integer-time cumulative profile;
its bounded scheduling/graph/expression results and non-isolated performance observations
are in [`docs/research/Q7_RESULTS.md`](docs/research/Q7_RESULTS.md).
Q8 campaign/restart results are in [`docs/research/Q8_RESULTS.md`](docs/research/Q8_RESULTS.md).
Q9 adds logical solver/precision qualification and a binary-route overlay closing the
nine remaining Q8 bounds: [`docs/research/Q9_RESULTS.md`](docs/research/Q9_RESULTS.md).
These finite slices do not complete arbitrary-arity XCSP3 coverage; see the
[methodological boundary](docs/research/Q9_CONCEPTUAL_REVIEW.md).

Q10 now makes the retained variants explicit as graphs of unary/binary operations
and ternary selection, with each operation exposed as square QUBO factors:
[atomic coverage and compositional proof](docs/research/Q10_ATOMIC_COVERAGE.md).
This is a constructive reference, not a search for the smallest matrices. Primary
encodings are preserved; internal value wires use one-hot. See the
[ICN error and neighbourhood contract](docs/research/Q10_ICN_INTERFACE.md) for the
distinction between projected error, conditioned energy and bit couplings.

Q11 closes the core constraint-form inventory (Syntax 15–57), with guarded
arithmetic, the distribute distinctness requirement, explicit semantic profiles,
and core parameter carriers. See the [complete form ledger and proof](docs/research/Q11_CORE_COMPLETENESS.md)
and [qualification results](docs/research/Q11_RESULTS.md). This backend contract
is distinct from XML import and the JuMP bridge maintained elsewhere.

Q11 alone did not establish that these constructions belonged to a learning
model. Q12 adds an explicit finite discrete DAG decision model, with operation,
connection, constant and output choices exposed to the optimizer. Its serialized
weight witnesses reconstruct the Q11 matrices through the normal generic decoder:
[representability contract](docs/research/Q12_REPRESENTABILITY.md),
[results and limits](docs/research/Q12_RESULTS.md), and
[API example](docs/src/atomic_components.md#learnable-discrete-compositions-q12).
This is an extension of the search architecture, not evidence that the historical
14-square/3-triangle motif model already had this capacity. Efficient discovery,
minimality and a single learned generator across arbitrary instance sizes remain
separate research goals.
