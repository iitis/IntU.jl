# Reviewer Report: IntU.jl — Symbolic Integration over Haar Measures

## Overall Assessment

The paper presents a well-engineered Julia package for symbolic Haar integration over a broad family of compact groups and random matrix ensembles. The scope is impressive — covering U(d), SU(d), O(d), Sp(d), circular/Gaussian/Ginibre ensembles, permutation groups, t-designs, Stiefel manifolds, and diagonal unitaries in a single unified framework. The writing is generally clear and the code examples are correct and reproducible. The paper is suitable for TOMS after addressing the following points.

## Major Points

**1. Missing comparison with existing tools (critical for TOMS)**
The introduction mentions RTNI (Mathematica) and Haarpy (Python) as competitors, but the benchmarks section contains no comparative timing data. A TOMS paper is expected to demonstrate advantage over alternatives quantitatively. Even if RTNI runs on a proprietary engine, Haarpy is open-source and directly comparable. I recommend adding a table comparing IntU.jl vs. Haarpy for at least the unitary integrals |U₁₁|^{2k} for k = 3, 4, 5 to substantiate the performance claims. If a direct comparison is infeasible, this should be explicitly discussed.

**2. Scalability limits deserve more discussion**
Section 5.1 notes that exact integration is feasible for degrees up to 2k ≈ 10-12, but the practical ceiling is not systematically explored. What is the maximum degree achievable before hitting memory limits? How does caching affect memory consumption? Adding a brief discussion of memory scaling (not just time) would strengthen the benchmarks.

**3. SU(d) treatment is incomplete and should be stated more carefully**
Section 2.1 (line 297-298) acknowledges that SU(d) integration is equivalent to U(d) only for "balanced" polynomials and that ε-tensor contractions are not covered. This is an important limitation that deserves more prominence — perhaps a dedicated "Limitations" subsection. A reader might assume full SU(d) support from the abstract, which lists SU(d) alongside U(d) without qualification.

**4. Numerical stability of HCIZ (Section 3.2)**
The text mentions a "robust perturbation mechanism to handle numerical degeneracies" (line 409-410) but provides no details. What perturbation is applied? How does it scale with the eigenvalue magnitude? For a software paper, users need to understand the error characteristics. A brief description of the perturbation strategy and its error bounds (even empirical ones) would be valuable.

## Minor Points

**5. Benchmark methodology**
- Table 1 reports "median of 30 runs" for most entries but notes "except for long-running symbolic cases." Which cases used fewer samples, and how many? This should be explicit.
- Table 2 mixes "Mean Time" (matrix integration) with "Median Time" (other tables). Using a consistent statistic across all tables would be cleaner.
- The GinUE entries show 0.01 ms and 0.00 ms, suggesting measurement at the resolution floor. Consider reporting in microseconds for sub-millisecond entries for clarity.

**6. Notation inconsistency for symplectic matrices**
In Table 1, the symplectic group entries use |S₁₁|^k, but S is also used for the COE/CSE symmetric matrices. In the code examples (line 641), `Sp[1,1]` is used. Consider using Sp₁₁ in the table to avoid ambiguity with the circular ensemble notation.

**7. Asymptotic expansion example is underdeveloped**
Section 4.5 shows a single asymptotic expansion. A more compelling example would demonstrate the utility for a problem where the exact answer is combinatorially complex but the leading-order behavior reveals physical insight — e.g., the 1/d expansion of Page's formula for average entanglement entropy (already referenced in the introduction).

**8. The "quantum information utilities" subsection (3.6) is thin**
It lists only `partial_trace`. The plural "utilities" and the text about "a suite of helper functions" (line 491) overpromises relative to the current content. Either add more utilities or adjust the language to accurately reflect the single function.

**9. Error handling and edge cases**
The paper does not discuss what happens when users provide invalid inputs (e.g., odd d for Sp(d), k > d for Stiefel manifolds, non-balanced polynomials for SU(d)). A brief paragraph on input validation and error messages would be appropriate for a TOMS paper, where reproducibility and robustness are key review criteria.

**10. Missing complexity analysis for orthogonal/symplectic Weingarten**
Section 5.1 gives the complexity for the unitary case (|S_k| = k!) but only states that pair partitions grow as (2k-1)!! without connecting this to the actual algorithmic cost. The Bareiss solver and exact rational summation are mentioned but their complexity is not discussed.

## Presentation

**11.** The abstract is long (17 lines). Consider tightening — the enumeration of all supported ensembles could be shortened to "a wide range of compact groups, circular and Gaussian ensembles, and discrete groups."

**12.** Line 853: backticks around `check_library` should use `\texttt{}` for consistency with the rest of the paper.

**13.** The paper would benefit from a software architecture diagram (even a simple one) showing the flow from user input through normalization → expansion → dispatch → Weingarten engine → result assembly.

**14.** Section 4 (Usage examples) is extensive (17 subsections) but many examples show only input/output without interpreting the physical or mathematical significance of the results. Examples O1, O2, S1, Sp1, res_a, res_b, res_c (lines 593-644) list code but don't show output or discuss what the results mean. Either show the outputs and discuss them or consolidate into fewer, more thoroughly explained examples.

**15.** The conclusion mentions "exceptional Lie groups" as future work (line 1080). This is a very ambitious direction — any specific groups in mind (G₂, F₄)? A brief mention would ground this claim.

## Bibliography

The 21 references are appropriate and well-chosen. I note the Symbolics.jl reference is listed as an arXiv preprint (2021) — if it has since been formally published, the citation should be updated.

## Summary

The package represents a significant contribution to the computational random matrix theory toolbox. The unified API across many groups, the symbolic dimension handling, and the trace logic abstraction are genuine innovations. Addressing the comparison with existing tools (point 1) and the SU(d) limitation clarity (point 3) would significantly strengthen the paper. The remaining points are straightforward to address.

**Recommendation**: Minor revision.
