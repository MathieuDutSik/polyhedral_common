This test code enumerates perfect forms in T-spaces for various discriminants.

The test uses the serial enumeration program PERF_SerialEnumeratePerfectCones to compute perfect forms in T-spaces defined by discriminants.

Test cases:
- For dimension n=3: discriminants d ∈ {-3, -4, -7, -8, -11, -15, -19, -20, -23, -24}
- For dimension n=4: discriminants d ∈ {-3, -4}

The GetFundamentalInfo function computes the RealImagSum and RealImagProd values for the T-space configuration based on the discriminant value:
- For d ≡ 0 (mod 4): eSum=0, eProd=-d/4 (with validity checks)
- For d ≡ 1 (mod 4): eSum=1, eProd=(1-d)/4

The test verifies that the serial enumeration program correctly identifies and enumerates perfect forms in these T-spaces.

TestGeneratorRelation.g checks the presentations by generators and relators of GL2(Z) and GL3(Z) computed by PERF_SerialPerfectComputation (QUERIES:FileGroupGeneratorRelation):
- The relators hold for the matrices of the generators.
- The abelian invariants are [2,2] for GL2(Z) and [2] for GL3(Z).
- The quotients by the normal closure of E12^N are GL_n(Z) / Gamma(N), of orders 12, 48, 96, 240 for n=2 and N=2..5, and 336, 11232 for n=3 and N=2, 3.
- The numbers of subgroups of index at most 8 of GL2(Z) are the ones of the amalgam D8 *_{C2 x C2} D12.
