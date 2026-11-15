# The pairs (L_2, zeta) of the hermitian genera over Z[zeta_d], d >= 3, of
# rank b over Z[zeta_d], positive definite, of scales between the inverse of
# the absolute different and N, whose trace lattice is even with det | N.
# Output: a list of Gram matrices and isometries (acting on the columns:
# g^T G g = G), in turn, for GENUS_GlueIsometry pieces.
# Usage: julia HermitianPieces.jl d b N OutFile
using Hecke
d = parse(Int, ARGS[1]); b = parse(Int, ARGS[2]); N = parse(Int, ARGS[3])
out = ARGS[4]
E, z = cyclotomic_field_as_cm_extension(d; cached = false)
K0 = base_field(E)
OE = maximal_order(E)
sig = Dict(v => 0 for v in real_places(K0))
# the trace form is integral when the scale is in the inverse of the
# absolute different D_{E/K0} D_{K0/Q}
DKE = sum(ideal(OE, OE(E(K0(a)))) for a in gens(different(maximal_order(K0))))
glist = hermitian_genera(E, b, sig; min_scale = inv(different(OE) * DKE),
                        max_scale = N * OE)
pieces = []
for g in glist
  L = representative(g)
  Lz, f = trace_lattice_with_isometry(L)
  G = gram_matrix(Lz)
  if !all(isinteger, G) || !all(i -> iseven(Int(G[i, i])), 1:nrows(G))
    continue
  end
  D = Int(det(G))
  if D == 0 || N % D != 0
    continue
  end
  # f acts on the ambient space; on the basis B of Lz: B f B^-1
  B = basis_matrix(Lz)
  f = B * f * inv(B)
  @assert f * G * transpose(f) == G
  push!(pieces, (G, f))
end
open(out, "w") do io
  println(io, 2 * length(pieces))
  for (G, f) in pieces
    n = nrows(G)
    println(io, n, " ", n)
    for i in 1:n
      println(io, join([string(Int(G[i, j])) for j in 1:n], " "))
    end
    println(io, n, " ", n)
    for i in 1:n
      println(io, join([string(Int(f[j, i])) for j in 1:n], " "))
    end
  end
end
println("pieces: ", length(pieces), " of ", length(glist), " genera")
