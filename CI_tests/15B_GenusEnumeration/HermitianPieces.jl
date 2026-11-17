# The pairs (L_2, zeta) of the hermitian genera over Z[zeta_d], d >= 3, of
# rank b over Z[zeta_d], positive definite, of scales between the inverse of
# the absolute different and E, whose trace lattice is even with det | N.
# Output: a list of Gram matrices and isometries (acting on the columns:
# g^T G g = G), in turn, for GENUS_GlueIsometry pieces and factors.
#
# For a piece P of a lattice L of determinant det L, glued to the rest along
# a glue killed by R = Res(Phi_d, rest): det P | R^min(rank P, rank rest)
# det L (that is N) and the exponent of D(P) divides R det L (that is E).
#
# Usage: julia HermitianPieces.jl d b N OutFile [E]   (E = N by default)
#        julia HermitianPieces.jl --jobs FileJobs     (lines "d b N E OutFile")
using Hecke

function hermitian_pieces(d::Int, b::Int, N::Int, Escale::Int, out::String)
  E, z = cyclotomic_field_as_cm_extension(d; cached = false)
  K0 = base_field(E)
  OE = maximal_order(E)
  sig = Dict(v => 0 for v in real_places(K0))
  # the trace form is integral when the scale is in the inverse of the
  # absolute different D_{E/K0} D_{K0/Q}
  DKE = sum(ideal(OE, OE(E(K0(a)))) for a in gens(different(maximal_order(K0))))
  glist = hermitian_genera(E, b, sig; min_scale = inv(different(OE) * DKE),
                           max_scale = Escale * OE)
  pieces = []
  n_failed = 0
  for g in glist
    # Hecke can fail on a representative (an overflow of Int64): the file
    # then misses that genus, which the log says
    L = try
      representative(g)
    catch err
      n_failed += 1
      println(out, ": no representative: ", first(sprint(showerror, err), 80))
      continue
    end
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
  println(out, ": ", length(pieces), " pieces of ", length(glist), " genera",
          n_failed > 0 ? ", INCOMPLETE: $(n_failed) without a representative" : "")
  flush(stdout)
end

if ARGS[1] == "--jobs"
  for line in eachline(ARGS[2])
    w = split(line)
    if isempty(w) || startswith(w[1], "#")
      continue
    end
    try
      hermitian_pieces(parse(Int, w[1]), parse(Int, w[2]), parse(Int, w[3]),
                       parse(Int, w[4]), String(w[5]))
    catch err
      println(w[5], ": FAILED: ", first(sprint(showerror, err), 80))
    end
  end
else
  N = parse(Int, ARGS[3])
  hermitian_pieces(parse(Int, ARGS[1]), parse(Int, ARGS[2]), N,
                   length(ARGS) >= 5 ? parse(Int, ARGS[5]) : N, ARGS[4])
end
