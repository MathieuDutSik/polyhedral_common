# the number of hermitian genera over Q(zeta_d) / Q(zeta_d)^+ of rank b,
# positive definite, scale between O_E and N O_E (hermitian_genera of Hecke)
# Usage: julia HermitianGenusCounts.jl OutFile (the data of hermitian_genus_counts_hecke.txt)
using Hecke
out = open(ARGS[1], "w")
for d in [3, 5, 7, 9, 11, 12, 13, 15, 20, 21, 24]
  E, z = cyclotomic_field_as_cm_extension(d; cached = false)
  K0 = base_field(E)
  OE = maximal_order(E)
  sig = Dict(v => 0 for v in real_places(K0))
  for b in 1:3, N in [1, 2, 3, 4, 5, 6, 7, 9, 12, 25, 27]
    c = try
      length(hermitian_genera(E, b, sig; min_scale = 1 * OE, max_scale = N * OE))
    catch err
      -1
    end
    println(out, d, " ", b, " ", N, " ", c)
    flush(out)
  end
end
close(out)
