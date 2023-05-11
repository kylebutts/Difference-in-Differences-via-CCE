# ------------------------------------------------------------------------------
# 10/18/2022: Preliminary sims for CCEDID paper to compare imputed CCE 
# estimator to TWFE. Currently there is no treatment effect

using Random, Distributions, Statistics, LinearAlgebra
using Plots, PrettyTables, CSV
using DataFrames, Chain, FixedEffectModels, Printf

include("helpers.module.jl")

rng = MersenneTwister(12);

# DGP:
# 1. No effect
# 2. Direct effect = 1
# 3. Direct effect = 1 + Mediated effect = 1
dgp = 2
parallel_trends = true

sims = DataFrame(;
  dgp=[2, 3, 2, 3], parallel_trends=[true, true, false, false]
)
# sims = DataFrame(;
#   dgp=[2], parallel_trends=[false]
# )

# parameters
N = 200
T = 9
T0 = 6
b = [1.0, 1.0]
K = 2
S = 250


for sim in eachrow(sims)
  
  dgp = sim.dgp
  parallel_trends = sim.parallel_trends

  twfe = zeros(Float32, T - T0, S)
  twfe_w_covariates = zeros(Float32, T - T0, S)
  imputed = zeros(Float32, T - T0, S)
  ccedid = zeros(Float32, T - T0, S)
  betas = zeros(Float32, S, K)
  println("Starting on: DGP = $(dgp), Parallel Trends = $(parallel_trends)")

  for s in 1:S
    (s % 250 == 0) &&  println(s)

    # Generate data ------------------------------------------------------------
    data = generate_data(rng, dgp, parallel_trends)
    # data[[1:3; 598:600], :]
    CSV.write("data.csv", data)

    # CCEDID -------------------------------------------------------------------
    (bhat, est) = est_ccedid(data)
    betas[s, :] = bhat
    ccedid[:, s] = est

    # TWFE Estimator -----------------------------------------------------------
    twfe[:, s] = est_twfe(data)

    # Imputated TWFE Estimator -------------------------------------------------
    imputed[:, s] = est_did2s(data)

    # TWFE w/ covariates Estimator ---------------------------------------------
    twfe_w_covariates[:, s] = est_twfe_w_covariates(data)

    # Synthetic control Estimator ----------------------------------------------
  end

  # Loop through each col
  results = fill("", 3, 1 + 2 * (T - T0))
  results[:, 1] = [
    "TWFE" "TWFE with Covariates" "CCEDID"
  ]

  true_te = 0 + 1 * (dgp == 2) + 2 * (dgp == 3)

  for i in 1:(T - T0)
    results[1, i * 2] = @sprintf "%.2f" mean(twfe[i, :] .- true_te)
    results[1, i * 2 + 1] = @sprintf "%.2f" sum((twfe[i, :] .- true_te) .^ 2) / S

    results[2, i * 2] = @sprintf "%.2f" mean(twfe_w_covariates[i, :] .- true_te)
    results[2, i * 2 + 1] = @sprintf "%.2f" sum((twfe_w_covariates[i, :] .- true_te) .^ 2) /
      S

    results[3, i * 2] = @sprintf "%.2f" mean(ccedid[i, :] .- true_te)
    results[3, i * 2 + 1] = @sprintf "%.2f" sum((ccedid[i, :] .- true_te) .^ 2) / S
  end

  pretty_table(results)
  tex = matrix_to_table(results)

  filename = "tables/simulation-$(dgp)-parallel_trends-$(parallel_trends).tex"

  if !isfile(filename)
    touch(filename)
  end

  open(filename, "w") do file
    write(file, tex)
  end
end
