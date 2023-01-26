function generate_data(rng, dgp, parallel_trends = false)
  
  f1 = ones(T)
  f2 = 1.0:T
  f = [f1 f2]
  p = size(f, 2)

  # Correlation matrix for errors
  C = zeros(T, T)
  for t = 1:T
    for s = 1:T
      C[t, s] = 0.75^abs(t - s)
    end
  end

  id = repeat(1:N, inner=T)
  t = repeat(1:T, outer=N)

  # factor loadings (G is for X, g is for y)
  G = zeros(N * p, K)
  g = zeros(N * p, 1)

  for i in 1:N
    draw = randn(rng, 4)

    idx = ((i-1)*p+1):(i*p)
    G[idx, :] = [
      1+draw[1] draw[2]
      draw[3] 1+draw[4]
    ]

    d = MvNormal(diag(G[idx, :]), I(p))
    g[idx, :] = rand(rng, d, 1)'
  end


  # Generate treatment, increasing probability in g2
  g2 = g[2:2:400]
  g2_range = maximum(g2) - minimum(g2)
  g2 = (g2 ./ g2_range)
  g2 = g2

  if(parallel_trends == false) 
    # Probability of treatment increasing in g2
    prob = 0.5 .+ 1 .* (g2 ./ g2_range)

    # Unconditional probability equal to 0.5
    prob = prob .* 0.5 ./ mean(prob)
  else
    prob = fill(0.5, N)
  end
  
  ever_treated = rand(rng, N) .>= prob  
  ever_treated = repeat(ever_treated, inner = T)
  treat = (ever_treated .== 1) .& (t .> T0)
  taus = tau * treat

  
  # Generating X
  d = MvNormal(zeros(T), C)
  V1 = rand(rng, d, N) |>
       x -> reduce(vcat, x)
  V2 = rand(rng, d, N) |>
       x -> reduce(vcat, x)

  # typical element: f * G[1:2, :] + [V1 V2][1:3, 1:2]
  X = kron(I(N), f) * G + [V1 V2]

  # DGP 3
  if(dgp == 3)
    X[:, 2] = X[:, 2] .+ 1 .* treat
  end


  # Generating y(0)
  d = MvNormal(zeros(T), C)
  u = rand(rng, d, N) |>
      x -> reduce(vcat, x)

  y0 = X * b + 4 * kron(I(N), f) * g + u

  # Generate y
  y = y0

  if(dgp == 2)
    y = y .+ (1 .* treat)
  end
  if(dgp == 3)
    y = y .+ (1 .* treat)
  end

  data = DataFrame(
    id=id, t=t, treat=treat, ever_treated=ever_treated,
    y0=vec(y0), y=vec(y),
    X1=X[:, 1], X2=X[:, 2]
  )

end

function est_ccedid(data)
  
  sort!(data, [:id, :t])

  X_vars = [:X1, :X2]
  X = [data.X1 data.X2]
  y = data.y
  id = data.id
  t = data.t
  ever_treated = data.ever_treated

  # Fhat = cross-sectional averages of X for control group
  Fhat = data |>
    (x -> x[(x.ever_treated.==0), :]) |>
    (x -> groupby(x, :t)) |>
    (x -> combine(x,
      X_vars .=> (x -> mean(x)) .=> X_vars
    )) |>
    (x -> sort(x, :t)) |>
    (x -> x[:, X_vars]) |>
    (x -> Matrix(x))

  Fpre = Fhat[1:T0, :]
  Fpost = Fhat[(T0+1):T, :]

  M_Fpre = I(T0) - Fpre * inv(Fpre' * Fpre) * Fpre'

  # Estimate CCEP for beta_hat
  B = zeros(K, K)
  A = zeros(K, 1)
  for i = 1:N
    Xi = Matrix(data[
      (data.id .== i) .& (data.t .<= T0), X_vars
    ])
    yi = data[
      (data.id .== i) .& (data.t .<= T0), :y
    ]
    B = B + Xi' * M_Fpre * Xi
    A = A + Xi' * M_Fpre * yi
  end

  # CCEP estimator using pre-treatment X's for all groups
  bhat = inv(B) * A

  # imputed factor loadings and covariates
  # (NOTE: The dimension depends on how many factor proxies. For now it's K because I only use Xbar)
  treated_ids = data |>
    x -> x[x.ever_treated .== 1, :id] |>
    x -> unique(x)
  
  ghat = zeros(Float32, N * K, 1)
  Xhat = zeros(Float32, N * T, K)
  y0hat = zeros(Float32, N * T)
  ccedid = zeros(Float32, T - T0)

  for i in treated_ids
    idx_ghat = ((i-1)*K+1):(i*K)
    idx_i_pre = ((i-1)*T+1):((i-1)*T+T0)
    idx_i_post = ((i-1)*T+T0+1):(i*T)

    # impute factor loadings
    ghat[idx_ghat, :] =
      inv(Fpre' * Fpre) * Fpre' *
      (y[idx_i_pre, :] - X[idx_i_pre, :] * bhat)

    # impute X(0)
    Xhat[idx_i_post, :] =
      Fpost * inv(Fpre' * Fpre) * Fpre' *
      X[idx_i_pre, :]

    # impute y(0)
    y0hat[idx_i_post] =
      Xhat[idx_i_post, :] * bhat +
      Fpost * ghat[((i-1)*K+1):(i*K), :]

    # estimate TE
    ccedid = ccedid + 
      (y[idx_i_post] - y0hat[idx_i_post]) / N1
  end

  return (bhat=bhat', ccedid=ccedid)
end

function est_twfe(data)
  data.rel_year = (data.t .- (T0 + 1)) .* (data.ever_treated .== 1)
  
  est = reg(
    data,
    @formula(y ~ rel_year + fe(id) + fe(t)),
    contrasts=Dict(:rel_year => DummyCoding(base=-1))
  )

  post = occursin.(r"rel\_year: [0-9]+", coefnames(est))

  return coef(est)[post]
end

function est_twfe_w_covariates(data)
  data.rel_year = (data.t .- (T0 + 1)) .* (data.ever_treated .== 1)
  
  est = reg(
    data,
    @formula(y ~ rel_year + fe(t) & X1 + fe(t) & X2 + fe(id) + fe(t)),
    contrasts=Dict(:rel_year => DummyCoding(base=-1)),
    save=:fe
  )

  post = occursin.(r"rel\_year: [0-9]+", coefnames(est))

  return coef(est)[post]
end

function est_did2s(data)
  y = data.y

  # Imputation matrix
  iota_Pre = [ones(T0, 1); zeros(T - T0, 1)]
  iota = ones(T, 1)
  Imp = iota * inv(iota_Pre' * iota_Pre) * iota_Pre'
  M = I(T) - Imp


  # untreated cross-sectional averages of y
  y0bar = data |>
          x -> x[(x.ever_treated.==0), :] |>
               x -> groupby(x, :t) |>
                    x -> combine(x,
    :y => (mean) => :y
  ) |>
                         x -> sort(x, :t) |>
                              x -> x[:, :y]

  # Create y1tilde
  data |>
  x -> groupby(x, :id) |>
       x -> transform!(x,
    :y => (y -> vec(M * (y - y0bar))) => :y1tilde
  )

  # Average y1tilde for treated units
  twfe = data |>
    (x -> x[(x.ever_treated.==1).&(x.t.>T0), :]) |>
    (x -> groupby(x, :t)) |>
    (x -> combine(x,
      :y1tilde => mean => :twfe
    )) |>
    (x -> sort(x, :t)) |>
    (x -> x[:, :twfe])

  # y1tilde = similar(y)
  # for i = 1:N
  #   y1tilde[data.id .== i] = 
  #     M * (y[data.id .== i] - y0bar); 
  # end
  # y1tilde
  # twfe = zeros(Float32, T - T0)
  # for i = (N0 + 1):N
  #   idx = ((i-1) * T + T0 + 1):(i * T);
  #   twfe = twfe + y1tilde[idx,:];
  # end
  # twfe = twfe / N1;

  return twfe
end

function matrix_to_table(mat::Matrix{String}) 
  (r, c) = size(mat)
  tex = ""

  for i in 1:r
    for j in 1:c
      sep = (j == 1) ? "" : " & "
      tex *= "$(sep)$(mat[i, j])"
    end
    tex *= " \\\\ \n"
  end

  return tex
end


