# ------------------------------------------------------------------------------
# 10/18/2022: Preliminary sims for CCEDID paper to compare imputed CCE 
# estimator to TWFE.
# direct effect = 2
# mediated effect = 2 (thru X2)


using Random, Distributions, Statistics, LinearAlgebra, Plots
using DataFrames, Chain
rng = MersenneTwister(12);

N = 200; N1 = 100; N0 = N - N1; 
T = 8; T0 = 5;
tau = 2; 
b = [1.0, 1.0];
K = 2;
S = 1000;

f1 = ones(T);
f2 = 1.0:T;
f = [f1 f2];
p = size(f, 2);

# Correlation matrix for errors
C = zeros(T, T); 
for t = 1:T
    for s = 1:T
        C[t, s] = .75^abs(t - s);
    end
end


s = 1;
ccedid = zeros(Float32, T-T0, S);
twfe = zeros(Float32, T-T0, S);
betas = zeros(Float32, S, K);


for s in 1:S
  
  # Generate data --------------------------------------------------------------
    
  id = repeat(1:N, inner = T)
  t = repeat(1:T, outer = N)

  treat = repeat([0.], N * T)
  treat[(t .> T0) .& (id .> N0)] .= 1
  ever_treated = [repeat([0], N0 * T); repeat([1], N1 * T)] |> vec;
  taus = tau * treat;


  # factor loadings (G is for X, g is for y)
  G = zeros(N * p, K);
  g = zeros(N * p, 1);

  for i in 1:N
    draw = randn(rng, 4);
    
    idx = ((i - 1) * p + 1):(i * p);
    G[idx, :] = [1 + draw[1] draw[2]; 
                 draw[3]     1 + draw[4]];

    d = MvNormal(diag(G[idx, :]), I(p));
    g[idx, :] = rand(rng, d, 1)';
    if i > N0
      g[idx, :] = 1 .+ g[idx, :]
    end
  end

  # Generating X
  d = MvNormal(zeros(T), C);
  V1 = rand(rng, d, N) |>
    x -> reduce(vcat, x);
  V2 = rand(rng, d, N) |> 
    x -> reduce(vcat, x);

  # typical element: f * G[1:2, :] + [V1 V2][1:3, 1:2]
  X = kron(I(N), f) * G + [V1 V2];
  
  # Mediated effect
  X[:, 2] = X[:, 2] + taus;


  # Generating y(0)
  d = MvNormal(zeros(T), C); 
  u = rand(rng, d, N) |>
    x -> reduce(vcat, x);
  
  y0 = X * b0 + 4 *  kron(I(N), f) * g + u;

  # Generate y
  y = y0 + tau * treat;

  # data = DataFrame(
  #   id = id, t = t, treat = treat, ever_treated = ever_treated,
  #   y0 = vec(y0), y = vec(y), 
  #   X1 = X[:, 1], X2 = X[:, 2]
  # );

  # data[[1:3; 598:600], :]


  # CCEDID ---------------------------------------------------------------------
  
  # Fhat = average of X and Y for control group
  Fhat = zeros(T, K) 
  for i = i:N0
    idx = (i-1)*T+1:i*T;
    Fhat = Fhat + X[idx, :] / N0;
  end
  Fpre = Fhat[1:T0, :];
  Fpost = Fhat[(T0+1):T, :];

  M_Fpre = I(T0) - Fpre * inv(Fpre' * Fpre) * Fpre';

  # Estimate CCEP for beta_hat
  B = zeros(K, K)
  A = zeros(K, 1)
  for i = 1:N
    idx_low = (i-1) * T + 1
    idx_high = i * T - T0 + 2
    B = B + X[idx_low:idx_high, :]' * M_Fpre * X[idx_low:idx_high, :];
    A = A + X[idx_low:idx_high, :]' * M_Fpre * y[idx_low:idx_high, :];
  end

  # CCEP estimator using pre-treatment X's for all groups
  bhat = inv(B) * A;
  betas[s, :] = bhat';


  # # imputed factor loadings 
  # # (NOTE: The dimension depends on how many factor proxies. For now it's K because I only use Xbar)
  ghat = zeros(N * K, 1);

  # imputed covariates
  Xhat = zeros(N * T, K);

  for i = (N0+1):N
    idx_ghat = ((i-1) * K + 1):(i * K);
    idx_i_pre = ((i-1) * T + 1):((i-1) * T + T0);
    idx_i_post = ((i-1) * T + T0 + 1):(i * T);

    ghat[idx_ghat, :] = 
      inv(Fpre' * Fpre) * Fpre' * 
        (y[idx_i_pre, :] - X[idx_i_pre, :] * bhat);

    Xhat[idx_i_post, :] = 
      Fpost * inv(Fpre' * Fpre) * Fpre' * 
        X[idx_i_pre, :];
  end


  y0hat = zeros(N * T);
  for i = (N0 + 1):N
    idx_low = (i-1) * T + T0 + 1
    idx_high = (i * T)

    y0hat[idx_low:idx_high] = 
      Xhat[idx_low:idx_high, :] * bhat + 
      Fpost * ghat[((i-1) * K + 1):(i * K), :] 

    ccedid[:, s] = ccedid[:, s] + ( 
      y[idx_low:idx_high] - y0hat[idx_low:idx_high]
    );
  end
  ccedid[:, s] = ccedid[:, s] / N1;



  # Imputated TWFE Estimator ---------------------------------------------------
  
  iota_Pre = [ones(T0, 1); zeros(T-T0, 1)]
  iota = ones(T, 1)
  M = I(T) - iota * inv(iota_Pre' * iota_Pre) * iota_Pre';

  y0bar = zeros(T, 1);
  for i = 1:N0
    y0bar = y0bar + y[(i-1)*T+1:i*T,];
  end
  y0bar = y0bar / N0;

  y1tilde = zeros(N*T)
  for i = 1:N
    idx = ((i-1) * T + 1):(i * T);
    y1tilde[idx] = M * (y[idx] - y0bar); 
  end

  for i = (N0 + 1):N
    idx = ((i-1) * T + T0 + 1):(i * T);
    twfe[:,s] = twfe[:,s]  + y1tilde[idx,:];
  end
  twfe[:,s] = twfe[:,s] / N1;

  
end


# median(betas, dims = 1)
# mean(betas, dims = 1)

results = DataFrame(
  median_ccedid = median(ccedid, dims = 2) |> vec, 
  mean_ccedid = mean(ccedid, dims = 2) |> vec,
  se_ccedid = std(ccedid, dims = 2) |> vec,
  median_twfe = median(twfe, dims = 2) |> vec,
  mean_twfe = mean(twfe, dims = 2) |> vec,
  se_twfe = std(twfe, dims = 2) |> vec
)




