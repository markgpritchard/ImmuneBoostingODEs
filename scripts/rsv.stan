


functions {
  real logistic(real x) {
    return 1 / (1 + exp(-x));
  }
  real beta_1(array[] real p) {
    return p[2];
  }
  real gamma(array[] real p) {
    return p[4];
  }
  real mu(array[] real p) {
    return p[5];
  }
  real psi(array[] real p) {
    return p[6];
  }
  real omega(array[] real p) {
    return p[7];
  }
  real calcbeta(vector u, array[] real p, real beta_0) {
    real x1 = u[6];
    return beta_0 * (1 + beta_1(p) * x1);
  }
  real calclambda(vector u, array[] real p, real beta_0) {
    real I = u[2];
    return calcbeta(u, p, beta_0) * I;
  }
  vector sirns(real t, vector u, array[] real p, real beta_0) {
    // lambda 
    real lambda = calclambda(u, p, beta_0);

    // compartments 
    real S = u[1];
    real I = u[2];
    real R1 = u[3];
    real R2 = u[4];
    real R3 = u[5];
    real x1 = u[6];
    real x2 = u[7];

    // differential equations
    vector[8] dudt;
    dudt[1] = 3 * omega(p) * R3 - lambda * S + mu(p) * (1 - S);  // S 
    dudt[2] = lambda * S - (mu(p) + gamma(p)) * I;  // I
    dudt[3] = gamma(p) * I + lambda * psi(p) * (R2 + R3) - (3 * omega(p) * mu(p)) * R1;  // R1
    dudt[4] = 3 * omega(p) * R1 - (3 * omega(p) + lambda * psi(p) + mu(p)) * R2;  // R2
    dudt[5] = 3 * omega(p) * R2 - (3 * omega(p) + lambda * psi(p) + mu(p)) * R3;  // R3
    dudt[6] = -2 * pi() * x2;  // x1
    dudt[7] = 2 * pi() * x1;  // x1 
    dudt[8] = lambda * S;  // cumulative cases 
    return dudt;
  }
  array[] real calcbetazeros(
    array[] real mobilityproportions,
    array[] real p, 
    int N_mobilitychanges
  ) {
    array[N_mobilitychanges] real beta_zeros;
    for (i in 1:N_mobilitychanges) {
      real adjustedreduction = (1 - mobilityproportions[i]) * p[9];
      beta_zeros[i] = p[8] * adjustedreduction;
    }
    return beta_zeros;
  }
  
}

data {
  int<lower=0> N;
  int<lower=0> N_mobilitychanges;
  array[N] int incidence;
  array[N_mobilitychanges] real mobilitytimeindexes;
  array[N_mobilitychanges] real mobilityproportions;
  real r0;
  real mu;
  real omega;
  real t0;

}

parameters {
  real logitbeta_1;
  real phi;
  real loggamma;
  real logpsi;
  real logitbetaprimemultiplier;
  real logitfinalbetaprime;
  real logitproportiondetected;
}

transformed parameters {
  array[11] real p_0;
  p[1] = r0 * (exp(loggamma) + mu);
  p[2] = logistic(logitbeta_1);
  p[3] = phi;
  p[4] = exp(loggamma);
  p[5] = mu;
  p[6] = exp(logpsi);
  p[7] = omega;
  p[8] = r0 * (exp(loggamma) + mu);
  array[N_mobilitychanges] real beta_zeros = calcbetazeros(
    mobilityproportions, p, N_mobilitychanges
  );
  vector[8] u0;
  u0[1] = S0;
  u0[2] = I0;
  u0[3] = R1;
  u0[4] = R2;
  u0[5] = R3;
  u0[6] = cos(2 * pi() * t0 - phi);
  u0[7] = sin(2 * pi() * t0 - phi);
  u0[8] = 0.0;
  array[T] vector[8] modelcases;
  modelcases[1:mobilitytimeindexes[1]] = ode_rk45(
    sirns, 
    u0, 
    t0, 
    ts[1:mobilitytimeindexes[1]], 
    p, 
    p[1]
  );
  for (i in 2:N_mobilitychanges) {
    modelcases[mobilitytimeindexes[i]] = ode_rk45(
      sirns, 
      modelcases(mobilitytimeindexes[i-1]), 
      ts[mobilitytimeindexes[1-i]], 
      { ts[mobilitytimeindexes[i]] }, 
      p, 
      beta_zeros[i]
    );
  }
  modelcases[mobilitytimeindexes[N_mobilitychanges]:T] = ode_rk45(
    sirns, 
    modelcases(mobilitytimeindexes[N_mobilitychanges]), 
    ts[mobilitytimeindexes[N_mobilitychanges]], 
    ts[mobilitytimeindexes[N_mobilitychanges]:T], 
    p, 
    beta_zeros[N_mobilitychanges]
  );
}

model {

}
