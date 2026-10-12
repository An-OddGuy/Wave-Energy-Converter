model WEC_Heave_Wide
  constant Real g  = 9.81;
  constant Real pi = Modelica.Constants.pi;

  parameter Real m   = 4.524;   // kg
  parameter Real K   = 443.8;   // N/m
  parameter Real Tp  = 1.0;     // s, peak period (swept by Python)
  parameter Real Hs  = 0.04;    // m
  parameter Real gam = 3.3;     // JONSWAP peak factor (1.0 = Pierson-Moskowitz)

  // ---- Capytaine buoy-only sweep ----
  parameter Integer NT = 21;
  parameter Real wT[NT] = {2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0};
  parameter Real AT[NT] = {3.96,3.90,3.84,3.78,3.68,3.57,3.43,3.29,3.15,3.04,2.94,2.87,2.82,2.80,2.78,2.79,2.80,2.81,2.84,2.86,2.89};
  parameter Real BT[NT] = {1.15,1.54,2.02,2.59,3.24,3.93,4.54,4.99,5.24,5.28,5.13,4.83,4.43,3.95,3.44,2.93,2.44,2.01,1.61,1.26,0.961};
  parameter Real FT[NT] = {409,392,371,347,321,294,267,241,216,192,169,148,129,111,95.1,81.0,68.5,58.0,48.4,40.0,32.7};

  function interp1
    input Real xt[:];
    input Real yt[:];
    input Real x;
    output Real y;
  protected
    Integer n = size(xt, 1);
    Integer i = 1;
  algorithm
    if x <= xt[1] then
      y := yt[1];
    elseif x >= xt[n] then
      y := yt[n];
    else
      while xt[i+1] < x loop
        i := i + 1;
      end while;
      y := yt[i] + (yt[i+1] - yt[i]) * (x - xt[i]) / (xt[i+1] - xt[i]);
    end if;
  end interp1;

  // ---- irregular wave components ----
  parameter Integer N = 40;
  parameter Real wmin = 2.0;
  parameter Real wmax = 16.0;
  parameter Real dw   = (wmax - wmin) / N;
  parameter Real wp   = 2 * pi / Tp;

  parameter Real wi[N] = {wmin + (i - 0.5) * dw
                          + 0.4 * dw * ((i * 0.7548776662 - floor(i * 0.7548776662)) - 0.5) for i in 1:N};
    parameter Real ph[N] = {2 * pi * (sin(i * 12.9898) * 43758.5453 - floor(sin(i * 12.9898) * 43758.5453)) for i in 1:N};

  parameter Real Sraw[N] = {g^2 * wi[i]^(-5) * exp(-1.25 * (wp / wi[i])^4)
                            * gam^exp(-(wi[i] - wp)^2 / (2 * (if wi[i] <= wp then 0.07 else 0.09)^2 * wp^2))
                            for i in 1:N};
  parameter Real sc    = (Hs / 4)^2 / (sum(Sraw) * dw);     // scales spectrum to the requested Hs
  parameter Real ai[N] = {sqrt(2 * sc * Sraw[i] * dw) for i in 1:N};
  parameter Real Fi[N] = {interp1(wT, FT, wi[i]) for i in 1:N};

  // radiation terms at the peak frequency
  parameter Real A = interp1(wT, AT, wp);
  parameter Real B = interp1(wT, BT, wp);

  input Real Bpto(start = 5.6);   // N s/m, PTO damping command
  Real z(start = 0);
  Real v(start = 0);
  Real Fe;
  Real eta;
  Real P;                         // instantaneous power, W
  Real E(start = 0);              // energy, J
equation
  Fe  = sum(Fi[i] * ai[i] * cos(wi[i] * time + ph[i]) for i in 1:N);
  eta = sum(ai[i] * cos(wi[i] * time + ph[i]) for i in 1:N);
  der(z) = v;
  (m + A) * der(v) = Fe - B * v - Bpto * v - K * z;
  P = Bpto * v^2;
  der(E) = P;
end WEC_Heave_Wide;
