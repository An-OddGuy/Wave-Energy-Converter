model WEC_Heave
  parameter Real m    = 4.524;   // kg
  parameter Real A    = 2.85;    // kg, added mass
  parameter Real B    = 4.38;    // N s/m, radiation damping at 8 rad/s
  parameter Real K    = 443.8;   // N/m
  parameter Real Fe0  = 127.85;  // N per m of wave amplitude
  parameter Integer N = 20;
  parameter Real wi[N] = {5.10000, 5.30000, 5.50000, 5.70000, 5.90000, 6.10000, 6.30000, 6.50000, 6.70000, 6.90000, 7.10000, 7.30000, 7.50000, 7.70000, 7.90000, 8.10000, 8.30000, 8.50000, 8.70000, 8.90000};
  parameter Real ai[N] = {0.00044, 0.00067, 0.00092, 0.00119, 0.00146, 0.00172, 0.00196, 0.00221, 0.00247, 0.00280, 0.00323, 0.00377, 0.00436, 0.00481, 0.00492, 0.00474, 0.00436, 0.00391, 0.00347, 0.00310};
  parameter Real Fi[N] = {260.98886, 250.45452, 239.92018, 229.81316, 219.70614, 209.85798, 200.26868, 190.67938, 181.65760, 172.63581, 163.92902, 155.53725, 147.14547, 139.42910, 131.71273, 124.34348, 117.32134, 110.29920, 103.97973, 97.66027};
  parameter Real ph[N] = {3.21587, 5.97194, 0.90578, 5.96054, 1.95929, 2.65984, 5.20061, 2.57107, 3.45320, 0.17316, 4.73446, 3.38125, 2.07177, 4.95384, 1.90503, 2.84941, 0.84221, 2.53283, 1.27835, 1.64816};
  input Real Bpto(start = 5.6);  // N s/m, PTO damping command    // N s/m, PTO damping (try 2.8, 5.6, 11.2)
  Real z(start = 0);
  Real v(start = 0);
  Real Fe;
  Real eta;                     
  Real P;                        // instantaneous power, W
  Real E(start = 0);             // energy, J
equation
  Fe  = sum(Fi[i] * ai[i] * cos(wi[i] * time + ph[i]) for i in 1:N);
  eta = sum(ai[i] * cos(wi[i] * time + ph[i]) for i in 1:N);
  der(z) = v;
  (m + A) * der(v) = Fe - B * v - Bpto * v - K * z;
  P = Bpto * v ^ 2;
  der(E) = P;
end WEC_Heave;
