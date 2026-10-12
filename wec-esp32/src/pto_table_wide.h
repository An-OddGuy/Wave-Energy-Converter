// Generated from wec-twin/results/pto_table_wide.csv (Hs = 0.04 m, Tp 0.5 .. 2 s)
// column used: best_Bpto. Put in include/ (or src/).
#pragma once

#define PTO_N 16

static const float PTO_TP[PTO_N]   = {0.50f, 0.60f, 0.70f, 0.80f, 0.90f, 1.00f, 1.10f, 1.20f, 1.30f, 1.40f, 1.50f, 1.60f, 1.70f, 1.80f, 1.90f, 2.00f};
static const float PTO_BPTO[PTO_N] = {31.1f, 16.7f, 11.0f, 8.3f, 11.9f, 19.8f, 27.8f, 35.8f, 44.9f, 51.7f, 59.5f, 68.5f, 74.6f, 83.5f, 88.3f, 96.1f};

// linear interpolation, clamped to the first/last entry outside the table
static inline float pto_lookup(float tp) {
  if (tp <= PTO_TP[0])         return PTO_BPTO[0];
  if (tp >= PTO_TP[PTO_N - 1]) return PTO_BPTO[PTO_N - 1];
  for (int i = 0; i < PTO_N - 1; i++) {
    if (tp < PTO_TP[i + 1]) {
      float f = (tp - PTO_TP[i]) / (PTO_TP[i + 1] - PTO_TP[i]);
      return PTO_BPTO[i] + f * (PTO_BPTO[i + 1] - PTO_BPTO[i]);
    }
  }
  return PTO_BPTO[PTO_N - 1];
}