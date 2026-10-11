#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "controller.h"

int main(void)
{
    ctrl_t c; ctrl_init(&c);
    const float dt = 0.01f;
    float phase = 0.0f;
    for (int i = 0; i < 18000; i++) {
        float t = i * dt;
        float tp = (t < 60) ? 0.8f : (t < 120) ? 1.1f : 0.7f;
        phase += 6.2831853f * dt / tp;
        float noise = 0.001f * ((float)rand() / RAND_MAX - 0.5f);
        float b = ctrl_step(&c, 0.02f * sinf(phase) + noise, dt);
        if (i % 1000 == 999)
            printf("t=%5.0f  true Tp=%.2f  est Tp=%.3f  rms=%.4f  Bpto=%.1f  surv=%d\n",
                   t + dt, tp, c.tp_est, c.amp_rms, b, c.survival);
    }
    return 0;
}