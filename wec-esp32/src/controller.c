#include "controller.h"
#include <math.h>

#define BPTO_MIN     3.0f
#define BPTO_MAX     40.0f
#define SLEW_RATE    5.0f   /* N s/m per second */
#define RMS_SURVIVE  0.06f   /* m, placeholder */

#define M_TOT    7.37f    /* kg, buoy mass + added mass */
#define B_RAD    4.38f    /* N s/m, radiation damping */
#define K_HYD    443.8f   /* N/m, hydrostatic stiffness */
#define TZ_TO_TP 1.0f

static const float TP_TAB[11] = {0.70f,0.75f,0.80f,0.85f,0.90f,0.95f,
                                 1.00f,1.05f,1.10f,1.15f,1.20f};
static const float BP_TAB[11] = {9,7,7,7,12,16,20,25,25,30,30};

static float lookup(float tp)
{
    if (tp <= TP_TAB[0])  return BP_TAB[0];
    if (tp >= TP_TAB[10]) return BP_TAB[10];
    for (int i = 0; i < 10; i++)
        if (tp < TP_TAB[i + 1]) {
            float f = (tp - TP_TAB[i]) / (TP_TAB[i + 1] - TP_TAB[i]);
            return BP_TAB[i] + f * (BP_TAB[i + 1] - BP_TAB[i]);
        }
    return BP_TAB[10];
}

void ctrl_init(ctrl_t *c)
{
    c->t = c->dc = c->ms = c->prev = c->last_cross = 0.0f;
    c->z1 = c->z2 = c->fms = 0.0f;
    c->armed = c->survival = 0;
    c->tp_est = 0.85f;
    c->amp_rms = 0.0f;
    c->bpto = 7.0f;
}

float ctrl_step(ctrl_t *c, float x, float dt)
{
    c->t += dt;

    c->dc += (dt / 5.0f) * (x - c->dc);
    float y = x - c->dc;
    c->ms += (dt / 20.0f) * (y * y - c->ms);
    c->amp_rms = sqrtf(c->ms);

    /* rebuild wave excitation force from heave, velocity, acceleration */
    float v = (y - c->z1) / dt;
    float a = (y - 2.0f * c->z1 + c->z2) / (dt * dt);
    c->z2 = c->z1;
    c->z1 = y;
    float f = M_TOT * a + (B_RAD + c->bpto) * v + K_HYD * y;
    c->fms += (dt / 20.0f) * (f * f - c->fms);
    float h = 0.1f * sqrtf(c->fms);
    if (h < 1e-3f) h = 1e-3f;

    /* upward zero crossings of the force -> period */
    if (f < -h) c->armed = 1;
    if (c->t > 2.0f && c->armed && c->prev < 0.0f && f >= 0.0f) {
        float tc = c->t - dt * f / (f - c->prev);
        if (c->last_cross > 0.0f) {
            float per = (tc - c->last_cross) * TZ_TO_TP;
            if (per > 0.3f && per < 3.0f)
                c->tp_est += 0.1f * (per - c->tp_est);
        }
        c->last_cross = tc;
        c->armed = 0;
    }
    c->prev = f;

    if (c->amp_rms > RMS_SURVIVE) c->survival = 1;
    else if (c->amp_rms < 0.8f * RMS_SURVIVE) c->survival = 0;

    float target = c->survival ? BPTO_MAX : lookup(c->tp_est);
    if (target < BPTO_MIN) target = BPTO_MIN;
    if (target > BPTO_MAX) target = BPTO_MAX;
    float step = SLEW_RATE * dt;
    float d = target - c->bpto;
    if (d >  step) d =  step;
    if (d < -step) d = -step;
    c->bpto += d;
    return c->bpto;
}

CTRL_API float ctrl_tp(ctrl_t *c) { return c->tp_est; }