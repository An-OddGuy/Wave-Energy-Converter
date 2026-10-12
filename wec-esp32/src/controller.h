#ifndef CONTROLLER_H
#define CONTROLLER_H

typedef struct {
    float t, dc, ms, prev, last_cross;
    int   armed, survival;
    float tp_est;
    float amp_rms;
    float bpto;
    float z1, z2, fms;
} ctrl_t;

#ifdef _WIN32
#define CTRL_API __declspec(dllexport)
#else
#define CTRL_API
#endif
CTRL_API void  ctrl_init(ctrl_t *c);
CTRL_API float ctrl_step(ctrl_t *c, float heave_m, float dt);


#endif