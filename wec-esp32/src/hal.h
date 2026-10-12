#ifndef HAL_H
#define HAL_H
#define PIN_IMU_SDA    21   /* I2C IMU (MPU6050) */
#define PIN_IMU_SCL    22
#define PIN_V_SENSE    34   /* ADC, generator voltage */
#define PIN_LOAD_PWM   25   /* PWM to load, duty follows Bpto */
#define USE_SIM_SERIAL 1    /* 1: sensor frames from Python, 0: real IMU + ADC */

#endif