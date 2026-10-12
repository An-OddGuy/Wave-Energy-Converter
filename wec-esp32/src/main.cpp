// WEC controller firmware (ESP32, Arduino framework / PlatformIO / Wokwi)
//
// Serial protocol with wec_realtime.py:
//   host -> "r\n"                 board -> "ok\n"            reset estimator
//   host -> "<t>,<z>,<p_mW>\n"    board -> "<bpto>,<tp>\n"   t = plant seconds, z = heave (m)
//   (p_mW is optional, only used for the OLED)
//
// LED  (GPIO25, PWM): brightness follows Bpto
// OLED (SSD1306, I2C 0x3C on GPIO21/22): Tp, Bpto, power
//
// platformio.ini:
//   lib_deps = adafruit/Adafruit SSD1306, adafruit/Adafruit GFX Library

#include <Arduino.h>
#include <Wire.h>
#include <Adafruit_GFX.h>
#include <Adafruit_SSD1306.h>
#include "hal.h"

// ---------------------------------------------------------------- tuning
#define B_INIT      7.0f      // damping at start / until a period is known
#define B_MIN       3.0f      // LED range (same as the CoppeliaSim meter)
#define B_MAX       40.0f
#define HYS         0.0003f   // m, zero-crossing hysteresis, raise if the period jumps around
#define TZ_TO_TP    1.03f     // heave zero-crossing period -> peak period, fitted on the WEC_Heave model
#define TP_LO      0.70f     // estimate is held inside the range the lookup table covers
#define TP_HI       1.20f
#define MEAN_TAU    10.0f     // s, slow mean removal on heave

// ---------------------------------------------------------------- HAL: LED + OLED
static Adafruit_SSD1306 oled(128, 64, &Wire, -1);
static bool oled_ok = false;
static uint32_t last_draw_ms = 0;

static void hal_init() {
  Wire.begin(PIN_IMU_SDA, PIN_IMU_SCL);
  Wire.setClock(400000);
#if ESP_ARDUINO_VERSION_MAJOR >= 3
  ledcAttach(PIN_LOAD_PWM, 5000, 8);
#else
  ledcSetup(0, 5000, 8);
  ledcAttachPin(PIN_LOAD_PWM, 0);
#endif
  oled_ok = oled.begin(SSD1306_SWITCHCAPVCC, 0x3C);
  if (oled_ok) {
    oled.setTextColor(SSD1306_WHITE);
    oled.clearDisplay();
    oled.setTextSize(1);
    oled.setCursor(0, 0);
    oled.println("WEC controller");
    oled.println("waiting...");
    oled.display();
  }
}

static void hal_set_load(float b) {
  float f = constrain((b - B_MIN) / (B_MAX - B_MIN), 0.0f, 1.0f);
  uint32_t duty = (uint32_t)(f * 255.0f);
#if ESP_ARDUINO_VERSION_MAJOR >= 3
  ledcWrite(PIN_LOAD_PWM, duty);
#else
  ledcWrite(0, duty);
#endif
}

static void hal_show(float tp, float b, float p_mw, bool force = false) {
  if (!oled_ok) return;
  uint32_t now = millis();
  if (!force && now - last_draw_ms < 200) return;   // max 5 redraws/s
  last_draw_ms = now;
  oled.clearDisplay();
  oled.setTextSize(1);
  oled.setCursor(0, 0);   oled.printf("Tp:   %.2f s", tp);
  oled.setCursor(0, 12);  oled.printf("Bpto: %.1f", b);
  oled.setCursor(0, 24);  oled.printf("P:    %.2f mW", p_mw);
  oled.display();
}

// ---------------------------------------------------------------- estimator
// Upward zero-crossings of the (mean-removed) heave, interpolated between
// time-stamped samples, so skipped samples do not break the period.
static float prev_t, prev_zc, zmean, last_cross_t, tp_est;
static bool  have_prev, have_cross, armed;
static float sum_T;                                // sum of the wave periods seen so far
static int   n_T;                                  // number of waves averaged

static void est_reset() {
  have_prev = have_cross = armed = false;
  prev_t = prev_zc = zmean = last_cross_t = 0.0f;
  tp_est = 0.0f;
  sum_T = 0.0f;
  n_T = 0;
}

static void est_update(float t, float z) {
  if (!have_prev) {
    zmean = z;
    prev_t = t;
    prev_zc = 0.0f;
    have_prev = true;
    return;
  }
  float dt = t - prev_t;
  if (dt <= 0.0f) return;                          // out-of-order or duplicate sample
  zmean += (dt / (MEAN_TAU + dt)) * (z - zmean);
  float zc = z - zmean;

  if (zc < -HYS) armed = true;
  if (armed && prev_zc < 0.0f && zc >= 0.0f) {
    float tc = prev_t + (0.0f - prev_zc) / (zc - prev_zc) * dt;
      if (have_cross) {                              // plain mean over every wave since the reset
        sum_T += tc - last_cross_t;
        n_T++;
        tp_est = (sum_T / (float)n_T) * TZ_TO_TP;
      }
      last_cross_t = tc;
      have_cross = true;
      armed = false;
  }
  prev_t = t;
  prev_zc = zc;
}

// ---------------------------------------------------------------- controller
// Best Bpto per Tp, from pto_table.csv (Hs = 0.04 m only). The sweep grid is coarse
// (3,5,7,9,12,16,...), so each value is a parabola fit around the best grid point.
// Outside 0.70-1.20 s the ends are held.
static const float LUT_TP[] = {0.70f, 0.75f, 0.80f, 0.85f, 0.90f, 0.95f, 1.00f, 1.05f, 1.10f, 1.15f, 1.20f};
static const float LUT_B[]  = {8.3f,  7.7f,  7.2f,  7.9f,  11.3f, 15.8f, 20.3f, 24.3f, 27.5f, 31.2f, 33.2f};
static const int   LUT_N    = 11;

static float controller(float tp) {
  if (tp <= 0.0f) return B_INIT;
  if (tp <= LUT_TP[0]) return LUT_B[0];
  for (int i = 1; i < LUT_N; i++) {
    if (tp <= LUT_TP[i]) {
      float f = (tp - LUT_TP[i - 1]) / (LUT_TP[i] - LUT_TP[i - 1]);
      return LUT_B[i - 1] + f * (LUT_B[i] - LUT_B[i - 1]);
    }
  }
  return LUT_B[LUT_N - 1];
}

// ---------------------------------------------------------------- serial
static float bpto = B_INIT;
static char  line[96];
static int   n = 0;

static void handle_line(const char* s) {
  if (s[0] == 'r' && (s[1] == '\0' || s[1] == '\r')) {
    est_reset();
    bpto = B_INIT;
    hal_set_load(bpto);
    hal_show(0.0f, bpto, 0.0f, true);
    Serial.println("ok");
    return;
  }
  float t, z, p = 0.0f;
  int got = sscanf(s, "%f,%f,%f", &t, &z, &p);
  if (got < 2) return;                              // ignore garbage

  est_update(t, z);
  float tp_rep = (tp_est > 0.0f) ? constrain(tp_est, TP_LO, TP_HI) : 0.0f;
  bpto = controller(tp_rep);

  Serial.printf("%.2f,%.3f\n", bpto, tp_rep);       // reply first, so the OLED never delays it
  hal_set_load(bpto);
  hal_show(tp_rep, bpto, p);
}

void setup() {
  Serial.begin(115200);
  hal_init();
  est_reset();
  hal_set_load(B_INIT);
}

void loop() {
  while (Serial.available()) {
    char c = (char)Serial.read();
    if (c == '\n') {
      line[n] = '\0';
      if (n > 0) handle_line(line);
      n = 0;
    } else if (n < (int)sizeof(line) - 1) {
      line[n++] = c;
    }
  }
}