// WEC controller firmware (ESP32, Arduino framework / PlatformIO / Wokwi)
//
// Serial protocol with wec_realtime.py:
//   host -> "r,<tp>\n"            board -> "ok\n"            reset estimator, start Bpto = pto_lookup(tp)
//   host -> "r\n"                 board -> "ok\n"            same, but start at B_INIT
//   host -> "<t>,<z>,<p_mW>\n"    board -> "<bpto>,<tp>\n"   t = plant seconds, z = heave (m)
//   (p_mW is optional, only used for the OLED)
//
// Control idea: no prediction. Each time a wave finishes (upward zero-crossing of
// heave), its period gives the latest Tp, pto_lookup() gives Bpto, and that Bpto is
// held for the next wave.
//
// LED  (GPIO25, PWM): brightness follows Bpto
// OLED (SSD1306, I2C 0x3C on GPIO21/22): Tp, Bpto, power
//
// platformio.ini:
//   lib_deps = adafruit/Adafruit SSD1306, adafruit/Adafruit GFX Library
// Put results/pto_table_wide.h in include/ (or src/).

#include <Arduino.h>
#include <Wire.h>
#include <Adafruit_GFX.h>
#include <Adafruit_SSD1306.h>
#include "hal.h"
#include "pto_table_wide.h" // PTO_TP[], PTO_BPTO[], pto_lookup(tp)

// ---------------------------------------------------------------- tuning
#define B_INIT      7.0f      // damping at reset only if the host sends no Tp ("r" without ",<tp>")
#define B_MIN       3.0f      // LED range
#define B_MAX       100.0f    // raised from 40: the wide table goes up to ~96 N s/m
#define HYS         0.0003f   // m, zero-crossing hysteresis, raise if the period jumps around
#define TZ_TO_TP    1.03f     // heave zero-crossing period -> peak period (tz_check.py on the Wide model: ~1.03 for Tp >= 0.9 s)
#define MEAN_TAU    10.0f     // s, slow mean removal on heave
#define TP_ALPHA    0.4f      // 1.0 = previous wave only; 0.3..0.6 = light smoothing over recent waves
#define T_WAVE_MIN  0.30f     // s, reject spurious crossings outside this range
#define T_WAVE_MAX  3.00f

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
// Returns true when a wave has just finished and tp_est was refreshed.
static float prev_t, prev_zc, zmean, last_cross_t, tp_est;
static bool  have_prev, have_cross, armed;

static void est_reset() {
  have_prev = have_cross = armed = false;
  prev_t = prev_zc = zmean = last_cross_t = 0.0f;
  tp_est = 0.0f;
}

static bool est_update(float t, float z) {
  bool new_wave = false;

  if (!have_prev) {
    zmean = z;
    prev_t = t;
    prev_zc = 0.0f;
    have_prev = true;
    return false;
  }
  float dt = t - prev_t;
  if (dt <= 0.0f) return false;                    // out-of-order or duplicate sample
  zmean += (dt / (MEAN_TAU + dt)) * (z - zmean);
  float zc = z - zmean;

  if (zc < -HYS) armed = true;
  if (armed && prev_zc < 0.0f && zc >= 0.0f) {
    float tc = prev_t + (0.0f - prev_zc) / (zc - prev_zc) * dt;
    if (have_cross) {
      float T = tc - last_cross_t;                 // the wave that just finished
      if (T >= T_WAVE_MIN && T <= T_WAVE_MAX) {
        float tp_new = T * TZ_TO_TP;
        tp_est = (tp_est > 0.0f) ? tp_est + TP_ALPHA * (tp_new - tp_est) : tp_new;
        new_wave = true;
      }
    }
    last_cross_t = tc;
    have_cross = true;
    armed = false;
  }
  prev_t = t;
  prev_zc = zc;
  return new_wave;
}

// ---------------------------------------------------------------- controller
// Bpto comes straight from the generated table (pto_table_wide.h, Hs = 0.04 m).
// It is held constant and only re-evaluated when a wave finishes.
static float controller(float tp) {
  if (tp <= 0.0f) return B_INIT;
  return pto_lookup(tp);
}

// ---------------------------------------------------------------- serial
static float bpto = B_INIT;
static char  line[96];
static int   n = 0;

static void handle_line(const char* s) {
  if (s[0] == 'r' && (s[1] == '\0' || s[1] == '\r' || s[1] == ',')) {
    est_reset();
    float tp0 = 0.0f;
    bool have_tp0 = (s[1] == ',' && sscanf(s + 2, "%f", &tp0) == 1 && tp0 > 0.0f);
    bpto = have_tp0 ? pto_lookup(tp0) : B_INIT;      // start the sea at its table value
    hal_set_load(bpto);
    hal_show(0.0f, bpto, 0.0f, true);
    Serial.println("ok");
    return;
  }
  float t, z, p = 0.0f;
  int got = sscanf(s, "%f,%f,%f", &t, &z, &p);
  if (got < 2) return;                              // ignore garbage

  if (est_update(t, z)) {                           // a wave just ended -> plan the next one
    bpto = controller(tp_est);
  }

  Serial.printf("%.2f,%.3f\n", bpto, tp_est);       // reply first, so the OLED never delays it
  hal_set_load(bpto);
  hal_show(tp_est, bpto, p);
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