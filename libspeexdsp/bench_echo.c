/* bench_echo.c: end-to-end benchmark of the public echo-canceller API
 * (speex_echo_*), the main speexdsp consumer of kiss_fft.  Toggles the
 * runtime RVV FFT dispatch flag between timed passes so one binary
 * measures the whole-API impact of the RVV butterfly kernels.
 *
 * Needs an optimized build to be meaningful (e.g. meson --buildtype=release);
 * the FFT backend must be kiss for the RVV comparison (-Dfft=kiss on
 * floating-point builds).  Usage: bench_echo [frames_per_pass]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <speex/speex_echo.h>

/* Runtime kill-switch for the RVV FFT kernels; only present in builds with
 * runtime dispatch, and hidden in shared-library builds, hence weak so the
 * bench links everywhere.  When it doesn't resolve, the A/B comparison
 * degrades to a single absolute-timing pass. */
#if defined(__GNUC__)
extern int spx_kf_rvv_enabled __attribute__((weak));
#define HAVE_RVV_TOGGLE (&spx_kf_rvv_enabled != NULL)
#else
static int spx_kf_rvv_enabled_dummy;
#define spx_kf_rvv_enabled spx_kf_rvv_enabled_dummy
#define HAVE_RVV_TOGGLE 0
#endif

static double now_sec(void)
{
   struct timespec ts;
   clock_gettime(CLOCK_MONOTONIC, &ts);
   return ts.tv_sec + 1e-9*ts.tv_nsec;
}

static unsigned lcg_state = 0x12345678;
static int lcg_noise(void)
{
   lcg_state = lcg_state*1103515245u + 12345u;
   return (int)(lcg_state >> 17) % 8192 - 4096;
}

typedef struct {
   const char *label;
   int rate;
   int frame;   /* samples per frame */
   int tail;    /* filter length in samples */
} bench_cfg;

/* Frame sizes chosen so the complex FFT the RVV butterflies run exercises
 * different radix mixes. */
static const bench_cfg cfgs[] = {
   { " 8 kHz, 16 ms frame, 128 ms tail (complex FFT 128: 4,4,4,2)",  8000, 128, 1024 },
   { "16 kHz, 10 ms frame,  80 ms tail (complex FFT 160: 4,4,2,5)", 16000, 160, 1280 },
   { "16 kHz, 16 ms frame, 128 ms tail (complex FFT 256: 4,4,4,4)", 16000, 256, 2048 },
   { "48 kHz, 10 ms frame,  80 ms tail (complex FFT 480: 4,4,2,3,5)", 48000, 480, 3840 },
};

/* One timed pass: feed `frames` frames of synthetic far-end noise plus a
 * delayed, attenuated echo so the adaptive filter does real work. */
static double run_pass(const bench_cfg *c, int frames)
{
   lcg_state = 0x12345678;

   SpeexEchoState *st = speex_echo_state_init(c->frame, c->tail);
   int rate = c->rate;
   speex_echo_ctl(st, SPEEX_ECHO_SET_SAMPLING_RATE, &rate);

   int hist_len = c->frame*4;
   spx_int16_t *play = malloc(sizeof(*play)*c->frame);
   spx_int16_t *rec  = malloc(sizeof(*rec)*c->frame);
   spx_int16_t *out  = malloc(sizeof(*out)*c->frame);
   spx_int16_t *hist = calloc(hist_len, sizeof(*hist));
   int hpos = 0, delay = c->frame + c->frame/2;

   double t0 = 0.0;
   int warmup = frames/10 + 1;
   for (int f = -warmup; f < frames; f++) {
      if (f == 0)
         t0 = now_sec();
      for (int i = 0; i < c->frame; i++) {
         play[i] = lcg_noise();
         hist[hpos] = play[i];
         int d = hpos - delay;
         if (d < 0) d += hist_len;
         rec[i] = hist[d]/4 + lcg_noise()/64;
         hpos = (hpos + 1) % hist_len;
      }
      speex_echo_cancellation(st, rec, play, out);
   }
   double dt = now_sec() - t0;

   free(play); free(rec); free(out); free(hist);
   speex_echo_state_destroy(st);
   return dt/frames;
}

int main(int argc, char **argv)
{
   int frames = argc > 1 ? atoi(argv[1]) : 2000;
   if (frames <= 0) frames = 2000;
   int reps = 3;

   int rvv_avail = 0;
   if (!HAVE_RVV_TOGGLE) {
      printf("no runtime RVV FFT dispatch in this build; single pass only\n");
   } else {
      /* the hwcap probe runs once, inside the first FFT alloc; trigger it
       * now so later writes to the flag stick */
      SpeexEchoState *probe = speex_echo_state_init(128, 1024);
      speex_echo_state_destroy(probe);
      rvv_avail = spx_kf_rvv_enabled;
      if (!rvv_avail)
         printf("note: RVV not usable on this machine; both passes are scalar\n");
   }

   printf("%d frames/pass, best of %d passes per mode\n\n", frames, reps);
   for (size_t i = 0; i < sizeof(cfgs)/sizeof(cfgs[0]); i++) {
      const bench_cfg *c = &cfgs[i];
      double best_on = 1e30, best_off = 1e30;
      for (int r = 0; r < reps; r++) {
         if (HAVE_RVV_TOGGLE) {
            spx_kf_rvv_enabled = 0;
            double t = run_pass(c, frames);
            if (t < best_off) best_off = t;
         }
         if (HAVE_RVV_TOGGLE)
            spx_kf_rvv_enabled = rvv_avail;
         double t = run_pass(c, frames);
         if (t < best_on) best_on = t;
      }
      double frame_ms = 1000.0*c->frame/c->rate;
      printf("%s\n", c->label);
      if (HAVE_RVV_TOGGLE) {
         printf("  scalar FFT: %8.1f us/frame  (%.1fx realtime)\n",
                1e6*best_off, frame_ms/(1000.0*best_off));
         printf("  RVV FFT:    %8.1f us/frame  (%.1fx realtime)\n",
                1e6*best_on, frame_ms/(1000.0*best_on));
         printf("  speedup: %.3fx  (%.1f%% of echo-cancel time saved)\n\n",
                best_off/best_on, 100.0*(best_off-best_on)/best_off);
      } else {
         printf("  %8.1f us/frame  (%.1fx realtime)\n\n",
                1e6*best_on, frame_ms/(1000.0*best_on));
      }
   }
   return 0;
}
