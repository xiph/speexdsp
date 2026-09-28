/* Copyright (C) 2026 Tristan Matthews */
/*
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   - Neither the name of the Xiph.org Foundation nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* bench_echo.c: end-to-end benchmark of the public echo-canceller +
 * preprocessor pipeline (speex_echo_* followed by speex_preprocess_run,
 * as a voice pipeline would chain them), timing the library as built.
 * To measure the RVV kernels, compare against a build with -Drvv=disabled.
 *
 * Needs an optimized build to be meaningful (e.g. meson --buildtype=release);
 * float builds exercise smallft by default, -Dfft=kiss the kiss kernels.
 * Usage: bench_echo [frames_per_pass]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <speex/speex_echo.h>
#include <speex/speex_preprocess.h>

static double now_sec(void)
{
#ifdef CLOCK_MONOTONIC
   struct timespec ts;
   clock_gettime(CLOCK_MONOTONIC, &ts);
   return ts.tv_sec + 1e-9*ts.tv_nsec;
#else
   return (double)clock()/CLOCKS_PER_SEC;
#endif
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
   SpeexPreprocessState *den = speex_preprocess_state_init(c->frame, c->rate);
   speex_preprocess_ctl(den, SPEEX_PREPROCESS_SET_ECHO_STATE, st);

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
      speex_preprocess_run(den, out);
   }
   double dt = now_sec() - t0;

   free(play); free(rec); free(out); free(hist);
   speex_preprocess_state_destroy(den);
   speex_echo_state_destroy(st);
   return dt/frames;
}

int main(int argc, char **argv)
{
   int frames = argc > 1 ? atoi(argv[1]) : 2000;
   if (frames <= 0) frames = 2000;
   int reps = 3;

   printf("%d frames/pass, best of %d passes\n\n", frames, reps);
   for (size_t i = 0; i < sizeof(cfgs)/sizeof(cfgs[0]); i++) {
      const bench_cfg *c = &cfgs[i];
      double best = 1e30;
      for (int r = 0; r < reps; r++) {
         double t = run_pass(c, frames);
         if (t < best) best = t;
      }
      double frame_ms = 1000.0*c->frame/c->rate;
      printf("%s\n  %8.1f us/frame  (%.1fx realtime)\n\n",
             c->label, 1e6*best, frame_ms/(1000.0*best));
   }
   return 0;
}
