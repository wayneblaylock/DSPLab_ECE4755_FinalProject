/* Question 6
 * 10-point DFT magnitude of a 256-sample excerpt, written to the DAC.
 * D0 released = speech, D0 pressed = impulse.
 * Excerpts are in Q6_signals.h (speech at 43968, impulse at 55680).
 */


#include <math.h>
#include "Q6_signals.h"

#define N            10
#define FIFO_LEN     256
#define PI           3.14159265358979323846f
#define DAC_PEAK     20000

volatile int *dac = (volatile int *)0x400002;
volatile unsigned int *dacFifoStatus = (volatile unsigned int *)0x400008;
volatile int *led = (volatile int *)0x600000;
volatile int *key = (volatile int *)0x600800;

static float seg[N];
static float Xr[N];
static float Xi[N];
static float mag[N];

short scope_trace[FIFO_LEN];
short bin_mag[N];
int   show_impulse;   /* 0 = speech, 1 = impulse */

static void dft10(const float *x)
{
    int k, n;

    for (k = 0; k < N; k++) {
        float acc_r = 0.0f;
        float acc_i = 0.0f;

        for (n = 0; n < N; n++) {
            float ang = -2.0f * PI * (float)k * (float)n / (float)N;
            acc_r += x[n] * cosf(ang);
            acc_i += x[n] * sinf(ang);
        }
        Xr[k] = acc_r;
        Xi[k] = acc_i;
        mag[k] = sqrtf(acc_r * acc_r + acc_i * acc_i);
    }
}

/* Hold each bin for 25 samples so the scope shows 10 bars. */
static void build_trace(const float *x, int len)
{
    int s, n, nSeg, i;
    float acc[N];
    float peak;
    int y;

    nSeg = len / N;
    for (n = 0; n < N; n++)
        acc[n] = 0.0f;

    for (s = 0; s < nSeg; s++) {
        for (n = 0; n < N; n++)
            seg[n] = x[s * N + n];
        dft10(seg);
        for (n = 0; n < N; n++)
            acc[n] += mag[n];
    }

    peak = 1.0e-6f;
    for (n = 0; n < N; n++) {
        acc[n] /= (float)nSeg;
        if (acc[n] > peak)
            peak = acc[n];
    }

    i = 0;
    for (n = 0; n < N; n++) {
        y = (int)(acc[n] / peak * (float)DAC_PEAK);
        if (y > 32767)
            y = 32767;
        bin_mag[n] = (short)y;
        for (s = 0; s < 25; s++)
            scope_trace[i++] = (short)y;
    }
    while (i < FIFO_LEN)
        scope_trace[i++] = 0;
}

void main(void)
{
    int i;
    int prev = -1;

    for (;;) {
        show_impulse = ((*key) & 0x1) ? 1 : 0;
        *led = show_impulse;            /* D0 on when showing impulse */

        if (show_impulse != prev) {
            if (show_impulse)
                build_trace(impulse_ex, NSEG_SAMPLES);
            else
                build_trace(speech_ex, NSEG_SAMPLES);
            prev = show_impulse;
        }

        /* bit 5 set means the DAC FIFO is empty */
        if (((*dacFifoStatus) & 0x20) != 0) {
            for (i = 0; i < FIFO_LEN; i++)
                *dac = scope_trace[i];
        }
    }
}
