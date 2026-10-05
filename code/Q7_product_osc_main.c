/* Question 7
 * Magnitude of X[k]*H[k] for each 10-sample pair, written to the DAC.
 * Same excerpts as Question 6.
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

static float seg[N];
static float Xr[N], Xi[N];
static float Hr[N], Hi[N];
static float Yr[N], Yi[N];

short scope_trace[FIFO_LEN];
short bin_mag[N];

static void dft10(const float *x, float *rr, float *ii)
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
        rr[k] = acc_r;
        ii[k] = acc_i;
    }
}

static void build_trace(void)
{
    int s, n, nSeg, i;
    float acc[N];
    float peak;
    int y;

    nSeg = NSEG_SAMPLES / N;
    for (n = 0; n < N; n++)
        acc[n] = 0.0f;

    for (s = 0; s < nSeg; s++) {
        for (n = 0; n < N; n++)
            seg[n] = speech_ex[s * N + n];
        dft10(seg, Xr, Xi);

        for (n = 0; n < N; n++)
            seg[n] = impulse_ex[s * N + n];
        dft10(seg, Hr, Hi);

        for (n = 0; n < N; n++) {
            float yr = Xr[n] * Hr[n] - Xi[n] * Hi[n];
            float yi = Xr[n] * Hi[n] + Xi[n] * Hr[n];
            Yr[n] = yr;
            Yi[n] = yi;
            acc[n] += sqrtf(yr * yr + yi * yi);
        }
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

    *led = 0x02;                 /* D1 on for this run */
    build_trace();

    for (;;) {
        if (((*dacFifoStatus) & 0x20) != 0) {
            for (i = 0; i < FIFO_LEN; i++)
                *dac = scope_trace[i];
        }
    }
}
