/* Question 8
 * IDFT of X[k]*H[k]. This is the time waveform, not the 10 bars from Q7.
 */


#include <math.h>
#include "Q6_signals.h"

#define N            10
#define FIFO_LEN     256
#define PI           3.14159265358979323846f
#define DAC_PEAK     16000

volatile int *dac = (volatile int *)0x400002;
volatile unsigned int *dacFifoStatus = (volatile unsigned int *)0x400008;
volatile int *led = (volatile int *)0x600000;

static float seg[N];
static float Xr[N], Xi[N];
static float Hr[N], Hi[N];
static float Yr[N], Yi[N];
static float y_time[NSEG_SAMPLES];

short scope_trace[FIFO_LEN];

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

static void idft10(float *y)
{
    int n, k;

    for (n = 0; n < N; n++) {
        float acc_r = 0.0f;

        for (k = 0; k < N; k++) {
            float ang = 2.0f * PI * (float)k * (float)n / (float)N;
            acc_r += Yr[k] * cosf(ang) - Yi[k] * sinf(ang);
        }
        y[n] = acc_r / (float)N;
    }
}

static void build_trace(void)
{
    int s, n, nSeg, i;
    float peak;
    int y;

    nSeg = NSEG_SAMPLES / N;

    for (s = 0; s < nSeg; s++) {
        for (n = 0; n < N; n++)
            seg[n] = speech_ex[s * N + n];
        dft10(seg, Xr, Xi);

        for (n = 0; n < N; n++)
            seg[n] = impulse_ex[s * N + n];
        dft10(seg, Hr, Hi);

        for (n = 0; n < N; n++) {
            Yr[n] = Xr[n] * Hr[n] - Xi[n] * Hi[n];
            Yi[n] = Xr[n] * Hi[n] + Xi[n] * Hr[n];
        }

        idft10(&y_time[s * N]);
    }

    peak = 1.0e-6f;
    for (i = 0; i < NSEG_SAMPLES; i++) {
        float a = y_time[i];
        if (a < 0.0f)
            a = -a;
        if (a > peak)
            peak = a;
    }

    for (i = 0; i < NSEG_SAMPLES; i++) {
        y = (int)(y_time[i] / peak * (float)DAC_PEAK);
        if (y > 32767)
            y = 32767;
        if (y < -32768)
            y = -32768;
        scope_trace[i] = (short)y;
    }
    for (; i < FIFO_LEN; i++)
        scope_trace[i] = 0;
}

void main(void)
{
    int i;

    *led = 0x04;                 /* D2 on for this run */
    build_trace();

    for (;;) {
        if (((*dacFifoStatus) & 0x20) != 0) {
            for (i = 0; i < FIFO_LEN; i++)
                *dac = scope_trace[i];
        }
    }
}
