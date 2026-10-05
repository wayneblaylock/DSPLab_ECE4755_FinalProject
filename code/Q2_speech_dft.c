/* Question 2
 * 10-point DFT of the speech samples.
 * Reads speech_f32.bin, writes speech_dft.txt.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 10
#define PI 3.14159265358979323846

static void dft10(const double *x, double *Xr, double *Xi)
{
    int k, n;
    for (k = 0; k < N; k++) {
        double acc_r = 0.0;
        double acc_i = 0.0;
        for (n = 0; n < N; n++) {
            double ang = -2.0 * PI * (double)k * (double)n / (double)N;
            acc_r += x[n] * cos(ang);
            acc_i += x[n] * sin(ang);
        }
        Xr[k] = acc_r;
        Xi[k] = acc_i;
    }
}

int main(void)
{
    FILE *in;
    FILE *out;
    float *x;
    double seg[N], Xr[N], Xi[N];
    long bytes, L;
    int nSeg, pad, s, n, base;

    in = fopen("speech_f32.bin", "rb");
    if (!in) {
        printf("Could not open speech_f32.bin\n");
        return 1;
    }
    fseek(in, 0, SEEK_END);
    bytes = ftell(in);
    fseek(in, 0, SEEK_SET);
    if (bytes <= 0 || (bytes % (long)sizeof(float)) != 0) {
        printf("Bad speech_f32.bin size %ld\n", bytes);
        return 1;
    }
    L = bytes / (long)sizeof(float);
    x = (float *)malloc((size_t)L * sizeof(float));
    if (!x || fread(x, sizeof(float), (size_t)L, in) != (size_t)L) {
        printf("read failed\n");
        return 1;
    }
    fclose(in);

    nSeg = (int)((L + N - 1) / N);
    pad = nSeg * N - (int)L;

    out = fopen("speech_dft.txt", "w");
    if (!out) {
        printf("Could not create speech_dft.txt\n");
        return 1;
    }
    fprintf(out, "# seg k Xr Xi mag\n");

    for (s = 0; s < nSeg; s++) {
        base = s * N;
        for (n = 0; n < N; n++) {
            if (base + n < (int)L)
                seg[n] = (double)x[base + n];
            else
                seg[n] = 0.0;
        }
        dft10(seg, Xr, Xi);
        for (n = 0; n < N; n++) {
            double mag = sqrt(Xr[n] * Xr[n] + Xi[n] * Xi[n]);
            fprintf(out, "%d %d %.8g %.8g %.8g\n", s, n, Xr[n], Xi[n], mag);
        }
    }
    fclose(out);
    free(x);

    printf("samples = %ld\n", L);
    printf("segments = %d\n", nSeg);
    printf("zero-pad on last segment = %d\n", pad);
    printf("wrote speech_dft.txt (%d rows)\n", nSeg * N);
    return 0;
}
