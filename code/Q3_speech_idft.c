/* Question 3
 * 10-point IDFT of speech_dft.txt.
 * Writes speech_idft.txt. Pad zeros from Q2 stay on the last segment.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 10
#define PI 3.14159265358979323846

static void idft10(const double *Xr, const double *Xi, double *xr, double *xi)
{
    int n, k;
    for (n = 0; n < N; n++) {
        double acc_r = 0.0;
        double acc_i = 0.0;
        for (k = 0; k < N; k++) {
            double ang = 2.0 * PI * (double)k * (double)n / (double)N;
            double c = cos(ang);
            double s = sin(ang);
            acc_r += Xr[k] * c - Xi[k] * s;
            acc_i += Xr[k] * s + Xi[k] * c;
        }
        xr[n] = acc_r / (double)N;
        xi[n] = acc_i / (double)N;
    }
}

int main(void)
{
    FILE *in;
    FILE *out;
    int seg, k, n, prev, count, nSeg;
    double Xr[N], Xi[N], mag, xr[N], xi[N];
    double max_abs_imag = 0.0;
    char line[256];

    in = fopen("speech_dft.txt", "r");
    if (!in) {
        printf("Could not open speech_dft.txt\n");
        return 1;
    }
    out = fopen("speech_idft.txt", "w");
    if (!out) {
        printf("Could not create speech_idft.txt\n");
        return 1;
    }
    fprintf(out, "# seg n xr xi\n");

    prev = -1;
    count = 0;
    nSeg = 0;
    while (fgets(line, sizeof line, in)) {
        if (line[0] == '#' || line[0] == '\n')
            continue;
        if (sscanf(line, "%d %d %lf %lf %lf", &seg, &k, &Xr[count], &Xi[count], &mag) != 5)
            continue;
        if (k != count) {
            printf("unexpected bin order at segment %d\n", seg);
            return 1;
        }
        count++;
        if (count == N) {
            idft10(Xr, Xi, xr, xi);
            for (n = 0; n < N; n++) {
                fprintf(out, "%d %d %.8g %.8g\n", seg, n, xr[n], xi[n]);
                if (fabs(xi[n]) > max_abs_imag)
                    max_abs_imag = fabs(xi[n]);
            }
            count = 0;
            nSeg++;
            prev = seg;
        }
    }
    fclose(in);
    fclose(out);

    printf("segments inverted = %d\n", nSeg);
    printf("last segment index = %d\n", prev);
    printf("max |imag(IDFT)| = %.3g\n", max_abs_imag);
    printf("wrote speech_idft.txt (%d rows)\n", nSeg * N);
    return 0;
}
