#include <stdio.h>
#include <complex.h>

#include <fftw3.h>



int main(int argc, char *argv[]) {
    double *in;
    fftw_complex *out;
    fftw_plan p;
    int N = 256;
    int i;

    printf("Alokuji...\n");
    out = fftw_malloc((N/2+1) * sizeof(fftw_complex));
    in = (double*) out;
    printf("Planuji...\n");
    p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);
    printf("Generuji vstupni data...\n");
    in[0]=1;
    in[1]=2;
    for(i=2; i<N; i++) {
        in[i]=0;
    }
    printf("Pocitam...\n");
    fftw_execute(p);
    printf("Spocteno:\n");
    fftw_destroy_plan(p);
    for(i=0; i<(N/2+1); i++) {
        printf("%03d: %+f %+fI\n", i, creal(out[i]), cimag(out[i]));
    }
    printf("Planuji zpetny prevod...\n");
    p = fftw_plan_dft_c2r_1d(N, out, in, FFTW_ESTIMATE);
    printf("Pocitam...\n");
    fftw_execute(p);
    printf("Spocteno:\n");
    fftw_destroy_plan(p);
    for(i=0; i<N; i++) {
        printf("%03d: %+f\n", i, in[i]/N);
    }

    fftw_free(out);
}
