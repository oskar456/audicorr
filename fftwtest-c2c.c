#include <stdio.h>
#include <complex.h>

#include <fftw3.h>



int main(int argc, char *argv[]) {
    fftw_complex *in;
    fftw_plan p;
    int N = 256;
    int i;

    printf("Alokuji...\n");
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    printf("Planuji...\n");
    p = fftw_plan_dft_1d(N, in, in, FFTW_FORWARD, FFTW_MEASURE);
    printf("Generuji vstupni data...\n");
    in[0]=1;
    in[1]=2;
    for(i=2; i<=255; i++) {
        in[i]=0;
    }
    printf("Pocitam...\n");
    fftw_execute(p);
    printf("Spocteno:\n");
    for(i=0; i<=255; i++) {
        printf("%03d: %f%+fI\n", i, creal(in[i]), cimag(in[i]));
    }
    fftw_destroy_plan(p);
    fftw_free(in);
}
