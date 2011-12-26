all: audicorr


audicorr: audicorr.c
	gcc -ggdb -Wall -pedantic --std=c99 -lm -lfftw3 -o audicorr audicorr.c

fftwtest: fftwtest.c
	gcc -ggdb -Wall -pedantic --std=c99 -lm -lfftw3 -o fftwtest fftwtest.c
