#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <math.h>
#include <complex.h>

#include <fftw3.h>

#define VERSION "9999"
#define program_invocation_short_name "audicorr"
#define PACKAGE "audicorr"

#pragma pack(1)
struct	WAV_HEADER
{
	uint8_t		RIFF[4];	    /* RIFF Header = "RIFF"	*/
	uint32_t	ChunkSize;	    /* RIFF Chunk Size = filesize-8 */
	uint8_t		WAVE[4];	    /* WAVE Header = "WAVE"	*/
	uint8_t		fmt[4];	    	/* FMT header  = "fmt "	*/
	uint32_t	Subchunk1Size;  /* Size of the fmt chunk (at least 16)			    */
	uint16_t	AudioFormat;    /* Audio format 1=PCM,6=mulaw,7=alaw, 257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM */
	uint16_t	NumOfChan;	    /* Number of channels 1=Mono 2=Sterio		    */
	uint32_t	SampleRate;     /* Sampling Frequency in Hz				    */
	uint32_t	bytesPerSec;    /* bytes per second */
	uint16_t	bytesPerSample; /* 2=16-bit mono, 4=16-bit stereo */
	uint16_t	bitsPerSample;  /* Number of bits per sample      */
};
struct WAV_HEADER2
{
	uint8_t		Subchunk2ID[4]; /* = "data" */
	uint32_t	Subchunk2Size;  /* Data length    */
}; 
#pragma pack()

enum loglevel {
	LOG_FATAL = 0, /* Always shown */
	LOG_ERROR,     /* Could be silenced */
	LOG_INFO,      /* Default verbosity */
	LOG_DEBUG
};

enum loglevel conf_verbosity = LOG_INFO;
long nfft;	/* FFT width for one iteration */
char *needle_fname;
double match_treshold = 0.95;
int match_end = 0;


/**
 * Logger function. Show the message if current verbosity is above
 * logged level.
 *
 * @param levem Message log level
 * @param format printf style format string
 * @returns Whatever printf returns
 */
int logger(enum loglevel level, const char *format, ...) {
	va_list ap;
	int r;
	if (conf_verbosity >= level) {
		va_start(ap, format);
		r=vfprintf(stderr,format, ap);
		va_end(ap);
		return r;
	}
	return 0;
}

void fatal(const char *message) {
	logger(LOG_FATAL, "%s\n", message);
	exit(EXIT_FAILURE);
}

void *ec_malloc(size_t size) {
	void *ptr;
	ptr = malloc(size);
	if (ptr == NULL) {
		fatal("Cannot allocate memory!");
	}
	return ptr;
}


void usage(FILE* f) {
	fprintf(f, 
"Audicorr - The audio correlator\n"
"\n"
"Version " VERSION "\n"
"Copyright 2011-12 Ondrej Caletka <ondrej@caletka.cz>\n"
"\n"
"This program is free software; you can redistribute it and/or modify\n"
"it under the terms of the GNU General Public License version 2\n"
"as published by the Free Software Foundation.\n"
"\n"
"Usage:	%s [options] <needle wave file> < <haystack wave file>\n"
"\n"
"Options:\n"
"\t-e --end			Match the end of needle, rather than beginning\n"
"\t-h --help --usage		Show this help\n"
"\t-V --version			Show program version\n"
"\t-t --treshold <num>		Match treshold (default 0.95)\n"
"\t-v --verbose			Increase verbosity\n"
"\t-q --quiet			Report only fatal errors\n",
	program_invocation_short_name);
}



void parseCmdLine(int argc, char *argv[]) {
	static const struct option longopts[] = {
		{ "help",	no_argument,		NULL,	'h' },
		{ "usage",	no_argument,		NULL,	'h' },
		{ "quiet",	no_argument,		NULL,	'q' },
		{ "verbose",	no_argument,		NULL,	'v' },
		{ "version",	no_argument,		NULL,	'V' },
		{ "treshold",	required_argument,	NULL,	't' },
		{ "end",	no_argument,		NULL,	'e' },
		{ 0,		0,			0,	0   }
	};
	static const char shortopts[] = "hqvVt:e";

	int option_index, opt;
	
	while ((opt = getopt_long(argc, argv, shortopts,
					longopts, &option_index)) != -1) {
		switch (opt) {
			case 0:
				break;
			case 'h':
				usage(stdout);
				exit(EXIT_SUCCESS);
				break;
			case 'q':
				conf_verbosity=0;
				break;
			case 'v':
				conf_verbosity++;
				break;
			case 'V':
				puts(PACKAGE " " VERSION);
				exit(EXIT_SUCCESS);
				break;
			case 't':
				match_treshold = atof(optarg);
				break;
			case 'e':
				match_end=1;
				break;

			default:
				usage(stderr);
				exit(EXIT_FAILURE);
		}
	}

	if (optind >= argc) {
		logger(LOG_FATAL, "Error, no needle wave found.\n");
		usage(stderr);
		exit(EXIT_FAILURE);
	}
	needle_fname = strdup((const char *) argv[optind]);
	if (++optind < argc) {
		logger(LOG_FATAL, "Error, garbage on command line.\n");
		usage(stderr);
		exit(EXIT_FAILURE);
	}

}


/*Check wave header correctness
 * Fails on a error
 * In case of success, return number of bytes to seek to reach header2*/
int check_wave_header(struct WAV_HEADER *wav_hdr) {
	int fail=0;
	if (strncmp((const char *)wav_hdr->RIFF, "RIFF", 4) != 0)
		fail=1;
	if (strncmp((const char *)wav_hdr->WAVE, "WAVE", 4) != 0)
		fail=1;
	if (strncmp((const char *)wav_hdr->fmt, "fmt ", 4) != 0)
		fail=1;
	if (wav_hdr->AudioFormat != 1) {
		logger(LOG_FATAL,"Unsupported data format, only PCM supported\n");
		fail=1;
	}
	if (wav_hdr->NumOfChan > 1) {
		logger(LOG_ERROR,"Warning: Multichannel input - considering only first channel\n");
	}
	logger(LOG_DEBUG, "Sample rate: %d\n", wav_hdr->SampleRate);
	logger(LOG_DEBUG, "Bits per sample: %d\n", wav_hdr->bitsPerSample);
	logger(LOG_DEBUG, "Bytes per sample: %d\n", wav_hdr->bytesPerSample);
	if (wav_hdr->bitsPerSample != 8 && wav_hdr->bitsPerSample != 16) {
		logger(LOG_FATAL, "Unsuported bit width, only 8 or 16 bits supported");
		fail=1;
	}
	if (fail != 0) {
		fatal("WAVE header check failed!\n");
	}
	return wav_hdr->Subchunk1Size - 16;
}

/**
 * Check for data header. If other header is present, 
 * returns number of bytes to seek the other header, otherwise returns 0
 */
int check_wave_header2(struct WAV_HEADER2 *wav_hdr2) {
	if (strncmp((const char *)wav_hdr2->Subchunk2ID, "data", 4) == 0)
		return 0;
	return wav_hdr2->Subchunk2Size + 8;
}

long samplesnumber(struct WAV_HEADER *wav_hdr, struct WAV_HEADER2 *wav_hdr2) {
	return wav_hdr2->Subchunk2Size / wav_hdr->bytesPerSample;
}

/**
 * Compute suitable FFT window size - a power of 2
 * and at least 4 times multiply of needle size.
 */
long compute_fft_size(long n_needle) {
	double npow;
	long nfft;
	npow = log2(n_needle);
	nfft = pow(2, 2 + round(npow));
	logger(LOG_DEBUG, "Computed FFT size: 2^(2 + %f) = %ld\n", npow, nfft);
	return nfft;
}

/*
 * Read at most nmax data samples
 * Returns characters actually read.
 */
long read_wav_data(FILE *wavfile, struct WAV_HEADER *wav_hdr, double *signal, long nmax) {
	long n;
	int r;
	char *sample;

	sample = (char *) ec_malloc(wav_hdr->bytesPerSample);
	for(n=0; n<nmax; n++) {
		r = fread(sample, wav_hdr->bytesPerSample, 1, wavfile);
		if (r != 1)
			break;
		if (wav_hdr->bitsPerSample == 8) {
			signal[n] = ( *((uint8_t *)sample) - 128) / 128.0;
		} else if (wav_hdr->bitsPerSample == 16) {
			signal[n] = *((int16_t *)sample) / 32768.0;
		}
	}
	free(sample);
	return n;
}

/*
 * Read at most nmax data samples. If less than nmax data samples are read,
 * remaining space is padded with zeroes.
 * Returns 0 if EOF, nmax else.
 */
long read_padded_wav_data(FILE *wavfile, struct WAV_HEADER *wav_hdr,
						double *signal, long nmax) {
	long read,n;
	read = read_wav_data(wavfile, wav_hdr, signal, nmax);
	if (read == 0) return 0;
	for (n=read; n<nmax; n++)
		signal[n]=0;
	return nmax;
}

int main(int argc, char *argv[])
{
	FILE *fneedle, *fhaystack = stdin;
	struct WAV_HEADER needle_hdr, haystack_hdr;
	struct WAV_HEADER2 needle_hdr2, haystack_hdr2;
	long n_needle, n, max_n, n0;
	fftw_plan plan, backplan;
	fftw_complex *needle_spec, *haystack_spec;
	double *needle_sig, *haystack_sig, *haystack_safe;
	double needle_energy, round_max_energy, round_match_time;
	double max_energy = 0, match_time = -1;
	int nseek;

	parseCmdLine(argc, argv);

	/* Prepare the needle */
	fneedle = fopen(needle_fname, "r");
	if (fneedle == NULL) {
		logger(LOG_FATAL, "Cannot open needle file!\n");
		exit(EXIT_FAILURE);
	}
	if (fread(&needle_hdr, sizeof(needle_hdr), 1, fneedle) != 1) {
		logger(LOG_FATAL, "Header read failure!\n");
		exit(EXIT_FAILURE);
	}
	nseek =  check_wave_header(&needle_hdr);
	fseek(fneedle, nseek, SEEK_CUR);
	do {
		if (fread(&needle_hdr2, sizeof(needle_hdr2), 1, fneedle) != 1) {
			logger(LOG_FATAL, "Header read failure!\n");
			exit(EXIT_FAILURE);
		}
		nseek = check_wave_header2(&needle_hdr2);
		fseek(fneedle, nseek, SEEK_CUR);
	} while (nseek > 0);

	
	n_needle = samplesnumber(&needle_hdr, &needle_hdr2);
	nfft = compute_fft_size(n_needle);

	needle_spec = (fftw_complex*) fftw_malloc((nfft/2+1) * sizeof(fftw_complex));
	if (needle_spec == NULL) {
		fatal("Allocation failed!\n");
	}
	needle_sig = (double *) needle_spec;

	logger(LOG_DEBUG, "Preparing needle FFT...\n");
	plan = fftw_plan_dft_r2c_1d(nfft, needle_sig, needle_spec, FFTW_ESTIMATE);
	memset((void *) needle_spec, 0, (nfft/2+1) * sizeof(fftw_complex));
	n_needle = read_wav_data(fneedle, &needle_hdr, needle_sig, nfft);
	needle_energy = 0;
	for (n=0; n<n_needle; n++) {
		needle_energy += needle_sig[n] * needle_sig[n];
	}
	logger(LOG_DEBUG, "Needle total energy; %f\n", needle_energy);
	logger(LOG_DEBUG, "Computing needle FFT\n");
	fftw_execute(plan);
	fftw_destroy_plan(plan);

	/* Prepare the complex conjugate of needle */
	for(n=0; n<(nfft/2+1); n++) {
		needle_spec[n] = conj(needle_spec[n]);
	}

	logger(LOG_DEBUG,"Needle length: %ld\nFFT length:%ld\n", n_needle, nfft);

	/* Prepare the haystack */
	if (fhaystack == NULL) {
		logger(LOG_FATAL, "Cannot open haystack file: %s\n", strerror(errno));
		exit(EXIT_FAILURE);
	}
	if (fread(&haystack_hdr, sizeof(haystack_hdr), 1, fhaystack) != 1) {
		logger(LOG_FATAL, "Header read failure: %s\n", strerror(errno));
		exit(EXIT_FAILURE);
	}
	nseek =  check_wave_header(&haystack_hdr);
	fseek(fhaystack, nseek, SEEK_CUR);
	do {
		if (fread(&haystack_hdr2, sizeof(haystack_hdr2), 1, fhaystack) != 1) {
			logger(LOG_FATAL, "Header read failure!\n");
			exit(EXIT_FAILURE);
		}
		nseek = check_wave_header2(&haystack_hdr2);
		fseek(fhaystack, nseek, SEEK_CUR);
	} while (nseek > 0);
	haystack_spec = (fftw_complex*) fftw_malloc((nfft/2+1) * sizeof(fftw_complex));
	if (haystack_spec == NULL) {
		fatal("Allocation failed!\n");
	}
	haystack_sig = (double *) haystack_spec;
	logger(LOG_DEBUG, "Preparing haystack FFT\n");
	plan = fftw_plan_dft_r2c_1d(nfft, haystack_sig, haystack_spec, FFTW_ESTIMATE);
	logger(LOG_DEBUG, "Preparing iFFT\n");
	backplan = fftw_plan_dft_c2r_1d(nfft, haystack_spec, haystack_sig, FFTW_ESTIMATE);
	memset((void *) haystack_spec, 0, (nfft/2+1) * sizeof(fftw_complex));

	/* haystack_safe is used to store last n_needle samples of haystack for
	 * overlaps between iterations
	 */
	haystack_safe = (double *) fftw_malloc(n_needle * sizeof(double));
	if (haystack_safe == NULL) {
		logger(LOG_FATAL,"Allocation failed!\n");
		exit(EXIT_FAILURE);
	}
	
	read_padded_wav_data(fhaystack, &haystack_hdr, haystack_sig, nfft);
	memcpy((void *) haystack_safe, (void *) (haystack_sig+nfft-1-n_needle),
			n_needle * sizeof(double));
	n0 = 0;
	while (1) {
		logger(LOG_DEBUG, "Computing haystack FFT\n");
		fftw_execute(plan);

		//xcorr(X,Y) = ifft( fft(X) * conj(fft(Y)) )
		for(n=0; n<(nfft/2+1); n++) {
			haystack_spec[n] *= needle_spec[n];
		}
		logger(LOG_DEBUG, "Computing iFFT\n");
		fftw_execute(backplan);
		max_n = -1;
		round_max_energy = 0;
		for(n=0; n<nfft-n_needle; n++) {
			if (round_max_energy < haystack_sig[n]) {
				round_max_energy = haystack_sig[n];
				max_n = n;
			}
		}
		/* Normalize maximum energy to 1.0 ~ same volume on haystack and needle */
		round_max_energy /= nfft * needle_energy; 
		round_match_time = ((double) max_n + n0 + (match_end?n_needle:0)) /
					needle_hdr.SampleRate;
		if (round_max_energy > max_energy) {
			max_energy = round_max_energy;
			match_time = round_match_time;
		}

		logger(LOG_DEBUG, "In sample from %ld to %ld:\n", n0,
				(n0+nfft-n_needle));
		logger(LOG_INFO, "In time from %f to %f:\n", (double) n0/needle_hdr.SampleRate,
				(double) (n0+nfft-n_needle) / needle_hdr.SampleRate);
		if (round_max_energy > match_treshold) {
			logger(LOG_INFO, "Match found, energy %f, time %f, sample %ld\n", round_max_energy,
					round_match_time, max_n+n0);
			break;
		} else {
			logger(LOG_INFO, "No match found, max energy %f, time %f\n", round_max_energy,
					round_match_time);
		}

		/* To another iteration:
		 *  - clean the signal
		 *  - restore last n_needle from safe
		 *  - read new samples from input
		 *  - store last n_needle to safe
		 *  - increment n0 with nfft - n_needle
		 */
		memset((void *) haystack_spec, 0, (nfft/2+1) * sizeof(fftw_complex));
		memcpy((void *) haystack_sig, haystack_safe, n_needle * sizeof(double));
		if (read_padded_wav_data(fhaystack, &haystack_hdr, haystack_sig+n_needle, nfft-n_needle) == 0)
			break;
		memcpy((void *) haystack_safe, (void *) (haystack_sig+nfft-n_needle),
				n_needle * sizeof(double));
		n0 += nfft-n_needle;
	}
	if (max_energy > match_treshold) {
		printf("%02d:%02d:%02d.%03d\n", ((int) match_time)/3600, ((int) match_time%3600)/60,
			((int) match_time % 60), ((int)(match_time*1000))%1000);
	} else {
		logger(LOG_INFO, "Finished, no match found. Max energy %f, time %f\n", max_energy,
				match_time);
	}

	return 0;
}




