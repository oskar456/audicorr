#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include <fftw3.h>

#pragma pack(1)
struct	WAV_HEADER
{
	uint8_t		RIFF[4];	    /* RIFF Header = "RIFF"	*/
	uint32_t	ChunkSize;	    /* RIFF Chunk Size = filesize-8 */
	uint8_t		WAVE[4];	    /* WAVE Header = "WAVE"	*/
	uint8_t		fmt[4];	    /* FMT header  = "fmt "	*/
	uint32_t	Subchunk1Size;  /* Size of the fmt chunk				    */
	uint16_t	AudioFormat;    /* Audio format 1=PCM,6=mulaw,7=alaw, 257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM */
	uint16_t	NumOfChan;	    /* Number of channels 1=Mono 2=Sterio		    */
	uint32_t	SampleRate;     /* Sampling Frequency in Hz				    */
	uint32_t	bytesPerSec;    /* bytes per second */
	uint16_t	bytesPerSample; /* 2=16-bit mono, 4=16-bit stereo */
	uint16_t	bitsPerSample;  /* Number of bits per sample      */
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
char needle_fname[] = "needle.wav";
double match_treshold = 0.95;



/**
 * Logger function. Show the message if current verbosity is above
 * logged level.
 *
 * @param levem Message log level
 * @param format printf style format string
 * @returns Whatever printf returns
 */
int logger(enum loglevel level, const char *format, ...) {
	va_list ap, aq;
	int r;
	if (conf_verbosity >= level) {
		va_start(ap, format);
		r=vfprintf(stderr,format, ap);
		va_end(ap);
		return r;
	}
}

void usage(FILE* f) {
	fprintf(f, 
"Audicorr - The audio correlator\n"
"\n"
"Version " VERSION "\n"
"Copyright 2011-12 Ondrej Caletka <ondrej.caletka@gmail.com>\n"
"\n"
"This program is free software; you can redistribute it and/or modify\n"
"it under the terms of the GNU General Public License version 2\n"
"as published by the Free Software Foundation.\n"
"\n"
"Usage:	%s [options] <needle wave file>\n"
"\n"
"Options:\n"
"\t-h --help --usage		Show this help\n"
"\t-V --version			Show program version\n"
"\t-v --verbose			Increase verbosity\n"
"\t-q --quiet			Report only fatal errors\n"
"\t-i --interface <ifname>	\tUse this interface\n"
"\t-s --source <ip or hostname>	Get stream from this source only\n"
"\t-o --output <filename>	\tWrite to file rather than stdout\n"
"\t-u --udponly			Stream is UDP-only, not RTP/UDP\n"
"\t-4 --inet			Force IPv4\n"
"\t-6 --inet6			Force IPv6\n",
	program_invocation_short_name);
}



void parseCmdLine(int argc, char *argv[]) {
	static const struct option longopts[] = {
		{ "help",	no_argument,		NULL,	'h' },
		{ "usage",	no_argument,		NULL,	'h' },
		{ "quiet",	no_argument,		NULL,	'q' },
		{ "verbose",	no_argument,		NULL,	'v' },
		{ "version",	no_argument,		NULL,	'V' },
		{ "interface",	required_argument,	NULL,	'i' },
		{ "source",	required_argument,	NULL,	's' },
		{ "output",	required_argument,	NULL,	'o' },
		{ "udponly",	no_argument,		NULL,	'u' },
		{ "inet",	no_argument,		NULL,	'4' },
		{ "inet6",	no_argument,		NULL,	'6' },
		{ 0,		0,			0,	0   }
	};
	static const char shortopts[] = "hqvVi:s:o:u46";

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
			case 'i':
				conf_interface = strdup(optarg);
				break;
			case 'o':
				conf_output = strdup(optarg);
				break;
			case 's':
				conf_source = strdup(optarg);
				break;
			case 'u':
				conf_udponly = 1;
				break;
			case '4':
				conf_family = AF_INET;
				break;
			case '6':
				conf_family = AF_INET6;
				break;
			default:
				usage(stderr);
				exit(EXIT_FAILURE);
		}
	}

	if (optind >= argc) {
		logger(LOG_FATAL, "Error, no IP address found.\n");
		usage(stderr);
		exit(EXIT_FAILURE);
	}
	conf_IP = strdup(argv[optind]);

	if (++optind < argc) {
		conf_port = strdup(argv[optind]);
	} else {
		conf_port = strdup("1234");
	}

	logger(LOG_DEBUG, "Vebosity: %d, UDPonly: %d, AF: %d\n"
			"Source: %s, Interface: %s,\n"
			"IP: %s, port: %s\n",
			conf_verbosity, conf_udponly, conf_family,
			conf_source, conf_interface,
			conf_IP, conf_port);
}

void check_wave_header(struct WAV_HEADER *wav_hdr) {
	int fail=0;
	if (strncmp((const char *)wav_hdr->RIFF, "RIFF", 4) != 0)
		fail=1;
	if (strncmp((const char *)wav_hdr->WAVE, "WAVE", 4) != 0)
		fail=1;
	if (strncmp((const char *)wav_hdr->fmt, "fmt ", 4) != 0)
		fail=1;
	if (strncmp((const char *)wav_hdr->Subchunk2ID, "data", 4) != 0)
		fail=1;
	if (wav_hdr->AudioFormat != 1) {
		puts("Unsupported data format, only PCM supported");
		fail=1;
	}
	if (wav_hdr->NumOfChan > 1) {
		puts("Warning: Multichannel input - considering only first channel");
	}
	printf("Sample rate: %d\n", wav_hdr->SampleRate);
	printf("Bits per sample: %d\n", wav_hdr->bitsPerSample);
	printf("Bytes per sample: %d\n", wav_hdr->bytesPerSample);
	if (wav_hdr->bitsPerSample != 8 && wav_hdr->bitsPerSample != 16) {
		puts("Unsuported bit width, only 8 or 16 bits supported");
		fail=1;
	}
	if (fail != 0) {
		puts("WAVE header check failed!");
		exit(EXIT_FAILURE);
	}
}

/* Returns next bigger power of 2 */
long next2pow(long n) { 
	double npow;
	npow = log2((double)n);
	return (long) pow(2, ceil(npow));
}

long samplesnumber(struct WAV_HEADER *wav_hdr) {
	return wav_hdr->Subchunk2Size / wav_hdr->bytesPerSample;
}

/*
 * Read at most nmax data samples
 * Returns characters actually read.
 */
long read_wav_data(FILE *wavfile, struct WAV_HEADER *wav_hdr, double *signal, long nmax) {
	int n,r;
	char *sample;

	sample = malloc(wav_hdr->bytesPerSample);
	if (sample == NULL) {
		perror("Allocation failed");
		exit(EXIT_FAILURE);
	}
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

int main(int argc, char *argv[])
{
	FILE *fneedle, *fhaystack = stdin;
	struct WAV_HEADER needle_hdr, haystack_hdr;
	long n_needle, n, max_n, n0;
	fftw_plan plan, backplan;
	fftw_complex *needle_spec, *haystack_spec;
	double *needle_sig, *haystack_sig, *haystack_safe;
	double needle_energy, max_energy, match_time;

	/* Prepare the needle */
	fneedle = fopen(needle_fname, "r");
	if (fneedle == NULL) {
		perror("Cannot open needle file");
		exit(EXIT_FAILURE);
	}
	if (fread(&needle_hdr, sizeof(needle_hdr), 1, fneedle) != 1) {
		perror("Header read failure");
		exit(EXIT_FAILURE);
	}
	check_wave_header(&needle_hdr);
	
	n_needle = samplesnumber(&needle_hdr);
	nfft = (long) pow(2, 18); /* TODO autodetect */

	needle_spec = (fftw_complex*) fftw_malloc((nfft/2+1) * sizeof(fftw_complex));
	if (needle_spec == NULL) {
		perror("Allocation failed");
		exit(EXIT_FAILURE);
	}
	needle_sig = (double *) needle_spec;

	puts("Preparing FFT");
	plan = fftw_plan_dft_r2c_1d(nfft, needle_sig, needle_spec, FFTW_ESTIMATE);
	memset((void *) needle_spec, 0, (nfft/2+1) * sizeof(fftw_complex));
	n_needle = read_wav_data(fneedle, &needle_hdr, needle_sig, nfft);
	needle_energy = 0;
	for (n=0; n<n_needle; n++) {
		needle_energy += needle_sig[n] * needle_sig[n];
	}
	printf("Needle total energy; %f\n", needle_energy);
	puts("Computing FFT");
	fftw_execute(plan);
	fftw_destroy_plan(plan);

	/* Prepare the complex conjugate of needle */
	for(n=0; n<(nfft/2+1); n++) {
		needle_spec[n] = conj(needle_spec[n]);
	}

	printf("Needle length: %ld\nFFT length:%ld\n", n_needle, nfft);

	/* Prepare the haystack */
	fhaystack = fopen("haystack.wav", "r");
	if (fhaystack == NULL) {
		perror("Cannot open haystack file");
		exit(EXIT_FAILURE);
	}
	if (fread(&haystack_hdr, sizeof(haystack_hdr), 1, fhaystack) != 1) {
		perror("Header read failure");
		exit(EXIT_FAILURE);
	}
	check_wave_header(&haystack_hdr);
	haystack_spec = (fftw_complex*) fftw_malloc((nfft/2+1) * sizeof(fftw_complex));
	if (haystack_spec == NULL) {
		perror("Allocation failed");
		exit(EXIT_FAILURE);
	}
	haystack_sig = (double *) haystack_spec;
	puts("Preparing FFT");
	plan = fftw_plan_dft_r2c_1d(nfft, haystack_sig, haystack_spec, FFTW_ESTIMATE);
	puts("Preparing iFFT");
	backplan = fftw_plan_dft_c2r_1d(nfft, haystack_spec, haystack_sig, FFTW_ESTIMATE);
	memset((void *) haystack_spec, 0, (nfft/2+1) * sizeof(fftw_complex));

	/* haystack_safe is used to store last n_needle samples of haystack for
	 * overlaps between iterations
	 */
	haystack_safe = (double *) fftw_malloc(n_needle * sizeof(double));
	if (haystack_safe == NULL) {
		perror("Allocation failed");
		exit(EXIT_FAILURE);
	}
	
	read_wav_data(fhaystack, &haystack_hdr, haystack_sig, nfft);
	memcpy((void *) haystack_safe, (void *) (haystack_sig+nfft-1-n_needle),
			n_needle * sizeof(double));
	n0 = 0;
	while (1) {
		puts("Computing FFT");
		fftw_execute(plan);

		//xcorr(X,Y) = ifft( fft(X) * conj(fft(Y)) )
		for(n=0; n<(nfft/2+1); n++) {
			haystack_spec[n] *= needle_spec[n];
		}
		puts("Computing iFFT");
		fftw_execute(backplan);
		max_n = -1;
		max_energy = 0;
		for(n=0; n<nfft; n++) {
			if (max_energy < haystack_sig[n]) {
				max_energy = haystack_sig[n];
				max_n = n;
			}
		}
		/* Normalize maximum energy to 1.0 */
		max_energy /= nfft * needle_energy; 
		match_time = ((double) max_n+n0) / needle_hdr.SampleRate;

		printf("In sample from %ld to %ld (edge %ld):\n", n0,
				(n0+nfft), (n0+nfft-n_needle));
		printf("In time from %f to %f (edge %f):\n", (double) n0/needle_hdr.SampleRate,
				(double) (n0+nfft) / needle_hdr.SampleRate,
				(double) (n0+nfft-n_needle) / needle_hdr.SampleRate);
		if (max_energy > match_treshold && max_n < (nfft-1-n_needle)) {
			printf("Match found, energy %f, time %f, sample %ld\n", max_energy, match_time, max_n+n0);
			break;
		} else {
			printf("No match found, max energy %f, time %f\n", max_energy, match_time);
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
		if (read_wav_data(fhaystack, &haystack_hdr, haystack_sig+n_needle, nfft-n_needle) == 0)
			break;
		memcpy((void *) haystack_safe, (void *) (haystack_sig+nfft-n_needle),
				n_needle * sizeof(double));
		n0 += nfft-n_needle;
	}

	return 0;
}




