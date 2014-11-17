#!/bin/bash

#Corecut - cut radio recording using audicorr - The Audio Correlator


INPUT="$1"
OUTPUT="$2"
START_JINGLE="./rozhrani_start.wav"
STOP_JINGLE="./rozhrani_stop.wav"
JINGLE_RATE="8000"

AUDICORR="../../audicorr"
AUDICORR_ARGS="--treshold 0.45"
LAME_ARGS="--abr 160"
TMP_FILE=$(mktemp -u /tmp/correcut.XXXXXX.wav)

if [[ $# -ne 2 ]]; then
	echo "Usage: $0 <input file> <output file>"
	exit
fi

#find the start time
start_time=$( avconv -i "$INPUT" -ac 1 -ar $JINGLE_RATE  -f wav -  </dev/null | 
	$AUDICORR $AUDICORR_ARGS "$START_JINGLE" )

if [[ "x$start_time" == "x" ]]; then
	echo "Cannot find start jingle! Aborting!"
	exit 1
fi

echo "Programme start found at $start_time"

#find the programme length
prog_length=$( avconv -i "$INPUT" -ss $start_time -ac 1 -ar $JINGLE_RATE  -f wav - </dev/null |
		$AUDICORR $AUDICORR_ARGS --end "$STOP_JINGLE" )


if [[ "x$prog_length" == "x" ]]; then
	echo "Cannot find stop jingle! Aborting!"
	exit 1
fi

echo "Programme length found to be $prog_length"

avconv -i "$INPUT" -ss $start_time -t $prog_length -f wav "$TMP_FILE"
normalize --peak "$TMP_FILE"
lame $LAME_ARGS "$TMP_FILE" "$OUTPUT"
rm "$TMP_FILE"
