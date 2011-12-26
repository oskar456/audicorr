   //c code to read wav header data and display
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
/* www.cpfaq.blogspot.com */

/* The wav files header information is found on the first 44 bytes.
This structure is used to recover the header information from the wav file
FOR WAVFORMAT */
#pragma pack(1)
typedef struct  WAV_HEADER
{
    uint8_t             RIFF[4];        /* RIFF Header      */ //Magic header
    uint32_t            ChunkSize;      /* RIFF Chunk Size  */
    uint8_t             WAVE[4];        /* WAVE Header      */
    uint8_t             fmt[4];         /* FMT header       */
    uint32_t            Subchunk1Size;  /* Size of the fmt chunk                                */
    uint16_t            AudioFormat;    /* Audio format 1=PCM,6=mulaw,7=alaw, 257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM */
    uint16_t            NumOfChan;      /* Number of channels 1=Mono 2=Sterio                   */
    uint32_t            SamplesPerSec;  /* Sampling Frequency in Hz                             */
    uint32_t            bytesPerSec;    /* bytes per second */
    uint16_t            blockAlign;     /* 2=16-bit mono, 4=16-bit stereo */
    uint16_t            bitsPerSample;  /* Number of bits per sample      */
    uint8_t             Subchunk2ID[4]; /* "data"  string   */
    uint32_t            Subchunk2Size;  /* Sampled data length    */
}wav_hdr; 
#pragma pack()

/* Function prototypes */
int getFileSize(FILE *inFile); 
int main(int argc,char *argv[])
{
    wav_hdr wavHeader;
    FILE *wavFile;
    int headerSize = sizeof(wav_hdr),filelength = 0;
    wavFile = fopen("test.wav","r");
    if(wavFile == NULL)
    {
        printf("Can not able to open wave file\n");
        exit(EXIT_FAILURE);
    }
    fread(&wavHeader,headerSize,1,wavFile);
    filelength = getFileSize(wavFile);
    fclose(wavFile);
    printf("File is %d bytes.\n",filelength);
    printf("RIFF header                           :%c%c%c%c\n",wavHeader.RIFF[0],wavHeader.RIFF[1],wavHeader.RIFF[2],wavHeader.RIFF[3]);
    printf("WAVE header                           :%c%c%c%c\n",wavHeader.WAVE[0],wavHeader.WAVE[1],wavHeader.WAVE[2],wavHeader.WAVE[3]);
    printf("FMT                                   :%c%c%c%c\n",wavHeader.fmt[0],wavHeader.fmt[1],wavHeader.fmt[2],wavHeader.fmt[3]);
    printf("Data size (based on bits used)        :%ld\n",wavHeader.ChunkSize);

    // Display the sampling Rate form the header
    printf("Sampling Rate                         :%ld\n",wavHeader.SamplesPerSec); //Sampling frequency of the wav file
    printf("Number of bits used                   :%d\n",wavHeader.bitsPerSample); //Number of bits used per sample
    printf("Number of channels                    :%d\n",wavHeader.NumOfChan);     //Number of channels (mono=1/sterio=2)
    printf("Number of bytes per second            :%ld\n",wavHeader.bytesPerSec);   //Number of bytes per second
    printf("DATA header                           :%c%c%c%c\n",wavHeader.Subchunk2ID[0],wavHeader.Subchunk2ID[1],wavHeader.Subchunk2ID[2],wavHeader.Subchunk2ID[3]);
    printf("data length                           :%ld\n",wavHeader.Subchunk2Size);
} 
/* find the file size */
int getFileSize(FILE *inFile)
{
 int fileSize = 0;
    fseek(inFile,0,SEEK_END);
    fileSize=ftell(inFile);
    fseek(inFile,0,SEEK_SET);
 return fileSize;
}
