#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <inttypes.h>
#include <unistd.h>

#include "kseq.h"

#define   OUTWRITES 1024

KSEQ_INIT(gzFile, gzread)

void printusage( )
{
  fprintf( stderr,  "usage: fqhar [-h] [-d] -l <seqlen>  -s <seedlen>           \n" );
  fprintf( stderr,  "\n" );
  fprintf( stderr,  "arguments:                                                 \n" );
  fprintf( stderr,  " -h  show this help message and exit                       \n" );
  fprintf( stderr,  " -l  required sequence length                              \n" );
  fprintf( stderr,  " -s  minimum length of shared string between two reads     \n" );
  fprintf( stderr,  " -d  die when find a sequence which length is smaller than \n" );
  fprintf( stderr,  "     SEQUENCELENGTH                                        \n" );
  fprintf( stderr,  " -v  (very) verbose output                                 \n" );
}

int main( int argc, char ** argv )
{

  gzFile        inf, ouf;
  int           reqlen     =0;
  unsigned int  reqseed    =0;
  int           verbose    =0;
  short int     dieonsmall =0;
  kseq_t       *seq;
  char          c;
  unsigned long rejc       =0, // Rejected reads
                newc       =0, // Produced reads
                oldc       =0; // Reads in original fastq

  while((c = getopt( argc, argv, "hvdl:s:")) >= 0)
    {
      switch(c)
        {
        case 'l': reqlen     = atoi(optarg); break;
        case 's': reqseed    = atoi(optarg); break;
        case 'v': verbose    = 1           ; break;
        case 'd': dieonsmall = 1;          ; break;
        case 'h': printusage()  ; return 1 ; break;
        }
    }

  if(optind == argc || reqseed >= reqlen)
    {
      printusage();
      return -1;
    }

  ouf = gzdopen( fileno(stdout), "wb");
  if(ouf == NULL){ perror("Can't open output file"); exit(1); }

  inf = strcmp(argv[optind], "-") ? gzopen(argv[optind], "r") : gzopen(fileno(stdin), "r");
  if(inf == NULL){ perror("Can't open file"); exit(1); }

  fprintf( stderr,  "Running fqhar with the following parameters \n");
  fprintf( stderr,  "Required length :\t%i\n", reqlen       );
  fprintf( stderr,  "Shared length   :\t%i\n", reqseed        );
  fprintf( stderr,  "Die on small seq:\t%i\n", dieonsmall     );

  seq = kseq_init(inf);

  char        outwrite[OUTWRITES]  ;
  short int   rserror            =0;

  while((kseq_read(seq) >= 0) && (!rserror) )
    {
      ++oldc;
      if(seq->seq.l < reqlen)
        {
          ++rejc;
          if(verbose)
            fprintf( stderr, "Sequence too small:\n\t@%s %s\n\t%s\n\t+\n\t%s\n",
                     seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);
          if(dieonsmall)
            {
              rserror = 1;
              fprintf( stderr, "Process will die.\n" );
            }
        }
      else
        {
          int   beg   =0;
          short lastseq =0;

          for( int sect=0; beg + reqlen <= seq->seq.l; ++sect )
            {
              ++newc;
              int i =1;
              int err;
              // Write seqname and seqinfos
              strncpy( outwrite, "@", 1);
              strncpy( outwrite +i, seq->name.s, seq->name.l);
              i += seq->name.l;
              strncpy( outwrite +i, seq->comment.s, seq->comment.l);
              i += seq->comment.l;
              i += snprintf( outwrite +i, OUTWRITES -i, " %d", sect);
              strncpy( outwrite +i++, "\n", 1);
              strncpy( outwrite +i, seq->seq.s +beg, reqlen);
              i += reqlen;
              strncpy( outwrite +i++, "\n", 1);
              strncpy( outwrite +i++, "+", 1);
              strncpy( outwrite +i++, "\n", 1);
              strncpy( outwrite +i, seq->qual.s +beg, reqlen);
              i += reqlen;
              strncpy( outwrite +i++, "\n", 1);
              strncpy( outwrite +i++, "\0", 1);
              outwrite[i] = '\0';
              if( gzwrite( ouf, outwrite, (unsigned)i-1 ) != i-1 ) error(gzerror(ouf, &err));
              if(!lastseq)
                {
                  if(beg + reqlen +(reqlen - reqseed) < seq->seq.l)
                    {
                      beg += (reqlen - reqseed);
                    }
                  else
                    {
                      if(seq->seq.l - reqlen != beg)
                        {
                          beg = seq->seq.l - reqlen;
                          lastseq = 1;
                        }
                      else
                        beg = seq->seq.l;
                    }
                }
              else
                beg = seq->seq.l;
            }
        }
    }

  kseq_destroy(seq);
  gzclose(inf);
  gzclose(ouf);
  if(!rserror)
    {
      fprintf( stderr, "Done harmonizing data.\n"    );
      fprintf( stderr, "Original read:  %lu\n", oldc );
      fprintf( stderr, "Produced read:  %lu\n", newc );
      fprintf( stderr, "Rejected read:  %lu\n", rejc );
      fprintf( stderr, "Success!\n");
      return 0;
    }
  else
    {
      fprintf( stderr, "Process aborted.\n");
      return -1;
    }
}
