# FASTQ sequences' length harmonizer
======

Harmonize lengths between sequences in a FASTQ (zipped) file.

```no-highlight
usage: fqhar [-h] [-d] [-r] -l <seqlen>  -s <seedlen>      

arguments:                                                 
 -h  show this help message and exit                       
 -l  required sequence length                              
 -s  minimum length of shared string between two reads     
 -d  die when find a sequence which length is smaller than 
     SEQUENCELENGTH                                        
 -r  add reverse and complement reads                      
 -v  (very) verbose output                                 
```

Input is stdin, output is stdout. If `-d` is set and a read with length smaller than <seqlen> is found the process will die, it `-d` is not set, the process will print an error to stderr.

`kseq.h` by [Heng Li](https://github.com/lh3/readfq).