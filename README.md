## Synopsis
RSICNV detects copy number variations based on read depths of high coverage whole genome sequencing using the robust segment identification algorithm with negative binomial transformations.

## Code Example
```	
rsicnv                                   #help
rsicnv rsi  -f $REF -b $BAM -o $CNVFILE  #detect CNV from a BAM file
rsicnv plot -f $REF -b $BAM -v $CNVFILE  #plot CNVs in $CNVFILE
```

## Installation
If you have GIT installed, you can download and compile the codes with
```
git clone https://github.com/yhwu/rsicnv.git
cd rsicnv 
make
```
or you can download the files directly with
```
wget https://github.com/yhwu/rsicnv/archive/master.zip -O rsicnv.zip
unzip rsicnv.zip 
cd rsicnv-master
make
```	

Note: you will need g++, gcc compilers, and -lm -lz libs to compile the codes. SAMTOOLS and ALGLIB are included in the download as rsicnv needs their libs. They are not needed afterward.
      
## Tests
This program requires a BAM file and the corresponding reference fasta file. The BAM file must be sorted and both must be indexed by samtools.

```
#1. download a test bam file
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12878/alignment/NA12878.chrom19.ILLUMINA.bwa.CEU.high_coverage.20100311.bam*

#2. download hg19(b37) reference
wget ftp://ftp.sanger.ac.uk/pub/1000genomes/tk2/main_project_reference/human_g1k_v37.fasta.gz
gunzip -c human_g1k_v37.fasta.gz > b37.fasta
samtools faidx b37.fasta

#3. detect CNV
rsicnv rsi  -f b37.fasta -b NA12878.chrom19.ILLUMINA.bwa.CEU.high_coverage.20100311.bam -o test.19.txt

#4. plot CNV
rsicnv plot -f b37.fasta -b NA12878.chrom19.ILLUMINA.bwa.CEU.high_coverage.20100311.bam -v test.19.txt

#5. check overlap with a reference CNV set
checkbp.pl        NA12878hg19_cleaned.txt test.19.txt
checkbp.pl -p 200 NA12878hg19_cleaned.txt test.19.txt
```

## Output fields
<table>
  <tr><td> 1 </td><td> CHROM </td><td> chromosome name </td></tr>
  <tr><td> 2 </td><td> START </td><td> start(lower) position, 1 based </td></tr>
  <tr><td> 3 </td><td> END </td><td>  end(higher) position, 1 based </td></tr>
<tr></tr>
  <tr><td> 4 </td><td> TYPE </td><td>  DEL or DUP </td></tr>
  <tr><td> 5 </td><td> SCORE </td><td> Phred quality score  </td></tr>
  <tr><td> 6 </td><td> LENGTH </td><td> length of variation </td></tr>
  <tr><td> 7 </td><td> CNV_MED </td><td> median of read depths in CNV region </td></tr>
  <tr><td> 7 </td><td> CNV_SD </td><td> standard deviation of read depths in CNV region </td></tr>
  <tr><td> 7 </td><td> NEIGHBOR_MED </td><td> median of read depths before and after CNV region, twice as long as CNV before and twice after </td></tr>
  <tr><td> 7 </td><td> NEIGHBOR_RUNMEANSD </td><td> standard deviation of read depths in neighboring regions </td></tr>
  <tr><td> 7 </td><td> CHR_MED </td><td> median of read depths of whole chromosome </td></tr>
  <tr><td> 7 </td><td> CHR_SD </td><td> standard deviation of read depths of whole chromosome </td></tr>
  <tr><td> 8 </td><td> RP </td><td> number of read pairs overlap with CNV region with at least 50% common to both, read pairs must also be outliers(3sd away) in terms of insert length </td></tr>
  <tr><td> 8 </td><td> Q0 </td><td> fraction of mapq=0 reads in CNV region </td></tr>
  <tr><td> 9 </td><td> METHOD </td><td> RSI </td></tr>
</table>

## System requirement
Linux. Max v_mem is less than 1GB. Running time should be 1 or 2 minutes longer than samtools mpile without the -f option.

## Help
```
[yhwu@debian rsicnv]$ ./rsicnv 
Usage:

1. detect CNV

   rsicnv rsi <options> [-b BAMFILE | -d RDFILE -c RNAME ] -f REFFILE 

Options:
   -m   INT  bin size, default=101
   -q   INT  minimum mapping quality, default=0
   -Q   INT  minimum base quality, default=10
   -cap INT  cap read depth at INT*median, dafault=4
             if INT<0, do not cap read depth
   -NOGC     do not adjust GC content, default=adjust
   -MED      only use median transformation 
   -NB       only use negative binomial transformation (default)
   -s        save read depths before capping and GC adjustment 
   -o   STR  output file, default=rsiout.txt 
   -p   STR  output plot folder, default=cnv_plots
   -np       do not plot CNV

Note:
   The reference file that was used for mapping is needed.
   Input can be given either as a BAM file or a read depth file. 
In the latter case, the chromosome name must be given with -c RNAME. 
A read depth file can be generated with samtools mpileup BAM | cut -f2,4  
keeping only the position and read depth fields. If -c RNAME is given with 
-b BAMFILE, only chromosome RNAME will be processed. By default, samtools 
mpileup does not count duplicated reads, bases with base quality less than 
13, or paired reads far apart. Rsicnv adopts the same rules except the last 
one.

Example:
    rsicnv rsi -f b37.fasta -b $BAM -o all.rsi
    rsicnv rsi -f b37.fasta -b $BAM -c $CHR -o $CHR.rsi
    rsicnv rsi -f b37.fasta -d $RD  -c $CHR -o $CHR.rsi

2. plot CNV

   rsicnv plot <options> [-b BAMFILE | -d RDFILE -c RNAME ] -f REFFILE -v CNVFILE 

Options:
   -p   STR  output folder, default=cnv_plots

Note:
Rsicnv requires gnuplot to plot the figures. The figures are saved in 
./cnv_plots directory in postscript format. If ImageMagic is available, 
the ps files will be converted to the png format. GC contents are not 
adjusted.

Example:
    rsicnv plot -f b37.fasta -b $BAM -v $CNV -p $PLOTFOLDER
```

## AUTHOR
Yinghua Wu, 
Department of Biostistics and Epidemiology, University of Pennsylvania, Philadelphia, PA 19104 

## License
GNU General Public License, version 3.0 (GPLv3)
