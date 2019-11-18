This pipeline is used to generate base pairing map from deep mutational sequencing data of RNA ribozyme, with the method introduced in the paper "Precise inference of RNA base-pairing structure by deep mutational scanning and covariation-induced deviation of activity". 

Before running this pipline, make sure these programs were correctly installed.

  SeqPrep (https://github.com/jstjohn/SeqPrep)
  python 2.7.15
  unpigz 2.4
  gsplit 8.29
  pigz 2.4
  samtools 0.0.18
  java 1.7.0_85

Usage:

  bash run.sh $DNAFILE1 $DNAFILE2 $RNAFILE1 $RNAFILE2 $FASTAFILE $OUTPUTPATH

Arguments:
  
  $DNAFILE1: first read DNA sequencing file, in gziped fastq format (fq.gz)
  $DNAFILE2: second read DNA sequencing file, in gziped fastq format (fq.gz)
  $RNAFILE1: first read RNA sequencing file, in gziped fastq format (fq.gz)
  $RNAFILE2: second read RNA sequencing file, in gziped fastq format (fq.gz)
  $FASTAFILE: base sequence file in fasta format (ATGC sequence)
  $OUTPUTPATH: all output files will be written here, empty file folder is recommended

Outputs:
  var.count: uncleaved and cleaved read number of each variant
  var.ra: organized relative activity of each variant
  var.pos.ra: organized relative activity of all mutants of each position pair
  var.msa_RA_0.5: sequence alignment of variants with relative activity higher than 0.5
  pred.mtx: ps score matrix
  pred.ss: 100 predicted seconary struture in bracket format with consensus prediction
