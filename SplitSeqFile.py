import sys
import os
import regex
from math import log10
from multiprocessing import Pool, cpu_count
import gzip
import struct

PROG_GUNZIP = "unpigz"
PROG_SPLIT = "gsplit"
PROG_GZIP = "pigz"


def _split_and_compress_fastq(fastq_fn, run_dir, split_prefix, split_lines, suffix_length):
    ret = os.system("%s --stdout %s | %s --numeric-suffixes=1 --suffix-length=%d -l %d --filter='%s > $FILE.fq.gz' - %s" % \
                    (PROG_GUNZIP,
                     fastq_fn,
                     PROG_SPLIT,
                     suffix_length,
                     split_lines,
                     PROG_GZIP,
                     run_dir + '/' + split_prefix))
    if ret != 0:
        return None
    else:
        fn_pattern = regex.compile(regex.escape(split_prefix) + "([0-9]{" + str(suffix_length) + "})\.fq\.gz")
        fn_matches = []
        for fn in os.listdir(run_dir):
            m = fn_pattern.match(fn)
            if m:
                fn_matches.append(m[1])
        fn_matches = sorted(fn_matches, key = lambda x: int(x))
        return [run_dir + '/' + split_prefix + fn_num + ".fq.gz" for fn_num in fn_matches]


def split_and_compress_fastq(run_dir, fastq_read1, fastq_read2, tag):
    if not os.path.isdir(run_dir):
        os.mkdir(run_dir)
    split_reads = 1000000
    split_lines = 4 * split_reads
    suffix_length = 4

    fn1 = _split_and_compress_fastq(fastq_read1, run_dir, tag+"1_", split_lines, suffix_length)
    fn2 = _split_and_compress_fastq(fastq_read2, run_dir, tag+"2_", split_lines, suffix_length)


if __name__ == "__main__":
    run_dir = sys.argv[1]
    fastq_read1 = sys.argv[2]
    fastq_read2 = sys.argv[3]
    tag = sys.argv[4]
    reads_list = split_and_compress_fastq(run_dir, fastq_read1, fastq_read2, tag)

