import os
import regex
import argparse
from collections import defaultdict
import gzip
from multiprocessing import Pool, cpu_count

def reverse_seq(dnaSeq):
    rseq = ""
    for i in range(len(dnaSeq)):
        c = dnaSeq[i]
        if c == 'A':
            rseq = 'T'+rseq
        elif c == 'T':
            rseq = 'A'+rseq
        elif c == 'G':
            rseq = 'C'+rseq
        elif c == 'C':
            rseq = 'G' + rseq
        else:
            rseq = 'N' + rseq
    return rseq


def readFastaFile(fastaFile):
    seq = ""
    with open(fastaFile) as f:
        for line in f:
            if(line.startswith('>')):
                continue
            seq += line.strip()
    return seq


p1 = 'TAAAACGACGGCCAGT'
p2 = 'AATACGACTCACTATAGGGA'
rP2 = 'TCCCTATAGTGAGTCGTATTA'
tag = 'AATCCTGAGTAACTCAAAT'


bc_pattern1 = "(?P<discard_1>." + p1 + "){s<=1}(?P<umi_1>.{15})(?P<discard_2>" + tag + "){s<=1}"\
                             +"(?P<target>(.+)){1}"

bc_pattern2 = "(?P<discard_1>." + p1 + "){s<=1}(?P<umi_1>.{15})(?P<discard_2>" + tag + "){s<=1}"\
                             +"(?P<target>(.+)){1}(?P<discard_3>"+reverse_seq(p2)+"){s<=1}"

umi_cluster = {}
cluster_umi = {}

fullLengthCount = []
partLengthCount = []
runPath = ""
refSeq = ""

def readUmiClusterFromFile(fileName):
    f = open(fileName, 'r')
    for line in f:
        spt = line.split('\t')
        umi = spt[0]
        id = int(spt[1])
        umi_cluster[umi] = spt[1]

        cluster_umi[id] = umi
        fullLengthCount = []*len(umi_cluster)
        partLengthCount = []*len(umi_cluster)


def readMergedFile(idTag):
    fq_read = gzip.open(runPath + "rm_" + idTag + ".fq.gz")
    pt1 = regex.compile(bc_pattern1)
    pt2 = regex.compile(bc_pattern2)
    of = open(runPath + "chop_"+idTag + ".txt", 'w')

    while True:
        try:
            read_lines = [next(fq_read).decode("ascii") for i in range(4)]
        except:
            break
        else:
            m1 = pt1.match(read_lines[1])
            m2 = pt2.match(read_lines[1])
            if m2:
                continue
            if m1:
                umi = m1.group('umi_1')
                if umi in umi_cluster:
                    id = int(umi_cluster[umi])
                    seq = m1.group('target')
                    of.write(str(id)+'\t'+seq+"\n")
    of.close()
    fq_read.close()


def readMergedFileAll():
    idListFile = runPath + "idListRNA"
    with open(idListFile) as f:
        idList = f.read().splitlines()
    p = Pool(cpu_count())
    id = 0
    for y in p.imap(readMergedFile, idList):
        id += 1
    return id


parser = None
def argParseInit():
    """
    :rtype: object
    """
    global parser
    parser = argparse.ArgumentParser(description='generate barcode map')
    parser.add_argument('--runPath', required=True, help='path to working directory')
    parser.add_argument('--refSeq', required=True, help='reference sequence in fasta format')


if __name__ == "__main__":

    argParseInit()
    args = parser.parse_args()
    runPath = args.runPath
    if not os.path.exists(runPath):
        os.makedirs(runPath)
    refSeq = readFastaFile(args.refSeq)
    readUmiClusterFromFile(runPath+'/bc.index')

    readMergedFileAll()
