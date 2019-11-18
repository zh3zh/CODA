
import sys
import os
import regex
import argparse
from collections import defaultdict

def readFastaFile(fastaFile):
    seq = ""
    with open(fastaFile) as f:
        for line in f:
            if(line.startswith('>')):
                continue
            seq += line.strip()
    return seq

def varToMutSeq(refSeq, varString):

    varList = varString.split(',')

    sList = list(refSeq)
    for var in varList:
        pos = int(var[1:len(var)-1])
        type = var[len(var)-1]
        sList[pos] = type
    return ''.join(sList)



def writeMSA(countFile, refSeq, raCutoff, outputFile):

    f = open(countFile, 'r')
    of = open(outputFile, 'w')
    lines = []

    for line in f:
        lines.append(line)

    wtActivity = 1.0
    for line in lines:
        spt = line.split(' ')
        if spt[0]=='0':
            fullCount = float(spt[2])
            partCount = float(spt[3])
            wtActivity = partCount/(partCount+fullCount)


    of.write('>wild type seq\n')
    of.write(refSeq+'\n')
    for line in lines:
        spt = line.split(' ')
        mutNum = int(spt[0])
        if mutNum == 0:
            continue
        fullCount = float(spt[2])
        partCount = float(spt[3])
        if (fullCount+partCount) < minCount:
            continue
        RA = partCount/(partCount+fullCount)/wtActivity
        if RA < raCutoff:
            continue
        of.write('>RA ' + str(round(RA,2)) + " mutNum :" + spt[0] + '\n')
        of.write(varToMutSeq(refSeq, spt[1])+'\n')

    f.close()
    of.close()


parser = None
def argParseInit():
    """
    :rtype: object
    """
    global parser
    parser = argparse.ArgumentParser(description='deduplication and restore real UMI')
    parser.add_argument('--runPath', default='./', help='path to working directory')
    parser.add_argument('--refSeq', required=True, help='reference sequence in fasta format')
    parser.add_argument('--RA', default='0.5', help='Relative Activity cutoff')
    parser.add_argument('--minCount', default='5', help='Relative Activity cutoff')


if __name__ == "__main__":


    argParseInit()
    args = parser.parse_args()
    run_dir = args.runPath
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)


    fasta = readFastaFile(args.refSeq)
    cutoff = float(args.RA)
    minCount = float(args.minCount)

    countFile = run_dir+'/var.count'
    outputFile = run_dir+'/var.msa_RA_'+str(round(cutoff,1))
    writeMSA(countFile, fasta, cutoff,outputFile)
