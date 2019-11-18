import os
import sys
from multiprocessing import Pool, cpu_count

runPath = ""

def mergeDNA(idTag):
    fqF = runPath + "/d1_" + idTag + ".fq.gz"
    fqR = runPath + "/d2_" + idTag + ".fq.gz"
    dis1 = runPath + "/dis_d1_" + idTag + ".fq.gz"
    dis2 = runPath + "/dis_d2_" + idTag + ".fq.gz"
    mg = runPath + "/dm_" + idTag + ".fq.gz"
    ret = os.system("SeqPrep -f %s -r %s -1 %s -2 %s -s %s" % (fqF, fqR, dis1, dis2, mg))
    return ret


def mergeRNA(idTag):
    fqF = runPath + "/r1_" + idTag + ".fq.gz"
    fqR = runPath + "/r2_" + idTag + ".fq.gz"
    dis1 = runPath + "/dis_r1_" + idTag + ".fq.gz"
    dis2 = runPath + "/dis_r2_" + idTag + ".fq.gz"
    mg = runPath + "/rm_" + idTag + ".fq.gz"
    ret = os.system("SeqPrep -f %s -r %s -1 %s -2 %s -s %s" % (fqF, fqR, dis1, dis2, mg))
    return ret


def mergeAllDNA():
    idListFile = runPath + "idListDNA"
    with open(idListFile) as f:
        idList = f.read().splitlines()
    p = Pool(cpu_count())
    id = 0
    for x in p.imap(mergeDNA, idList):
        id += 1
    return id


def mergeAllRNA():
    idListFile = runPath + "idListRNA"
    with open(idListFile) as f:
        idList = f.read().splitlines()
    p = Pool(cpu_count())
    id = 0
    for y in p.imap(mergeRNA, idList):
        id += 1
    return id


if __name__ == "__main__":
    runPath = sys.argv[1]
    mergeAllDNA()
    mergeAllRNA()
