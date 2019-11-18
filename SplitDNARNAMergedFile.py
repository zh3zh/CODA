import os
import sys
from multiprocessing import Pool, cpu_count
import gzip

runPath = ""

def splitFile(idTag):
    fqFile = runPath + "/tt_" + idTag + ".fq.gz"
    fq_read = gzip.open(fqFile)
    out1 = open(runPath + "/dm_"+idTag + ".fq",  'w')
    out2 = open(runPath + "/rm_" + idTag + ".fq", 'w')
    while True:
        try:
            read_lines = [next(fq_read).decode("ascii") for i in range(4)]
        except:
            break
        else:
            seq = read_lines[1]
            s = seq[len(seq)-4:len(seq)-1]
            if seq[len(seq)-4 : len(seq)-1] == "TTA":
                for k in range(4):
                    out1.write(read_lines[k])
            if seq[len(seq)-4 : len(seq)-1] == "CCC":
                for k in range(4):
                    out2.write(read_lines[k])

    out1.close()
    out2.close()
    os.system("gzip " + runPath + "/dm_"+idTag + ".fq")
    os.system("gzip " + runPath + "/rm_"+idTag + ".fq")
    return 1


def runAll():
    idListFile = runPath + "idListDNA"
    with open(idListFile) as f:
        idList = f.read().splitlines()
    p = Pool(cpu_count())
    id = 0
    for x in p.imap(splitFile, idList):
        id += 1
    return id


if __name__ == "__main__":
    runPath = sys.argv[1]
    runAll()