import sys
from collections import defaultdict

input = sys.argv[1]
fasta = sys.argv[2]
seqLen = len(fasta)
output = sys.argv[3]
f = open(input,'r')

lines = []



pairToLineID = {}
for pos1 in range(seqLen):
    for pos2 in range(pos1+1,seqLen):
        key = str(pos1)+','+str(pos2)
        pairToLineID[key] = []


posToPairIndexes = [set() for i in range(seqLen)]


lineID = 0
lines = []
for s in f:
    lines.append(s.strip())
    spt = s.strip().split()
    mutNum = int(spt[0])
    if (mutNum == 1):
        pos = int(spt[1])
        for i in range(pos):
            key = str(i)+','+str(pos)
            pairToLineID[key].append(lineID)
        for i in range(pos+1,seqLen):
            key = str(pos)+','+str(i)
            pairToLineID[key].append(lineID)
    if (mutNum == 2):
        pos1 = int(spt[1])
        pos2 = int(spt[2])
        key = str(pos1)+','+str(pos2)
        pairToLineID[key].append(lineID)
    else :
        posList = []
        for i in range(mutNum):
            posList.append(spt[i+1])
        for i in range(mutNum):
            for j in range(i+1, mutNum):
                key = posList[i]+','+posList[j]
                pairToLineID[key].append(lineID)
    lineID += 1


of = open(output,'w')
for pos1 in range(seqLen):
    for pos2 in range(pos1+1,seqLen):
        key = str(pos1)+','+str(pos2)
        idList = pairToLineID[key]
        of.write("#"+key+"\n")
        for id in idList:
            of.write(lines[id]+'\n')





