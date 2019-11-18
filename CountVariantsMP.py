
import sys
import os
import argparse
from collections import defaultdict

#barcode id to seq
sBcMap = {}
#barcode id + partPos to seqList
mBcMap = {}
mBcMapWeight = {}
mBcToPosList = {}
revBase = {}
revBase['A'] = 'T'
revBase['T'] = 'A'
revBase['G'] = 'C'
revBase['C'] = 'G'

varSet = {}

varFullMap = defaultdict(int)
varPartMap = defaultdict(int)
runPath = ""

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


def mutSeq(refSeq,readSeq):

    if(len(refSeq) != len(readSeq)):
        return 'length not equal'
    mutPosNum = 0
    varList = []
    for i in range(len(refSeq)):
        a = refSeq[i]
        b = readSeq[i]
        if a != b:
            mutPosNum += 1
            varList.append(a+str(i)+b)
    return str(mutPosNum)+' '+','.join(varList)



def mutSeq2(refSeq,readSeq):

    if(len(refSeq) != len(readSeq)):
        return 'length not equal'
    mutPosNum = 0
    varList = []
    for i in range(len(refSeq)):
        a = refSeq[i]
        b = readSeq[i]
        if a != b:
            mutPosNum += 1
            varList.append(a+str(i))
    return str(mutPosNum)+' '+','.join(varList)


def generateNewBarcodeLib(bcFile, fullSeqLen, shortSeqLen, refSeq, varOut):
    varCount = 0
    f = open(bcFile, 'r')
    of = open(varOut, 'w')
    for line in f:
        if line.startswith('bc'):
            continue
        spt = line.split('\t')
        bcID = spt[0]
        varNum = int(spt[1])
        if varNum == 0:
            continue
        elif varNum == 1:
            seqSpt = spt[2].split('_')
            mutInfo = mutSeq(refSeq, seqSpt[0])
            # mutInfo2 = mutSeq2(refSeq, seqSpt[0])
            sBcMap[bcID] = mutInfo
            infoSpt = mutInfo.split(' ')
            varSet[mutInfo] = int(infoSpt[0])
            varCount += 1
        else:
            for i in range(2,2+varNum):
                seqSpt = spt[i].split('_')
                if len(seqSpt) != 4:
                    print (line)
                    continue
                mutInfo = mutSeq(refSeq, seqSpt[0])
                # mutInfo2 = mutSeq2(refSeq, seqSpt[0])
                infoSpt = mutInfo.split(' ')
                varSet[mutInfo] = int(infoSpt[0])
                varCount += 1
                posList = seqSpt[3].split(',')
                posIndexList = []
                key = bcID
                proportion = float(seqSpt[1])
                for posInfo in posList:
                    revPos = fullSeqLen - int(posInfo[1:len(posInfo)]) -1
                    if revPos >= shortSeqLen:
                        continue

                    posIndexList.append(str(revPos))
                    key += '_'+revBase[posInfo[0]]+str(revPos)
                mBcToPosList[bcID] = ','.join(posIndexList)
                if key in mBcMap:
                    mBcMap[key].append(mutInfo)
                    mBcMapWeight[key].append(proportion)
                else:
                    mBcMap[key] = [mutInfo]
                    mBcMapWeight[key] = [proportion]
    sBcNum = len(sBcMap)
    mBcNum = len(mBcMap)
    mBcCount = [0]*10
    varInfoCount = [0]*100


    of.write('var: ' + str(varCount))
    for key in mBcMap:
        n = len(mBcMap[key])
        mBcCount[n] += 1

    for key in varSet:
        n = varSet[key]
        varInfoCount[n] += 1

    for i in range(10):
        of.write(str(i) + " " + str(varInfoCount[i]) + "\n")
    of.close()


def matchSeq(revRefSeq, readSeq):
    '''
    :param revRefSeq:
    :param readSeq:
    if match full seq: return 0
    if match part seq: return 2
    if not match: return 1
    :return:
    '''
    n1 = len(revRefSeq)
    n2 = len(readSeq)

    n = min(n1,n2)
    diff = 0
    matchLen = 0
    for i in range(n):
        if revRefSeq[i] != readSeq[i]:
            diff += 1
        if diff < 5:
            matchLen += 1

    if matchLen == n1:
        return 0
    elif matchLen < n2 - 8:
        return 1
    else:
        return 2


def readFullShortCounts(fullLen, shortLen, revRefSeq, output, discard):
    df = open(discard, 'w')
    for mutInfo in varSet:
        varFullMap[mutInfo] = 0.0
        varPartMap[mutInfo] = 0.0

    totRead = 0
    matchedRead = 0
    fullLenRead = 0
    fullLenMatchBcRead = 0
    fullLenNotMatchedBcRead = 0
    shortLenRead = 0
    shortLenMatchBcRead = 0
    shortLenNotMatchedBcRead = 0

    idListFile = runPath + "idListRNA"
    with open(idListFile) as f:
        idList = f.read().splitlines()

    for idTag in idList:
        f = open(runPath + "chop_" + idTag + ".txt", 'r')
        for line in f:
            totRead += 1
            spt = line.split('\t')
            diff = 0
            seqLen = len(spt[1])
            testLen = min(seqLen - 10, 30)
            for i in range(testLen):
                if revRefSeq[i] != spt[1][i]:
                    diff += 1
            if diff > 10:
                continue

            matchedRead += 1
            matchType = matchSeq(revRefSeq, spt[1])

            if (seqLen > fullLen - 2 and seqLen < fullLen + 2) or matchType == 0:
                fullLenRead += 1
                seq = spt[1][0:fullLen]
                bcID = spt[0]
                if bcID in sBcMap:
                    mutInfo = sBcMap[bcID]
                    varFullMap[mutInfo] += 1
                    fullLenMatchBcRead += 1
                elif bcID in mBcToPosList:
                    fullLenRead += 1
                    posList = []
                    if len(mBcToPosList[bcID]) > 0:
                        posList = mBcToPosList[bcID].split(',')
                    varList = [bcID]
                    for pos in posList:
                        varList.append(seq[int(pos)] + str(pos))
                    varKey = '_'.join(varList)
                    if varKey in mBcMap:
                        fullLenMatchBcRead += 1
                        vars = mBcMap[varKey]
                        # pList = mBcMapWeight[varKey]
                        for mutInfo in vars:
                            varFullMap[mutInfo] += 1
                    else:
                        fullLenNotMatchedBcRead += 1
            elif (seqLen > shortLen - 2 and seqLen < shortLen + 2) or matchType == 2:
                shortLenRead += 1
                seq = spt[1][0:shortLen]
                bcID = spt[0]
                if bcID in sBcMap:
                    shortLenMatchBcRead += 1
                    mutInfo = sBcMap[bcID]
                    varPartMap[mutInfo] += 1
                elif bcID in mBcToPosList:
                    posList = []
                    if len(mBcToPosList[bcID]) > 0:
                        posList = mBcToPosList[bcID].split(',')
                    varList = [bcID]
                    for pos in posList:
                        varList.append(seq[int(pos)] + str(pos))
                    varKey = '_'.join(varList)
                    if varKey in mBcMap:
                        shortLenMatchBcRead += 1
                        vars = mBcMap[varKey]
                        # pList = mBcMapWeight[varKey]
                        for mutInfo in vars:
                            varPartMap[mutInfo] += 1
                    else:
                        shortLenNotMatchedBcRead += 1

    df.write("total read: " + str(totRead) + " matchedRead: " + str(matchedRead) + "\n")
    df.write("full length: " + str(fullLenRead) + " full length matched: " + str(fullLenMatchBcRead) + "\n")
    df.write("short length: " + str(shortLenRead) + " short length matched: " + str(shortLenMatchBcRead) + "\n")
    df.close()
    print (totRead, matchedRead, fullLenRead, fullLenMatchBcRead, shortLenRead, shortLenMatchBcRead)
    of = open(output,'w')
    for mutInfo in varSet:
        fCount = varFullMap[mutInfo]
        pCount = varPartMap[mutInfo]
        of.write(mutInfo+' '+str(fCount)+' '+str(pCount)+'\n')
    of.close()
    f.close()


def countLen(fileName, revRefSeq):
    f = open(fileName,'r')
    countStat = defaultdict(int)
    count = 0
    for line in f:
        spt = line.split('\t')
        diff = 0
        testLen = min(len(spt[1]), 30)

        for i in range(testLen):
            if revRefSeq[i] != spt[1][i]:
                diff += 1
        if diff > 6:
            continue
        countStat[str(len(spt[1]))] += 1
        count += 1
        if count == 10000:
            break;

    f.close()

    st = sorted(countStat.iteritems(), key = lambda (k,v):(v,k), reverse=True)


    print (st[0])
    print (st[1])
    print (st[3])
    print (st[4])
    print (st[5])
    print (st[6])


parser = None
def argParseInit():
    global parser
    parser = argparse.ArgumentParser(description='count variants')
    parser.add_argument('--runPath', default='./', help='path to working directory')
    parser.add_argument('--refSeq', required=True, help='reference sequence in fasta format')


if __name__ == "__main__":

    argParseInit()
    args = parser.parse_args()
    runPath = args.runPath
    if not os.path.exists(runPath):
        os.makedirs(runPath)

    fasta = readFastaFile(args.refSeq)
    bcFileName = runPath + '/bc2Variants.txt'
    varOutput = runPath + '/var.stat'

    output1 = runPath + '/var.count'
    discard1 = runPath + '/var.discard'

    nRef = len(fasta)
    nShort = 0

    x = 0
    y = 0
    if nRef == 81:
        nShort = 72
        x = 88
        y = 79
    elif nRef == 146:
        nShort = 79
        x = 153
        y = 86
    elif nRef == 140:
        nShort = 31
        x = 147
        y = 34
    elif nRef == 69:
        nShort = 53
        x = 76
        y = 60
    elif nRef == 58:
        nShort = 5
        x = 64
        y = 8

    generateNewBarcodeLib(bcFileName, nRef, nShort, fasta, varOutput)
    readFullShortCounts(x, y, reverse_seq(fasta), output1, discard1)

