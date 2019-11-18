import sys
import os
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


contactList = {}

def isWCPair(a,b):
    x = a+b
    if x == 'AT' or x == 'TA' or x == 'GC' or x == 'CG':
        return True
    return False

def getContactMap(contactFile):
    lines = []
    contactList = defaultdict(int)
    with open(contactFile) as f:
        for line in f:
            lines.append(line)
    ref = lines[1].strip()
    ss = lines[2].strip()
    posInfoList = ['-']*len(ref)
    ssList = list(ss)
    n = len(ss)
    for i in range(n):
        a = ssList[i]
        if (a == ')'):
            for j in range(i-1, -1, -1):
                if (ssList[j] == '('):

                    contactList[str(i)+','+str(j)] = 1
                    contactList[str(j)+','+str(i)] = 1
                    ssList[j] = '.'
                    if isWCPair(ref[j], ref[i]):
                        posInfoList[i] = 'W'
                        posInfoList[j] = 'W'
                    else:
                        posInfoList[i] = 'N'
                        posInfoList[j] = 'N'
                    break
        if (a == ']'):
            for j in range(i-1, -1, -1):
                if (ssList[j] == '['):
                    contactList[str(i)+','+str(j)] = 1
                    contactList[str(j)+','+str(i)] = 1
                    ssList[j] = '.'
                    if isWCPair(ref[j], ref[i]):
                        posInfoList[i] = 'W'
                        posInfoList[j] = 'W'
                    else:
                        posInfoList[i] = 'N'
                        posInfoList[j] = 'N'
                    break
        if (a == '}'):
            for j in range(i-1, -1, -1):
                if (ssList[j] == '{'):
                    contactList[str(i)+','+str(j)] = 1
                    contactList[str(j)+','+str(i)] = 1
                    ssList[j] = '.'
                    if isWCPair(ref[j], ref[i]):
                        posInfoList[i] = 'W'
                        posInfoList[j] = 'W'
                    else:
                        posInfoList[i] = 'N'
                        posInfoList[j] = 'N'
                    break
    # print (contactList)

    return contactList

def readContact(contactFile):
    lines = []
    with open(contactFile) as f:
        for line in f:
            lines.append(line)
    ref = lines[1].strip()
    ss = lines[2].strip()
    posInfoList = ['-']*len(ref)
    ssList = list(ss)
    n = len(ss)
    for i in range(n):
        a = ssList[i]
        if (a == ')'):
            for j in range(i-1, -1, -1):
                if (ssList[j] == '('):
                    contactList[str(i)+','+str(j)] = 1
                    contactList[str(j)+','+str(i)] = 1
                    ssList[j] = '.'
                    if isWCPair(ref[j], ref[i]):
                        posInfoList[i] = 'W'
                        posInfoList[j] = 'W'
                    else:
                        posInfoList[i] = 'N'
                        posInfoList[j] = 'N'
                    break
        if (a == ']'):
            for j in range(i-1, -1, -1):
                if (ssList[j] == '['):
                    contactList[str(i)+','+str(j)] = 1
                    contactList[str(j)+','+str(i)] = 1
                    ssList[j] = '.'
                    if isWCPair(ref[j], ref[i]):
                        posInfoList[i] = 'W'
                        posInfoList[j] = 'W'
                    else:
                        posInfoList[i] = 'N'
                        posInfoList[j] = 'N'
                    break
        if (a == '}'):
            for j in range(i-1, -1, -1):
                if (ssList[j] == '{'):
                    contactList[str(i)+','+str(j)] = 1
                    contactList[str(j)+','+str(i)] = 1
                    ssList[j] = '.'
                    if isWCPair(ref[j], ref[i]):
                        posInfoList[i] = 'W'
                        posInfoList[j] = 'W'
                    else:
                        posInfoList[i] = 'N'
                        posInfoList[j] = 'N'
                    break
    posPairInfo = ''.join(posInfoList)
    return posPairInfo


def varToMutSeq(refSeq, varString):

    varList = varString.split(',')

    sList = list(refSeq)
    for var in varList:
        pos = int(var[1:len(var)-1])
        type = var[len(var)-1]
        sList[pos] = type
    return ''.join(sList)



def readRA(countFile, contactFile, refSeq):

    posPairInfo = '-'*len(refSeq)
    if (contactFile):
        posPairInfo = readContact(contactFile)

    contactList = defaultdict(int)
    if(contactFile):
        contactList = getContactMap(contactFile)

    f = open(countFile, 'r')
    lines = []

    for line in f:
        lines.append(line)

    wtActivity = 1.0
    for line in lines:
        spt = line.split()
        if spt[0]=='0':
            fullCount = float(spt[1])
            partCount = float(spt[2])
            wtActivity = partCount/(partCount+fullCount)
            break

    sVarRAMap = defaultdict(float)

    atgc = ['A', 'T', 'G', 'C']
    baseToInt = {}
    baseToInt['A'] = 0
    baseToInt['T'] = 1
    baseToInt['G'] = 2
    baseToInt['C'] = 3

    n = len(refSeq)
    posRAList = [[-1.0 for i in range(4)] for j in range(n)]

    for i in range(n):
        base = refSeq[i]
        posRAList[i][baseToInt[base]] = 1.0

    print ('0\t'+str(round(wtActivity,3)))

    m1List = []
    for line in lines:
        spt = line.split(' ')
        if spt[0]=='1':
            fullCount = float(spt[2])
            partCount = float(spt[3])
            totCount  = fullCount + partCount
            if (totCount < 5):
                continue
            var = spt[1]
            a = var[0]
            b = var[len(var)-1]

            pos = int(var[1:len(var)-1])
            RA = partCount/totCount/wtActivity
            posRAList[pos][baseToInt[b]] = RA
            sVarRAMap[var] = RA
            outList = ['1', str(pos), b, posPairInfo[pos], str(round(fullCount,1)), str(round(partCount,1)), str(round(RA,2))]
            m1List.append(outList)

    m1List = sorted(m1List, key = lambda x : int(x[1]))
    for m1 in m1List:
        print ('\t'.join(m1))


    pairRAMap = defaultdict(float)

    m2List = []
    for line in lines:
        spt = line.split(' ')
        if spt[0] == '2':
            fullCount = float(spt[2])
            partCount = float(spt[3])
            totCount  = fullCount + partCount
            if (totCount < 5):
                continue
            varList = spt[1].split(',')
            var = varList[0]
            a1 = var[0]
            b1 = var[len(var)-1]
            pos1 = int(var[1:len(var)-1])
            var = varList[1]
            a2 = var[0]
            b2 = var[len(var)-1]
            pos2 = int(var[1:len(var)-1])
            key = str(pos1)+','+str(pos2)


            pairInfo = "none"
            if (key in contactList):
                pairInfo = "pair"

            ra1 = '****'
            ra2 = '****'
            if (varList[0] in sVarRAMap):
                ra1 = str(round(sVarRAMap[varList[0]],2))
            if (varList[1] in sVarRAMap):
                ra2 = str(round(sVarRAMap[varList[1]],2))

            RA = partCount/totCount/wtActivity
            pairRAMap[spt[1]] = RA
            outList = ['2', str(pos1), str(pos2),b1 , b2, spt[2], spt[3].strip(),  pairInfo, ra1, ra2, str(round(RA,2))]
            m2List.append(outList)
            # if (pairInfo == 'none'):
            #     print (str(ra1) + "\t" + str(ra2)+ "\t" + str(RA))

            '''
            if(pairInfo == 'pair'):
                coveredPair += 1
                cPairs.append('\t'.join(outList))

            if (RA >  (ra1 +  ra2)*0.7 and RA > 0.8):
                print ('\t'.join(outList))
            '''
    m2List = sorted(m2List, key = lambda x : int(x[2]))
    m2List = sorted(m2List, key = lambda x : int(x[1]))

    for m2 in m2List:
        print ('\t'.join(m2))
    # print ("total coveredPair: " + str(coveredPair))
    # for p in cPairs:
    #     print (p)

    m3List = []
    for line in lines:
        spt = line.split(' ')
        mutNum = int(spt[0])
        if mutNum > 2:
            fullCount = float(spt[2])
            partCount = float(spt[3])
            totCount  = fullCount + partCount
            if (totCount < 5):
                continue
            varList = spt[1].split(',')
            n = len(varList)
            interPairInfo = ""
            interPairRA = []
            posList = []
            baseList = []
            baseRAList = []
            for i in range(n):
                var1 = varList[i]
                a1 = var1[0]
                b1 = var1[len(var1)-1]
                pos1 = int(var1[1:len(var1)-1])
                posList.append(str(pos1))
                baseList.append(b1)
                sRA = '****'
                if (var1 in sVarRAMap):
                    sRA = str(round(sVarRAMap[var1],2))
                baseRAList.append(sRA)

                for j in range(i+1,n):
                    var2 = varList[j]
                    a2 = var2[0]
                    b2 = var2[len(var2)-1]
                    pos2 = int(var2[1:len(var2)-1])
                    key = str(pos1)+','+str(pos2)
                    pairInfo = "n"
                    if (key in contactList):
                        pairInfo = "p"
                    interPairInfo += pairInfo
                    pRA = '****'
                    key = var1 + ',' + var2
                    if (key in pairRAMap):
                        pRA = str(round(pairRAMap[key],2))
                    interPairRA.append(pRA)

            RA = partCount/totCount/wtActivity

            outList = [spt[0]]
            for pos in posList:
                outList.append(pos)
            outList.append(''.join(baseList))
            outList.append(str(round(fullCount,1)))
            outList.append(str(round(partCount,1)))
            for sRA in baseRAList:
                outList.append(sRA)
            outList.append(interPairInfo)
            outList.append(','.join(interPairRA))
            outList.append(str(round(RA,2)))
            m3List.append(outList)
    m3List = sorted(m3List, key = lambda x : int(x[3]))
    m3List = sorted(m3List, key = lambda x : int(x[2]))
    m3List = sorted(m3List, key = lambda x : int(x[1]))
    m3List = sorted(m3List, key = lambda x : int(x[0]))

    for m3 in m3List:
        print('\t'.join(m3))


parser = None
def argParseInit():
    """
    :rtype: object
    """
    global parser
    parser = argparse.ArgumentParser(description='deduplication and restore real UMI')
    parser.add_argument('--runPath', default='./', help='path to working directory')
    parser.add_argument('--refSeq', required=True, help='reference sequence in fasta format')
    parser.add_argument('--contact', default = None, help='reference sequence in fasta format')
    parser.add_argument('--count', default = None, help='count file')

if __name__ == "__main__":


    argParseInit()
    args = parser.parse_args()
    run_dir = args.runPath
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)

    conFile = args.contact
    fasta = readFastaFile(args.refSeq)
    countFile = args.count
    readRA(countFile, conFile, fasta)
