import sys
from sklearn.svm import SVR
from math import log
from math import sqrt
from math import exp
from collections import defaultdict


wcPairs = {}
wcPairs ['AT'] = 1
wcPairs['TA'] = 1
wcPairs['GC'] = 1
wcPairs['CG'] = 1

atgc = ['A', 'T', 'G', 'C']
atgcToInt = {'A':0, 'T':1, 'U':1, 'G':2, 'C':3}



def trainingPredictor(dataFile):
    X = []
    Y = []

    wtRA = 0.0

    f = open(dataFile, 'r')
    for line in f:
        spt = line.split('\t')
        if (spt[0] == '0'):
            wtRA = float(spt[1])
        if (spt[0] == '2'):
            if (spt[8] == '****' or spt[9] == '****' or spt[10] == '****'):
                continue
            r1 = float(spt[8])
            r2 = float(spt[9])
            ra = float(spt[10])

            count = float(spt[5]) + float(spt[6])
            if count < 10.1:
                continue
            if ra > r1 and ra > r2:
                continue
            X.append([r1,r2])
            Y.append(ra)
    f.close()

    #print ("SVR")
    # Fit regression model
    svr = SVR(kernel='rbf', C=2e3, gamma=2.0)
    svr.fit(X,Y)
    predY = svr.predict(X)
    #print("Mean squared error: %.3f "
    #      % mean_squared_error(Y, predY))
    #print('Variance score: %.3f ' % r2_score(Y, predY))
    return svr

def printRegreesionResult(dataFile, outputFile):
    X = []
    Y = []
    f = open(dataFile, 'r')
    for line in f:
        spt = line.split('\t')
        if (spt[0] == '0'):
            wtRA = float(spt[1])
        if (spt[0] == '2'):
            if (spt[8] == '****' or spt[9] == '****' or spt[10] == '****'):
                continue
            r1 = float(spt[8])
            r2 = float(spt[9])
            ra = float(spt[10])

            count = float(spt[5]) + float(spt[6])
            if count < 10.1:
                continue
            if ra < r1 and ra < r2:
                continue
            X.append([r1,r2])
            Y.append(ra)
    f.close()

    #print ("SVR")
    # Fit regression model
    svr = SVR(kernel='rbf', C=1e3, gamma=0.1)
    svr.fit(X,Y)
    predY = svr.predict(X)
    of = open(outputFile,'w')
    for y,py in zip(Y,predY):
        of.write(str(round(y,3))+'\t'+str(round(py,3)))
    of.close()


def generateFeatureMatrix(posA, posB, infoList, svr):
    pairRAList = [99 for i in range(16)]
    tripleMutRAList = [[] for i in range(16)]

    for line in infoList:
        spt = line.strip().split('\t')
        mutNum = int(spt[0])
        if mutNum == 2:
            if spt[8] == "****" or spt[9] == '****' or spt[10] == '****':
                continue
            ra1 = float(spt[8])
            ra2 = float(spt[9])
            predRA = svr.predict([[ra1,ra2]])[0]
            RA = float(spt[10])
            deltaRA = RA - predRA
            base1 = spt[3]
            base2 = spt[4]
            pairType = atgcToInt[base1]*4 + atgcToInt[base2]
            pairRAList[pairType] = deltaRA/(predRA+0.2)
            # pairRAList [pairType+16] = deltaRA

        if mutNum == 3:
            if (spt[7] == '****' or spt[8] == '****' or spt[9] == '****' or spt[12] == '****'):
                continue

            pos1 = int(spt[1])
            pos2 = int(spt[2])
            pos3 = int(spt[3])
            base = spt[4]
            countFull = float(spt[5])
            countPart = float(spt[6])
            if countFull+countPart < 5:
                continue

            ra1 = float(spt[7])
            ra2 = float(spt[8])
            ra3 = float(spt[9])
            RA = float(spt[12])
            rePred = 0.0

            if pos1 == posA and pos2 == posB:
                raThird = ra3
                RA = RA/raThird
                pairType = atgcToInt[base[0]]*4 + atgcToInt[base[1]]
                if raThird < 0.5:
                    continue
                raPred = svr.predict([[ra1, ra2]])[0]
                deltaRA = (RA - raPred)/(raPred + 0.2)
                tripleMutRAList[pairType].append(deltaRA)
            elif pos1 == posA and pos3 == posB:
                raThird = ra2
                RA = RA/raThird
                pairType = atgcToInt[base[0]]*4 + atgcToInt[base[2]]
                if raThird < 0.5:
                    continue
                raPred = svr.predict([[ra1, ra3]])[0]
                deltaRA = (RA - raPred) / (raPred + 0.2)
                tripleMutRAList[pairType].append(deltaRA)
            elif pos2 == posA and pos3 == posB:
                raThird = ra1
                RA = RA/raThird
                pairType = atgcToInt[base[1]]*4 + atgcToInt[base[2]]
                if raThird < 0.5:
                    continue
                raPred = svr.predict([[ra2, ra3]])[0]
                deltaRA = (RA - raPred) / (raPred + 0.2)
                tripleMutRAList[pairType].append(deltaRA)
    for i in range(16):
        if(len(tripleMutRAList[i]) > 0):
            meanDeltaRA = sum(tripleMutRAList[i]) / len(tripleMutRAList[i])
            if len(tripleMutRAList[i]) == 1 and pairRAList[i] == 99 and i == 4:
                continue
            if(pairRAList[i] == 99):
                pairRAList[i] = meanDeltaRA
    return pairRAList


def readFastaFile(fastaFile):
    seq = ""
    with open(fastaFile) as f:
        for line in f:
            if(line.startswith('>')):
                continue
            seq += line.strip()
    return seq


def normalDistribution(x, mean, sd):
    sd2 = sd*sd
    p = (2*3.14159265*sd2)**(-0.5)*exp(-1.0*(x-mean)*(x-mean)/2/sd2)
    return p


'''
paired AT GC  delta RA: mean 1.273 sd 0.77
unpaired delta RA:     mean -0.07 sd 0.21
'''

def pairedDeltaRAToScore(ra, len, mean1, sd1, mean2, sd2):
    if ra > 98:
        return 0
    pRA_unPaired = normalDistribution(ra, mean1, sd1)
    pRA_paired = normalDistribution(ra, mean2, sd2)
    p1 = 0.67/(len-1)
    pRA = pRA_paired*p1 + pRA_unPaired*(1-p1)
    return log(pRA_paired/pRA)

def featureToContactScore(pairFeatures, len, mean1, sd1, mean2, sd2):
    score = 0
    rAT = pairFeatures[1]
    rTA = pairFeatures[4]
    rGC = pairFeatures[11]
    rCG = pairFeatures[14]

    score += pairedDeltaRAToScore(rAT, len, mean1, sd1, mean2, sd2)
    score += pairedDeltaRAToScore(rTA, len, mean1, sd1, mean2, sd2)
    score += pairedDeltaRAToScore(rGC, len, mean1, sd1, mean2, sd2)
    score += pairedDeltaRAToScore(rCG, len, mean1, sd1, mean2, sd2)
    return score


def meanValue(data, startID, endID):
    tot = 0
    for i in range(startID, endID):
        tot += data[i]
    return tot/(endID-startID)


def standardDeviation(data, startID, endID):
    mean = meanValue(data, startID, endID)
    sdTot = 0
    for i in range(startID, endID):
        sdTot += (data[i]-mean)*(data[i]-mean)
    return sqrt(sdTot/(endID-startID))


def predictContact(RAFile, pairRAFile, seqLen, featureOutputFile, outputFile):
    svr = trainingPredictor(RAFile)
    f = open(pairRAFile, 'r')
    featureOut = open(featureOutputFile, 'w')
    features = []
    featureIDA = []
    featureIDB = []
    deltaRAList = []
    matrix = [[0.0 for i in range(seqLen)] for j in range(seqLen)]
    lines = []
    posA = int(0)
    posB = int(0)
    for s in f:
        if s.startswith("#"):
            if len(lines) == 0:
                spt = s[1:len(s)].split(',')
                posA = int(spt[0])
                posB = int(spt[1])
                continue
            pairFeatures = generateFeatureMatrix(posA, posB, lines, svr)
            features.append(pairFeatures)
            featureIDA.append(posA)
            featureIDB.append(posB)

            outList = []
            for i in range(16):
                pf = pairFeatures[i]
                if (pf < 98):
                    outList.append(str(round(pf, 2)))
                else:
                    outList.append("-")
            featureOut.write(str(posA) + "," + str(posB) + "\t" + "\t".join(outList))
            featureOut.write("\n")

            lines = []
            spt = s[1:len(s)].split(',')
            posA = int(spt[0])
            posB = int(spt[1])
        else:
            lines.append(s)
    pairFeatures = generateFeatureMatrix(posA, posB, lines, svr)
    features.append(pairFeatures)
    featureIDA.append(posA)
    featureIDB.append(posB)

    outList = []
    for i in range(16):
        pf = pairFeatures[i]
        if(pf < 98):
            outList.append(str(round(pf,2)))
        else:
            outList.append("-")
    featureOut.write(str(posA)+","+str(posB)+"\t" + "\t".join(outList))
    featureOut.write("\n")
    featureOut.close()

    for i in range(len(features)):
            score1 = features[i][1]
            score2 = features[i][4]
            score3 = features[i][11]
            score4 = features[i][14]
            if score1 < 98:
                deltaRAList.append(score1)
            if score2 < 98:
                deltaRAList.append(score2)
            if score3 < 98:
                deltaRAList.append(score3)
            if score4 < 98:
                deltaRAList.append(score4)

    mean0 = meanValue(deltaRAList, 0, len(deltaRAList))
    sd0 = standardDeviation(deltaRAList, 0, len(deltaRAList))
    
    listA = []
    listB = []
    for i in range(len(deltaRAList)):
        score = deltaRAList[i]
        if score < 3*sd0 + mean0:
            listA.append(score)
        else:
            listB.append(score)
    lenA = len(listA)
    lenB = len(listB)
    print ("lenA: "+str(lenA)+" lenB:"+str(lenB)) 

    mean1 = meanValue(listA, 0, len(listA))
    sd1 = standardDeviation(listA, 0, len(listA))
    mean2 = meanValue(listB, 0, len(listB))
    sd2 = 2*standardDeviation(listB, 0, len(listB))

    print("mean1: " + str(round(mean1, 3)))
    print("sd1: " + str(round(sd1, 3)))
    print("mean2: " + str(round(mean2, 3)))
    print("sd2: " + str(round(sd2, 3)))

    for i in range(len(features)):
        score = featureToContactScore(features[i], seqLen, mean1, sd1, mean2, sd2)
        matrix[featureIDA[i]][featureIDB[i]] = score
        matrix[featureIDB[i]][featureIDA[i]] = score

    f = open(outputFile, 'w')
    for i in range(seqLen):
        list = []
        for j in range(seqLen):
            list.append(str(round(matrix[i][j], 2)))
        f.write('\t'.join(list)+'\n')
    f.close()


if __name__ == "__main__":

    raFile = sys.argv[1]
    pairRAFile = sys.argv[2]
    seq = sys.argv[3]
    seqLen = len(seq)
    featureOutputFile = sys.argv[4]
    outputFile = sys.argv[5]
    predictContact(raFile, pairRAFile, seqLen, featureOutputFile, outputFile)

