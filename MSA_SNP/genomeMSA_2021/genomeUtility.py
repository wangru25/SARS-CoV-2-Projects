# -*- coding: utf-8 -*- 
""" 
* @Author: Changchuan Yin 
* @Date: 2020-04-06 08:37:05 
* @Last Modified by:   Rui Wang 
* @Last Modified time: 2020-04-06 08:37:05  
* @Contact: cyin1@uic.edu
* @Contact: wangru25@msu.edu 
"""

import os
import sys
SubDate = sys.argv[1]   # 0405: Represent for submission date: April 05
count = int(sys.argv[2])
numN = int(sys.argv[3])  # 10: Extract sequences with at least 10 descendants. A nonnegative integer.


# Function to check if a list is a 2-D lists
# Find if [0,2,6] is there. 
# [[0, 2, 6], [1, 2, 5], [2, 2, 4], [0, 2, 6]]
def isInList(query,variants):
    query.sort()
    inList = False
    for variant in variants:
        variant.sort()
        if query == variant:
            inList = True
            break
            # print(inList)
    return inList
'''
query = [0, 2, 7] 
datedVariants=[[0, 2, 6], [1, 2, 5], [2, 2, 4], [0, 2, 6]]
print(isInList(query,datedVariants))
'''

def Jaccard(listA,listB):
    unionAB = list(set(listA).union(listB))
    interAB = list(set(listA).intersection(listB))
    dist = 1-len(interAB)/len(unionAB)
    return dist

# Directed Jaccard distance
def DrJaccard(listA,listB):
    
    dist = Jaccard(listA,listB)
    if len(listA) > len(listB):
        dist = dist*(-1)
    return dist

def isInSet(A,B):
    A = set(A)
    B = set(B)
    
    if len(A)<len(B):
        return A.issubset(B)
    else:
        return False

# Check if B is the direct generations (children, grand children, etc) of A
def JaccardNextGens(A,B):
    isNext = False
    dist = DrJaccard(A,B)
    if isInSet(A,B):
        isNext = True
    return dist, isNext
'''
A = [100,150]
B = [100,150,170,180]
dist, isNext = JaccardNextGens(A,B)
print('test2',dist,isNext) 
'''
# Check if B is the direct generations (children, grand children, etc) of A
def JaccardRelatives(A,B):
    isRelative = False
    dist = DrJaccard(A,B)
    isRelative = True
    return dist, isRelative

#Test
'''
A = [100,150]
B = [100,120,170,180]

print('Directed',DrJaccardDr(A,B))
print('Directed',DrDrJaccard(B,A))


B = [100,120]
C = [100,120,150,180]
D = [100,150]
E = [170]
F = [170,200,300]
G = [200,250]
'''

# Get the largest mutations in virus SNP dictionary
def getMaxMutations(virusDict):
    #virusDict ={'CN|WU01|ID1|2020-03-10':[0,2,6,8],'CN|WU01|ID2|2020-02-10':[2,2,4,9,10,12],'CN|WU01|ID3|2020-02-01':[1,2,5],'CN|WU01|ID4|2019-12-01':[0,2,6]}
    lenX = []
    
    for key in virusDict:
        n = len(virusDict[key])
        lenX.append(n)
    return max(lenX)

# Get the date sorted variants list of the viruses, the input is like this:
# virusDict ={'CN|WU01|ID1|2020-03-10':[0,2,6],'CN|WU01|ID2|2020-02-10':[2,2,4],'CN|WU01|ID3|2020-02-01':[1,2,5],'CN|WU01|ID4|2019-12-01':[0,2,6]}
# viruseDict is a dictionary (GenomeMSA_SNP_03022020.py) from cluster Omega MSA using data from GISAID database
def getDateSortedViruses(virusDict):
    #virusDict ={'CN|WU01|ID1|2020-03-10':[0,2,6],'CN|WU01|ID2|2020-02-10':[2,2,4],'CN|WU01|ID3|2020-02-01':[1,2,5],'CN|WU01|ID4|2019-12-01':[0,2,6]}
    
    #Repopulate each dict item into a list after adding a key name 'date' to the date values
    snpRecords = [] #convert the snp dictionary to a list for sorting
    
    for key in virusDict:
        # print(virusDict[key])
        record = {}
        keys= key.split('|')
        record['date'] = keys[3]
        record['variants'] = virusDict[key]
        snpRecords.append(record)
        
    # print('new records:\n',snpRecords)
    # snpRecords= [{'date': '2019-12-01', 'variants': [1, 2, 6]}, {'date': '2019-12-01', 'variants': [1, 2, 6]}, {'date': '2019-12-01', 'variants': [1, 2, 6]}, {'date': '2019-12-01', 'variants': [1, 2, 6]}]
    snpRecords.sort(key = lambda x:x['date'])
    # print('Sorted:\n',snpRecords)
    
    # For the sorted list, we then get the variants according to the date time order.
    datedVariants = []
    for record in snpRecords: 
        datedVariants.append(record['variants'])
    # print(datedVariants)
    return snpRecords, datedVariants

# TEST
'''
virusDict ={'CN|WU01|ID1|2020-03-10':[0,2,6],'CN|WU01|ID2|2020-02-10':[2,2,4],'CN|WU01|ID3|2020-02-01':[1,2,5],'CN|WU01|ID4|2019-12-01':[0,2,6]}
snpRecords, datedVariants = getDateSortedViruses(virusDict)
print(datedVariants)
'''

def getUniqueVariants(datedVariants):
    # datedVariants = [[0, 2, 6], [1, 2, 5], [2, 2, 4], [0, 2, 6]]
    uVariants = []
    for variant in datedVariants:
        if isInList(variant,uVariants) == False:
            uVariants.append(variant)
    # print(uVariants)
    # print(len(uVariants))
    return uVariants
# TEST
'''
datedVariants = [[0, 2, 6], [1, 2, 5], [2, 2, 4], [0, 2, 6]]
print(getUniqueVariants(datedVariants))    
'''

# snpRecords is returned by function: getDateSortedViruses(virusDict)
# snpRecords = [{'date': '2020-03-01', 'variants': [1, 2, 6]}, {'date': '2020-01-01', 'variants': [1, 2, 6]}, {'date': '2019-12-01', 'variants': [1, 2, 6]}, {'date': '2019-12-01', 'variants': [1, 2, 6]}]
def getVariantsByDates(snpRecords, dateS, dateE):
    # Get variants at specific date range
    snpRecordsDT = []
    for record in snpRecords:
        if (record['date'] >= dateS) and (record['date'] <= dateE):
            snpRecordsDT.append(record)
            print('found',record)
    
    return snpRecordsDT

# uniqueDatedVariants = set(datedVariants)
# print(uniqueDatedVariants)
# virusA = {'date': '2019-12-01', 'variants': [1, 2, 6]}

# virusDict = {'A':[100],'A1':[100,120],'B1':[100,120,150,180],'B2':[100,120,150,180,190],'B21':[100,120,150,180,190,200],'C':[100,120],'P':[100,120,150,180],'D':[100,150],'E':[170],'F': [170,200,300],'G': [200,250]}
# snpList =[100,120]
def getAllChildren(snpList,virusDict):
    children = {}
    
    for key,value in virusDict.items():
        dist, isNext = JaccardNextGens(snpList,value)
        # print(dist,isNext)
        if isNext:
            children[key] = value
    
    # childrenList = list(children.items())
    dists = []
    for key ,vlaue in children.items():
        dist = Jaccard(snpList,value)
        dists.append(dist)
        
    # mDist = min(dists)
    # idx =dists.index(mDist)
    # print('Y',idx)
    # firstChild = childrenList[idx] 
    # print(childrenList[idx])
    return children#,firstChild

# virusA = {'date': '2019-12-01', 'variants': [100,120]}
'''
virusA = {'A':[100,120]}
virusDict = {'A':[100],'A1':[100,120],'B1':[100,120,150,180],'B2':[100,120,150,180,190],'B21':[100,120,150,180,190,200],'C':[100,120],'P':[100,120,150,180],'D':[100,150],'E':[170],'F': [170,200,300],'G': [200,250]}
children,firstChild = getAllChildren(virusA,virusDict)
print(children.items())
print(firstChild)
'''
def getSingleMutations(virusDict):
    singles = {}
    
    for key, value in virusDict.items():
        if len(value)==1:
            singles[key]=value   
    # print(singles)
    return singles
'''
virusA = {'A':[100,120]}
virusDict = {'A':[100],'A1':[100,120],'B1':[100,120,150,180],'B2':[100,120,150,180,190],'B21':[100,120,150,180,190,200],'C':[100,120],'P':[100,120,150,180],'D':[100,150],'E':[170],'F': [170,200,300],'G': [200,250]}
singles = getSingleMutations(virusDict)
'''

# print(singles) 
# Get subset of the dictionary which only contains single variants such as 'A':[100], and 'E':[170]
# combine two or more dictionaries using **kwargs
# virusDict  = {**virusDict_US , **virusDict_EU,**virusDict_EA}

#-------------TEST-----------------------------------------------------------
'''
virusDict = {'CN|WU01|ID1|2020-03-10':[0,2,6],'CN|WU01|ID2|2020-02-10':[2,2,4],'CN|WU01|ID3|2020-02-01':[1,2,5],'CN|WU01|ID4|2019-12-01':[0,2,6]}
snpRecords, datedVariants = getDateSortedViruses(virusDict)
print('Total sorted:',snpRecords)

uniqueVariants = getUniqueVariants(datedVariants)
print('Unique mutations',len(uniqueVariants))

#Get virus mutations per date range
# inout: [{'date': '2019-12-01', 'variants': [0, 2, 6]}, {'date': '2020-02-01', 'variants': [1, 2, 5]}, {'date': '2020-02-10', 'variants': [2, 2, 4]}, {'date': '2020-03-10', 'variants': [0, 2, 6]}]
# output: [{'date': '2019-12-01', 'variants': [0, 2, 6]}, {'date': '2020-02-01', 'variants': [1, 2, 5]}]
dateS = '2019-12-01'
dateE = '2020-02-01'
snpRecords2 = getVariantsByDates(snpRecords, dateS, dateE)
print(snpRecords2) 
'''
# This method depends on MSA data (clusteral W)
from Bio import AlignIO
#--------------------------------------------
# what am I doing here??
# if we know the positin in an MSA file, what is the positinin the actual genome?
# for example, position in MSA file entry is 8781, then the correpsonding actual genome (reference genome)

# seqId UK|PHW1|EPI_ISL_413555|2020-02-27 240
# id TW|3|EPI_ISL_411926|2020-01-24 the target genome length: 0
# MSA file GenomesGISAID_SARS-CoV-2_03062020_CN_MSA.txt

# BUT THIS DOES NOT CONSIDER THE GAPS!!
# posMSA is the position on MSA file aligned sequences.
# posMSA is not actual genome sequence position because there are dashes of the MSA positions
# posMSA is returned from alignIO parsing the MSA genomes
# 17858
def mapMSA2RefGenome(msaName, posMSA):
    # msaName = 'GenomesGISAID_2019-nCoV_03022020_EA_MSA.txt'
    # seqId = '2019-nCoV|WH01|NC_045512|2020-01-05'
    # posMSA = 8781
    
    # Fixed RefID, no change
    refSeqId = '2019-nCoV|WH01|NC_045512|2020-01-05'

    # seqId = seqId.strip()
    
    align = AlignIO.read(msaName, "clustal")
    # print('seqId',seqId,posMSA)
    refSeq = ''
    idd = ''
    # seqId is the query seqID, and the idds are those recorde that are in MSA file can align to this id
    
    # UK|PHW1|EPI_ISL_413555|2020-02-27|  CN|P0020|EPI_ISL_413866|2020-02-05
    # UK|PHW1|EPI_ISL_413555|2020-02-27  CN|XN4373-P0039|EPI_ISL_413851|2020-01-30
    # UK|PHW1|EPI_ISL_413555|2020-02-27  CN|P0023|EPI_ISL_413855|2020-02-07

    # print('MSA file',msaName)
    for i in range (0,len(align)):
        idd = align[i].id
        # print('seqId,idd:',seqId, idd)
        if (idd == refSeqId):
            refSeq = align[i].seq
            break
    # print(posMSA)
    # Using MSA positin to get the real position in the genome, starting from 0
    preSeq = refSeq[0: posMSA + 1]  #The sequence before the mismatching position in MSA entries. There are dashed in the sequences (gaps)
    print(preSeq)
    # print('StartSeq',startSeq)
    # Count how many dashes in front of the aligned sequence MSA records
    lastDash = 0
    for ch in preSeq:
        if ch == '-':
            lastDash = lastDash + 1
    # print('id',idd,'the target genome length:',len(seq))
    # print('MSA file',msaName)
    # print('Initial position:',posMSA)
    # print('Geneome seq',seq)
    # print('id',idd,'the reference genome, dashes:',lastDash)
    # Map the MSA position to the real position in the genome
    ntRefGenome = refSeq[posMSA]     #-lastDash] #Python list start from 0
    posRefGenome = posMSA-lastDash+1 #actual position in a genome start from 1, python start from 
    
    #==========================
    # Get the mutated NT on the reference based on the seqId (target virus genome) and position on MSA

    # print('actual position in the genome',posGenome) 
    return [posRefGenome,ntRefGenome]

def mapMSA2VirusGenome(msaName,seqId,posMSA):
    # msaName = 'GenomesGISAID_2019-nCoV_03022020_EA_MSA.txt'
    # seqId = '2019-nCoV|WH01|NC_045512|2020-01-05'
    # posMSA = 8781

    seqId = seqId.strip()
    
    align = AlignIO.read(msaName, "clustal")
    # print('seqId',seqId,posMSA)
    seq = ''
    idd = ''
    # seqId is the query seqID, and the idds are those recorde that are in MSA file can align to this id
    
    # print('MSA file',msaName)
    for i in range (0,len(align)):
        idd = align[i].id
    
        if (idd == seqId):
            seq = align[i].seq
            # print('seqId,idd:',seqId,idd)
            break
        
    # Using MSA positin to get the real position in the genome, starting from 0
    preSeq = seq[0:posMSA+1]  #The sequence before the mismatching position in MSA entries. There are dashed in the sequences (gaps)
    # print('StartSeq',startSeq)
    # Count how many dashes in front of the aligned sequence MSA records
    lastDash = 0
    for ch in preSeq:
        if ch == '-':
            lastDash = lastDash+1
    # print('id',idd,'the target genome, dashes:',lastDash)
    posGenome = posMSA-lastDash+1 #actual position in a genome start from 1, python start from 
    ntGenome = seq[posMSA] 
    return [posGenome,ntGenome]

def mismatchMSA(refSeq,virusSeq):
    cntMatch = 0
    cntTotal = 0
    SNPs = {}
    for i in range(0,len(refSeq)):
        a = refSeq[i]
        b = virusSeq[i]
        # For exact matching only consider non-gaps alignments
        if a != '-' and b != '-':
            if a == b:
                cntMatch=cntMatch+1
            cntTotal=cntTotal+1 #do not count the both are deleted in the MSA file   
        
        # For SNPs, consider both non-gaps and gaps alignments
        if a != b and (a != '-' and b != '-' and a != 'N' and b != 'N'):
            SNPs[i] = a+':'+b
    return SNPs
'''
msaName = 'GenomesGISAID_SARS-CoV-2_US_MSA.txt'
seqId = 'US|WA-UW73|EPI_ISL_415601|2020-03-10'
posMSA = 18060 #28144
# Virus genome0
align = AlignIO.read(msaName, "clustal")
refId =align[36].id
refSeq = align[36].seq
print('reference',refId)
print('reference genome',refSeq[0:250])
# Virus genome0
virusId=align[0].id
virusSeq =align[0].seq
print('virusId',virusId)
print(virusSeq[0:250])
# now compare them and get the mismatches!!
SNPs = mismatchMSA(refSeq,virusSeq)
print(SNPs)    
posMSA = 23403

[posGenomeRef,ntRef] = mapMSA2RefGenome(msaName,posMSA)
print('Reference2',posGenomeRef,ntRef)
[posGenome,ntGenome] = mapMSA2VirusGenome(msaName,seqId,posMSA)    
print('NTGenome2',ntGenome)
'''

def msaNameByStrainId(seqId):
    msaFilePrefix = 'GenomesGISAID_SARS-CoV-2_03062020_'
    #seqId ='US|TX1|EPI_ISL_411956|2020-02-11'
    seqIds = seqId.split('|')
    countryCode = seqIds[0]
    
    if countryCode in ['IT','FR','DE','SE','EN','CH']:
        group = 'EU'
    elif countryCode in ['SG','AU','HT','VN']:
        group = 'SA'
    elif countryCode in ['JP','KR']:
        group = 'EA'
    elif countryCode in ['US','CA']:
        group = 'US'
    elif countryCode in ['CN','TW','HK']:
        group = 'CN'
    else:
        group ='CN'
        
    msaName = msaFilePrefix+group+'_MSA.txt'
    return msaName

# colorsStrains,labelStrains  = getStrainColors(strainNames)
# Examples headers of the genome fasta records
strainNames = ['USA|CA1|ELS1234|2020-02-12','USA|WA1|ELS1234|2020-0223','USA|IL2|ELS1234|2020-0223','USA|CA2|ELS1234|2020-0223','USA|IL1|ELS1234|2020-0223']
def getCityName(strainName):
    # strainName = 'USA|CA1|ELS1234|2020-0223'
    strainNameS = strainName.split('|')
    print(strainName)
    
    labelName = strainNameS[0]+ '|'+strainNameS[1]#strainNameS[3]
    # print(labelName)
    
    cityName = strainNameS[1][0:2] # The first two characters of cities are the same
    # print(cityName)
    return cityName,labelName

# Color a strain sequence in the same color if the strain is in the same city.
def getStrainColors(strainNames):
    cityNames = []
    colorsCity ={}
    cnt = 0
    for strainName in strainNames:
        cityName,labelName = getCityName(strainName)
        # print(cityName,labelName)
        
        if cityName not in cityNames:
            cnt = cnt+1
            cityNames.append(cityName)
            c = 'C{}'.format(cnt%10) # only 10 colors are supported
            colorsCity[cityName]= c
            
    
    # label each strain by a color
    colorsStrains = []
    labelStrains =[]
    for strainName in strainNames:
        cityName,labelName = getCityName(strainName)
        colorsStrains.append(colorsCity[cityName])
        labelStrains.append(labelName)
        
    # print(colorsStrains)
    return colorsStrains,labelStrains

#--------------------------------------------------------------------
# Find the common SNPs among relatives
def getCommon(union):
    common = set(union[0])
    for s in union[1:]:
        common.intersection_update(s)
    common = list(common)
    return common
# TEST
#---------------------------------------------------------------------
'''
n = 2
snps0 = [10,15,16,28]
snps1 = [5,15,16,40]
snps2 = [7,15,16,20]
#snps3 = [8,14,17,30]
#The out should be common=[15,16] in snp0, snp1, and snp2, excluding snp3
union = []
union.append(snps0)
union.append(snps1)
union.append(snps2)
#union.append(snps3)
common = getCommon(union)
print(common)
'''


import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
import networkx as nx
import random
from sklearn.manifold import MDS

def getCloseNext(seqId_Q,virusDict):

    nexts = {}
    snpQ = virusDict[seqId_Q]
    
    for seqId,snp in virusDict.items():
        dist, isNext = JaccardNextGens(snpQ,snp)
        # print(dist,isNext)
        if isNext:
            nexts[seqId] = snp
    
    nextList = list(nexts.items())
    # print('all descendants',nextList)
    dists = []
    
    for seqId ,snp in nexts.items():
        #print('snp',snp)
        #print('virus A',snpA)
        dist = Jaccard(snpQ,snp)
        dists.append(dist)
        
    # print('query',seqId)    
    # print('dist:',dists)
    if len(dists) > 0:  
        mDist = min(dists)
        idx = dists.index(mDist)
        closeNextSeqId = nextList[idx][0]
    else:
        closeNextSeqId=seqId_Q #If there is no cloesed child, it is self, lonely guy
        mDist = 0
    return closeNextSeqId, mDist

# TEST
'''
for seqId,snp in virusDict.items():
closeNextSeqId,mDist = getCloseNext(seqId,virusDict)
print(seqId,closeNextSeqId)
'''
#=============================================================================
# https://jakevdp.github.io/PythonDataScienceHandbook/05.10-manifold-learning.html
# virusDict = {'Y':[500],'A0':[100],'B2':[100,120,150,180,190],'A1':[100,120],'B1':[100,120,150,180],'B21':[100,120,130,150,180,190,200],'C':[100,120,130],'D':[100,150],'E':[100],'F': [170,200,300],'G': [200,250]}
def snpMDS(virusDict):
    """
    Returns the embedded points by MDS.
    Parameters
    ----------
    features: numpy.ndarray
        contains the input feature vectors.
    n_components: int
        number of components to transform the features into

    Returns
    -------
    embedding: numpy.ndarray
        points that the feature vectors have been transformed into
    """
    seqIds = list(virusDict.keys())
    n = len(seqIds)
    # print('Seq', len(seqIds))
    distM = np.zeros([n, n])
    distV = []
    for i in range(0,len(seqIds)):
        seqIdA = seqIds[i]
        posA = virusDict[seqIdA]
        for j in range (0,n):
            seqIdB = seqIds[j]
            posB = virusDict[seqIdB]
            dist = Jaccard(posA ,posB)
            distV.append(dist)
            distM[i,j]= dist
            
    random.seed(9000)
    mds = MDS(n_components = 2, max_iter = 500, dissimilarity = "precomputed", n_jobs = 1)
    coords = mds.fit(distM).embedding_
    mdsCoords = {}
    for idd in range(0,len(seqIds)):
        seqId = seqIds[idd]
        coord = np.array(coords[idd])
        mdsCoords[seqId] = coord
    return mdsCoords

# virusDict = {'Y':[500],'A0':[100],'B2':[100,120,150,180,190],'A1':[100,120],'B1':[100,120,150,180],'B21':[100,120,130,150,180,190,200],'C':[100,120,130],'D':[100,150],'E':[100],'F': [170,200,300],'G': [200,250]}
# mdsCoords = snpMDS(virusDict)
# print('POS',mdsCoords)

# ============================================================================
# when plotting newwork, if two nodes are the same from virusDict, the text are overlapped
# I need to remove the duplicated.
# virusDict = {'Y':[500],'A2':[100,150],'B2':[100,150]}
def getUniqueSNPs(virusDict):  
    snps=[]
    virusDict2={}
    for seqId,snp in virusDict.items():
        if isInList(snp,snps)==False:
            snps.append(snp)
            virusDict2[seqId]=snp  
    # print(virusDict2)    
    return virusDict2
        
#============================================================================
# virusDict = {'Y':[500],'A0':[100],'B2':[100,120,150,180,190],'A1':[100,120],'B1':[100,120,150,180],'B21':[100,120,130,150,180,190,200],'C':[100,120,130],'D':[100,150],'E':[100],'F': [170,200,300],'G': [200,250]}
# pos:  {'Y': array([-0.54224562, -0.59649685]), 'A0': array([0.19024191, 0.45440071]), 'B2': array([ 0.11928144, -0.38031776])...]
def snpNetwork(virusDict, pos = None):
    """
    Input is the virus dictionary (seqId and its snp), and 2-D positions of each viru node
    Pos is a dictionary of coordinates for each virus node
    Output is the network visualization of the virus relationship based on Jaccard distance
    """
    virusDict = getUniqueSNPs(virusDict) 
    
    # 1. get the node colors based on the city of seqID
    nodeNames = []
    for seqId,snp in virusDict.items():
        nodeNames.append(seqId)
    nodeColors,nodeLabels = getStrainColors(nodeNames)

    # 2.get the edges
    snpEdges = []
    for seqId,snp in virusDict.items():
        closeNextSeqId,mDist = getCloseNext(seqId,virusDict)
        edge = (seqId,closeNextSeqId)
        snpEdges.append(edge)
        # print(closeNextSeqId,mDist)
        
    # print(snpEdges)
    
    # 3. create network
    G = nx.DiGraph() # Directed graph
    G.add_edges_from(snpEdges)
    posG = nx.spring_layout(G) # positions
    # print('Node position:',pos)
    
    if pos is None:
        pos = posG
        # print('using default network position')
    else:
        print('Using MDS position') # come from the function paramter
    plt.figure(figsize=(10, 10))
    plt.axis('off')

    nx.draw_networkx_nodes(G, pos, node_color = nodeColors,node_size = 300)#,alpha = 0.50)
    nx.draw_networkx_labels(G, pos, font_size =7,font_family='times',font_weight = 'bold')
    nx.draw_networkx_edges(G, pos, edgelist = snpEdges, edge_color='blue',width = 3,style = 'dashed',arrows=True)
    # nx.draw_networkx(G, with_label = True, node_size=700, node_color = 'red',edge_color = 'blue',width = 3,style = 'dashed',pos = pos)

# mdsCoords = snpMDS(virusDict) 
# snpNetwork(virusDict,pos=mdsCoords)
    
def plotSNPFrequencies(virusDict):
    # virusDict = {'Y':[190],'A0':[100],'B2':[100,120,150,180,190],'A1':[100,120],'B1':[100,120,150,180],'B21':[100,120,130,150,180,190,200],'C':[100,120,130],'D':[100,150],'E':[100],'F': [170,200,300],'G': [200,250]}
    genomeSNPs = {}
    # keys = virusDict.keys()
    for i in range(0,31000): #length of 2019-nCoV
        genomeSNPs[i] = 0
    for seqId,snp in virusDict.items():
        for pos in snp:
            genomeSNPs[pos] = genomeSNPs[pos]+1 #Python starts from 0, whereas genome starts from 1
    # print(genomeSNPs)
    frequencyList = sorted(genomeSNPs.items()) # sorted by key, return a list of tuples
    positions, frequencies = zip(*frequencyList) # unpack a list of pairs into two tuples
    
    frequencies = list(frequencies)
    positions = list(positions)
    snpFreqDict = {}
    
    for i in range(0,len(frequencies)):
        freq = frequencies[i]
        pos = positions[i]
        
        if freq >= 2:
            snpFreqDict[pos] = freq
    plt.figure(figsize = (10,8))
    [markerline, stemlines, baseline] = plt.stem(positions, frequencies, '-.')
    plt.setp(markerline, 'markerfacecolor', 'b')
    plt.setp(baseline, 'color', 'r', 'linewidth', 2)
    plt.xticks(fontsize=14, fontweight = 'bold')
    plt.yticks(fontsize=14, fontweight = 'bold')
    plt.xlabel('Nucleotide position',fontsize = 14,fontweight = 'bold')
    plt.ylabel('Frequency',fontsize = 14,fontweight = 'bold')
    plt.show()

    return snpFreqDict#positions,frequencies

# virusDict = {'Y':[190],'A0':[100],'B2':[100,120,150,180,190],'A1':[100,120],'B1':[100,120,150,180],'B21':[100,120,130,150,180,190,200],'C':[100,120,130],'D':[100,150],'E':[100],'F': [170,200,300],'G': [200,250]}
# positions,frequencies = getSNPFrequencies(virusDict)

# virusDict = {'A':[100],'A1':[100,120],'B1':[100,120,150,180],'B2':[100,120,150,180,190],'B21':[100,120,150,180,190,200],'C':[100,120],'P':[100,120,150,180],'D':[100,150],'E':[170],'F': [170,200,300],'G': [200,250]}
# virusA = {'A':[100,120]}
def getAllNexts(seqId_Q,virusDict):
    children = {}
    # A = virusA['A']
    snpQ = virusDict[seqId_Q]
    for seqId,snp in virusDict.items():
        dist, isNext = JaccardNextGens(snpQ,snp)
        # print(dist,isNext)
        if isNext:
            children[seqId] = snp  
    
    return children

# WHAT IS THIS 
# records=[{'}  ]               
def getAllNextsFromRecords(records):

    for recordQ in records:
        # seqId =  recordQ['seqId']
        snpQ = recordQ['snps']
        cnt=0
        for recordP in records:
            snp = recordP['snps']
            dist, isNext = JaccardNextGens(snpQ,snp)
            # print(dist,isNext)
            if isNext:
                cnt=cnt+1
        recordQ['numNexts'] = cnt
        
    return records

# Get all nexts of recordQ and save the nexts as dictionary
def nextsFromRecords(recordQ,records):
    virusDict = {}
    seqIdQ = recordQ['seqId']
    snpQ = recordQ['snpsRef']
    virusDict[seqIdQ] = snpQ #include itself
    
    for record in records:
        seqId = record['seqId']
        snp = record['snpsRef']
        dist, isNext = JaccardNextGens(snpQ,snp)
        # print(dist,isNext)
        if isNext:
            virusDict[seqId] = snp

    return virusDict

# Convert long records to virus dictionary
def getVirusDictFromRecords(records):
    virusDict = {}
    
    for record in records:
        seqId = record['seqId']
        snp = record['snpsRef']
        virusDict[seqId]=snp
    
    return virusDict

# Convert long records to virus dictionary
def getVirusDictFromRecordsByRegion(records,region):
    virusDict = {}
    #region = 'US'
    for record in records:
        seqId = record['seqId']
        snp = record['snpsRef']
        
        if region in seqId:
            virusDict[seqId]=snp
    
    return virusDict

def mapMutations(genomeRef_PosNTs,genomeVirus_PosNTs):
    #genomeRef_PosNTs={241: 'C', 3037: 'C', 23403: 'A'} 
    #genomeVirus_PosNTs={187: 'T', 2983: 'T', 23349: 'G'}
    
    K1 = list(genomeRef_PosNTs.keys())
    K2 = list(genomeVirus_PosNTs.keys())
    # print(K1)
    
    mutationDict = {}
    for i in range(0,len(K1)):
        r = K1[i]
        NT_r = genomeRef_PosNTs[r]
        v = K2[i]
        NT_v = genomeVirus_PosNTs[v]
        NTInfo = NT_r+'->'+NT_v
        # print(NTInfo)    
        mutationDict[r] = NTInfo
    # print(mutaionDict)
    
    return mutationDict

#============================================================================
# virusDict = {'Y':[500],'A0':[100],'B2':[100,120,150,180,190],'A1':[100,120],'B1':[100,120,150,180],'B21':[100,120,130,150,180,190,200],'C':[100,120,130],'D':[100,150],'E':[100],'F': [170,200,300],'G': [200,250]}
# pos:  {'Y': array([-0.54224562, -0.59649685]), 'A0': array([0.19024191, 0.45440071]), 'B2': array([ 0.11928144, -0.38031776])...]
def snpNetworkNexts(virusDict, pos = None):
    """
    Input is the virus dictionary (seqId and its snp), and 2-D positions of each viru node
    Pos is a dictionary of coordinates for each virus node
    Output is the network visualization of the virus relationship based on Jaccard distance
    """
    virusDict = getUniqueSNPs(virusDict) 
    
    #1. get the node colors based on the city of seqID
    nodeNames = []
    for seqId,snp in virusDict.items():
        nodeNames.append(seqId)
        
    nodeColors,nodeLabels = getStrainColors(nodeNames)
    
    #2.get the edges
    snpEdges = []
    for seqId,snp in virusDict.items():      
        children = getAllNexts(seqId,virusDict)
        for seqIdN,snpN in children.items():
            edge = (seqId,seqIdN) # connect current node to all of its nexts
            snpEdges.append(edge)

    # print(snpEdges)
    
    # 3. create network
    G = nx.DiGraph() # Directed graph
    G.add_edges_from(snpEdges)
    posG = nx.spring_layout(G) # positions
    # print('Node position:',pos)

    if pos is None:
        pos = posG
        # print('using default network position')
    else:
        print('Using MDS position') # come from the function paramter
    plt.figure(figsize=(10, 10))
    plt.axis('off')
    nx.draw_networkx_nodes(G, pos, node_color = nodeColors,node_size = 250)#,alpha = 0.50)
    nx.draw_networkx_labels(G, pos, font_size =7,font_family='times',font_weight = 'bold')
    nx.draw_networkx_edges(G, pos, edgelist = snpEdges, edge_color='blue',width = 1,style = 'dashdot',arrows=True)
    # nx.draw_networkx(G, with_label = True, node_size=700, node_color = 'red',edge_color = 'blue',width = 3,style = 'dashed',pos = pos)

# mdsCoords = snpMDS(virusDict) 
# snpNetwork(virusDict,pos=mdsCoords)
    
def isInGroup(seqId,group):
    isEU = False
    
    if group =='EU':
        EU=['IT','FR','DE','SE','EN','CH','IR','FL','UK']  #CH->Swithwland, SE-Sweden, DE, Germen,
    if group =='EU2':
        EU=['IE','FI','CZ','NL','GE','DK','HU']  #CH->Swithwland, SE-Sweden, DE, Germen,
    if group =='SA':
        EU=['SG','AU','TH','VN','IN','NZ']
    if group =='EA':
        EU=['JP','KR']
    if group =='US':
        EU=['US']
    if group =='US2':
        EU=['MX','CL','BR']
    if group =='CN':
        EU=['CN','TW','HK']
    
    for country in EU:
        #if country in header:
        if seqId.startswith(country):
            isEU = True
    return isEU
'''
seqId = 'NL|NA8|EPI_ISL_415497|2020-03-10'
group = 'EU2'
if isInGroup(seqId,group):
    print('YES')
'''

def getCountryCodes(countries):
    # countries = ['CN_Hangzhou','IT','CN_Beijing','IT','USA','CN_Shangdong']
    # countries = countries.sort()
    codes = []
    for code in countries:
        if 'CN_' in code:
            codes.append('CN')
        else:
            codes.append(code)
    # print(codes) 
    # codes = codes.sort()
    codes.sort()
    codes = str(set(codes))
    yLabel = codes.replace('\'','')
    # print(yLabel)  
    return yLabel   



#=============================================================================
    # -*- coding: utf-8 -*-
"""
# Helper programs
#
# Changchuan Yin
# Dept. of Mathematics, Statistics and Computer Science
# University of Illinois at Chicago
# Chicago, IL 60607
# USA
#
# Email cyin1@uic.edu
# Last update 02/16/2020
#
"""
# DO NOT remove the usage examples
# from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
from Bio import pairwise2
# from Bio.pairwise2 import format_alignment
# https://pypi.python.org/pypi/fuzzysearch/0.3.0#downloads
from fuzzysearch import find_near_matches


#import GenomeGraph as gg
#==============================================================================
#https://en.wikipedia.org/wiki/Jaccard_index
def distJaccard(listA,listB):
    
    unionAB=list(set(listA).union(listB))
    interAB=list(set(listA).intersection(listB))
    dist = 1-len(interAB)/len(unionAB)
    
    return dist

#==============================================================================
def removeNonATCG(seq):
    seq = str(seq)
    for ch in seq:
        if ch not in ['A','T','C','G']:
            seq = seq.replace(ch,'')
    return seq

#==============================================================================
def reverseComplement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
    return "".join([seq_dict[base] for base in reversed(seq)])

#==============================================================================
# To get DNA sequence from a fasta file
# return list of sequences from a fasta file
# outputs: lists of headers and sequences
#==============================================================================
def getDNASequenceFasta(fastaFileName):

    for record in SeqIO.parse(fastaFileName, "fasta"):
        seq=record.seq
        seq=str(seq) 

    return seq

#==============================================================================
# To get DNA sequence from a fasta file
# return list of sequences from a fasta file
# outputs: lists of headers and sequences
#==============================================================================
def getDNASequencesFasta(fastaFileName):
    seqs=[]
    headers=[]
    for record in SeqIO.parse(fastaFileName, "fasta"):
        #print(record.id) #header
        headers.append(record.id)
        seq=record.seq
        seq=str(seq) 
        seqs.append(seq)
    return [headers,seqs]
    
# Usage example:
''' 
fastaFileName='ar_genes_cds.fasta'
fastaFileName =os.path.join('./DNAMobile', fastaFileName)   
seqs= getDNASequencesFasta(fastaFileName)
print('sequence3:',seqs[2])
'''

#==============================================================================
# To write a fasta file of a list of headers and sequences:
# Inputs: list headers, list sequences, and file name
# OutputS: fasta file containing the headers and sequences
#==============================================================================
def writeFastaFileOneRecord(header,seq,handle):
    #handle=open(fileName, 'w')
    handle.write('>'+header+'\n')
    width=70
    for i in range(0, len(seq), width):
        seqT=seq[i:i+width]
        handle.write(seqT + '\n')

#==============================================================================
# To write a fasta file of a list of headers and sequences:
# Inputs: list headers, list sequences, and file name
# OutputS: fasta file containing the headers and sequences
#==============================================================================
def writeFastaFile(headers,seqs,handle):
    #handle=open(fileName, 'w')

    for j in range(0,len(seqs)):
        header=headers[j]
        handle.write('>'+header+'\n')
        seq=seqs[j]

#Usage example:
'''
headers=[]
seqs=[]
headers.append('Gene1')
seqs.append('ATGGCTAAGGCAACAGGTAGGTACAACTTGGTTTCACCTAAAAAGGACCTCGAGAGGGGGCTTGTTTTGAGTGATTTGTGCACGTTTTTAGTTGATCAGACTATCCAGGGGTGGCGGGTGACTTGGGTTGGGATTGAAT')
headers.append('Gene2')
seqs.append('AAACTGCTGACTTCGCTCCTGCATGGTCGATGACAAGGAATTTATTTCCTCATTTATTTCAAAATTCAAATTCTACTATTGAGTCTCCCCTCTGGGCATTACGAGTGATTCTGGCATACCATTATCATCACAAGAACCAACATGGGTT')

fileName='YINSEQ.fasta'
handle = open(fileName,'w')
writeFastaFile(headers,seqs,handle)
'''

#==============================================================
def getComplement(seq):
    #seq=Seq(seq)
    complementSeq=seq.reverse_complement()
    #complementSeq=str(complementSeq)
    return complementSeq

def replaceAll(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j)
    return text

# Usage example:
'''
myText = 'Hello everybody.'
# The dictionary with our key:values,we want to replace 'H' with '|-|''e' with '3' and 'o' with '0'
reps = {'H':'|-|', 'e':'3', 'o':'0'}
txt = replaceAll(myText, reps)
print(txt)   # it prints '|-|3ll0 3v3ryb0dy'
'''
def getSlidingWindows(seq,win,step):
    n=len(seq)
    lastPos=n-win+step
    seqs=[]
    for i in range(0,lastPos,step):
        startPos=i
        endPos=i+win
        seqW=seq[startPos:endPos]
        seqs.append(seqW)
        #print('Seq:',seqW)
    #print('Last seq',seqW)
    return seqs

# Usage example:
'''
seq='AAAAABBBBBCCCCCDDDDDEEEEEFFFF'
win=5
step=2
seqs=getSlidingWindows(seq,win,step)
print('Window seqs:',seqs)

for seq in seqs:
    print('Winoow seq:\n',seq)
'''

#=============================================================================
# Sequence searching and matching
#Search forward subsequence in seqT that matches query sequence seqQ
def searchSequence(seqQ,seqT,maxDist):
    seqQ = seqQ.upper()
    seqT = seqT.upper()
    
    isAlmost = False
    seqFound =''
    distList = []
    startList = []
    endList = []
    minDist = 100
    start = 0
    end = 0
    
    matches = find_near_matches(seqQ, seqT,max_l_dist=maxDist)
    #print('Matches',matches)
    
    numMatches = len(matches)
    if numMatches>0:
        isAlmost = True
    for m in matches:
        startList.append(m.start)
        endList.append(m.end)
        distList.append(m.dist)
    
    minDist = min(distList)
    idx = distList.index(minDist)
    start = startList[idx]
    end = endList[idx]
    
    seqFound = seqT[start:end]
    
    return [isAlmost,seqFound,start,end,minDist] 

# Usage example:
'''
seqQ='AAAAAGGCCAGTCACAATGG' 
seqT='cgtgacgtagtgtgacgtAAAccGGCCAGTCACAATGGggttacgtatcgcgtgtaagtgacgtaAAAAAtttCAGTCACAATGGagtgaccgta'
maxDist=5 #how many can be not matched, When it is zero, isAlmost is the same as isExact
[isAlmost,seqFound,start,end,minDist]  = searchSequence(seqQ,seqT,maxDist)
print('IsAlmost',isAlmost)
print('1.Start position:',start,',End postion:',end)
print('2.Mached target sequence:',seqFound)
print('3.Mismatched NTs:',minDist)
seqFound = seqT[start:end]
print('matched',seqFound.upper())
'''

#pip install biopython 
#from Bio import SeqIO

#Function to get each record in a fasta file (protein or DNA)

def getDNASeqFastaX(seqFasta,idx):
    seqs = []
    ids = []
    for record in SeqIO.parse(seqFasta, "fasta"):
        #print(record.id)  
        idd = record.id
        idd = idd.replace('|','_')
        ids.append(idd)
        seq = record.seq
        seqs.append(seq)
    seq = removeNonATCG(seq) 
    seq = seqs[idx]
    idd = ids[idx]
    return idd,seq

#Test
'''
seqFasta ='SARS-CoV-2_proteome.fasta' #This is proteom  of SARS-CoV-2
k = 0 #First record
idd,seq = getDNASeqFastaX(seqFasta,k) 
'''


#-------------------- ------------------NCBI GenBank--------------------------------------------
# http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_superposition/
#==============================================================================
# To retrieve list of sequences, ids from a genbank file
# Inputs: gbkFileName
# OutputS: GenBankIds,seqs
#==============================================================================

# retrieve list of sequences, ids from a genbank file
def retrieveSequenceGBFile(gbkFileName):
    # https://www.ncbi.nlm.nih.gov/nuccore/?term=txid1914540
    # gbkFilename='HMI.gb'
    inputHandle  = open(gbkFileName, "r")

    genBankIds = []
    seqs = []
    for seqRecord in SeqIO.parse(inputHandle, "genbank"):
        # print('Processing GenBank record', seqRecord.id)
        print('Seq',seqRecord.seq)
        genBankIds.append(seqRecord.id)
        seq = str(seqRecord.seq)
        seqs.append(seq)
        
    print('Seq2:',seqs[2])  
    return [genBankIds,seqs]

# Usage example # https://www.ncbi.nlm.nih.gov/nuccore/FJ349556
'''
gbkFileName = 'HMI.gb'
[genBankIds,seqs] = retrieveSequenceGBFile(gbkFileName)

'''

#==============================================================================
# Function returns description and DNA sequenes for a given genBank ID
# Inputs: a genome accessID,
# OutputS: DNA sequence of the gene ID
#==============================================================================
def getDNASequence(genBankId):
    Entrez.email = 'cyin1@uic.edu'
    handle = None
    seq = ''
    # handle = Entrez.efetch(db = 'nuccore', rettype = 'gbwithparts', id = genBankId,retmode = 'text')
    handle = Entrez.efetch(db = 'nuccore',  rettype = "gb",id = genBankId, retmode = 'text')
    for seqRecord in SeqIO.parse(handle, 'genbank'):
        genBankId = seqRecord.id
        desc = seqRecord.description
        seq = str(seqRecord.seq)
        
    return [desc,seq]  
# rettype = "gb"
# usage example

#==============================================================================
# To retrieve gene genBankId for a given protein Id:
# Inputs: protein Id
# OutputS: gene genBankId
#==============================================================================
def getGenBankByProteinId(ProteinAccessId):
    # accessId = prots[3] 
    reps = {'<':'', '>':'',')':''}
    
    Entrez.email = 'cyin1@uic.edu'
    handle = None
    
    handle = Entrez.efetch(db = 'protein', rettype = "gb", id = ProteinAccessId, retmode = "text")
    for seqRecord in SeqIO.parse(handle, "genbank"):
        # description = seq_record.description
        # print('Record:',seqRecord)
        desc = seqRecord.description
        # print('Desc:',desc)
        aaSeq = seqRecord.seq
        print('proteinSeq retrieved by proteinId:\n',aaSeq)
        
        features = seqRecord.features
        # print('All features',features)
        
        codeBy= ''
        strand = 1
        for f in features:
            tp=f.type
            if tp == 'CDS':
                cdsFeatures = f
                cdsQualifiers = cdsFeatures.qualifiers
                #print('cdsQualifiers:',cdsQualifiers)
                codeBy= cdsQualifiers['coded_by']
            
        codeBy = codeBy[0]
        codeBy = str(codeBy)
        # print('codeBy:',codeBy) 
        codeBy = replaceAll(codeBy, reps)
        
        if codeBy.find('complement')==0: # 0: it is complement, -1, forward
            genBankIdPos = codeBy.split(':')
            ids=genBankIdPos[0].split('(')
            # print('genBankIdPos:',genBankIdPos )
        
            strand = -1
            cdsDirection='complement'
            # print('direction:',cdsDirection)
            
            genBankId = ids[1]
            # print('GenBankId:',genBankId)
            
            pos = genBankIdPos[1]
            startEnd = pos.split('..')
            startPos = startEnd[0]
            # print('Pos start:',startPos)
            
            endPos = startEnd[1]
            # print('Pos end:',endPos)
        
        else:
            strand = 1;      
            cdsDirection = 'forward'
            # print('direction:',cdsDirection)
            
            genBankIdPos = codeBy.split(':')
            genBankId = genBankIdPos[0]
            # print('genBankId:',genBankId )
            
            startEnd = genBankIdPos[1].split('..')
            startPos = startEnd[0]
            # print('Pos start:',startPos)
            
            endPos = startEnd[1]
            # print('Pos end:',endPos)
    
    return [genBankId,cdsDirection,startPos,endPos,desc]  
    #usage
    '''
    proteinAccessId='AAA03550'
    nBankId,strand,startPos,endPos,desc] = getGenBankByProteinId(proteinAccessId)
    print('protein access Id:',proteinAccessId)
    print('gene access Id:',proteinAccessId)
    print('cds direction',strand)
    print('cds starts:',startPos)
    print('cds ends:',endPos)
    print('descs:',desc)
    '''
    # Need to verify it!!!
    def getGenBankCDS(genBankId,strand,startPos,endPos): 
        seq = ''
        cds = ''
        Entrez.email = 'cyin1@uic.edu'
        handle = None

        handle = Entrez.efetch(db='nuccore', rettype='gbwithparts', id=genBankId,retmode='text')
        for seqRecord in SeqIO.parse(handle, 'genbank'):
            # print('Record:',seqRecord)
            # name = seqRecord.description
            seq = seqRecord.seq
            # print('seq len:',len(seq))
            cds=seq[int(startPos)-1:int(endPos)]
            
            if strand == -1:#=='complement':
                cds=cds.reverse_complement()
            
            aa = cds.translate()
            print('proteinSeq translated from genBankId:\n',aa)
        
    cds=str(cds)
    return cds 


#==============================================================================
# Function retrieves gene and protein sequence for a given gene name
# Inputs: a genome accessID, and gene name
# OutputS: proteinId, gene sequence and protein sequence of the gene name
#==============================================================================
def getGeneSequence(genBankId, gene):
    # gene = 'vanRM'
    # genBankId = 'FJ349556'
    startPos = 0
    endPos = 0
    strand = 0
    geneSeq = ''
    AASeq = ''
    # AASeq2 = ''
    proteinId = ''
    # printeinSeq = ''
    # key = 'translation'
    keys= {'gene', 'translation'}# <= set(some_dict)
    Entrez.email = 'cyin1@uic.edu'
    handle=None

    handle = Entrez.efetch(db='nuccore', rettype='gbwithparts', id=genBankId,retmode='text')
    for seqRecord in SeqIO.parse(handle, 'genbank'):
        genBankId = seqRecord.id
        # desc = seqRecord.description
        seq = seqRecord.seq
        
        # print('seqRecord:\n',seqRecord)
        # print('Id:', genBankId)
        # print('desc:',desc) 
        # print('Seq',seq)
        
        for seqFeature in seqRecord.features :
            # print('sequence feature:\n',seqFeature)
            if seqFeature.type == 'CDS':
                if keys <= set(seqFeature.qualifiers): #Check if both gene and transaction exist in feature
                    geneNames = seqFeature.qualifiers['gene']
                    geneName = geneNames[0]
                    geneNames = seqFeature.qualifiers['gene']
                    # product = seqFeature.qualifiers['product']
                if geneName.upper() == gene.upper():
                    loc = seqFeature.location
                    startPos = loc.start
                    endPos = loc.end
                    strand = loc.strand
                    AASeq = seqFeature.qualifiers['translation']
                    proteinId = seqFeature.qualifiers['protein_id']
                
                    # print('ProteinID:',proteinId)
                    # print('Loc:',loc)
                    # print('start',loc.start)
                    # print('end',loc.end)
                    # print('strand',loc.strand) #1=forward strand, -1=complement strand
                    # print('AA_GenBank',AASeq)
                
                if strand == 1:
                    geneSeq=seq[startPos:endPos]
                    # AASeq2 = Seq.translate(geneSeq)
                    # print('AA_TranslatedGene',AASeq2)
                else:
                    geneSeq = seq[startPos:endPos]
                    # print('geneSeq:\n',geneSeq)
                    geneSeq = geneSeq.reverse_complement()
                    # tb=seqFeature.qualifiers['transl_table']
                    # AASeq2=Seq.translate(geneSeq)
                    # print('AA_TranslatedGene',AASeq2)
                break
    return [proteinId,str(geneSeq),str(AASeq)]    

#Usage example 1 
#https://www.ncbi.nlm.nih.gov/nuccore/FJ349556
'''
gene = 'vanXM' #Strand  is 1
gene = 'tnp' #Strand is -1
genBankId = 'FJ349556'
#gene = 'blaZ'
#genBankId = 'BX571856'
[proteinId,geneSeq,AASeq] = getGeneSequence(genBankId,gene)
print('geneSeq: ',geneSeq)
'''
# Usage example 2 (Ebola virus)
'''
genomeGenBankId ='KJ660346'
gene = 'VP24'
[proteinId,geneSeq,AASeq] = getGeneSequence(genomeGenBankId,gene)
print('ProteinID:',proteinId)
print('Gene sequence: ', geneSeq)
print('CDS sequence: ', AASeq)
'''

# Four parameters (2,-1,-1,-1) in alignment:Identical characters are given 2 points, 1 point is deducted for each non-identical character.
# 1 point is deducted when opening a gap, and 1 point is deducted when extending it.

def alignLocal(seqQuery,seqFound):
    # print('Seq1',seqQuery)
    # print('Seq2',seqFound)
    seqQuery = str(seqQuery) #It can take Biopython Seq object
    seqFound = str(seqFound)
    alignments = pairwise2.align.localms(seqQuery,seqFound, 2, -1, -1, -1)

    match = []
    score=0
    for a, b in zip(alignments[0][0],alignments[0][1]):
            if a == b:
                    match.append('|')
            else:
                    match.append('*')
                    score=score+1

    m = "".join(match)
    s = []
    s.append(alignments[0][0] + '\n')
    s.append(m + '\n')
    s.append(alignments[0][1])
    # s = str(s)
    alignedSeqs = "".join(s)
    
    return [alignedSeqs,score]

#=========================================TEST===========================================
# Key viruses for transform study
genes = {
'MN908947.3':'2019-nCoV/Wuhan-Hu-1',          
'MN996532.1':'human-SCoV/RaTG13',
'AY274119':'human-SCoV/Tor2',
'NC_019843.3':'MERS-CoV',
'NC_005831':'human-CoV/NL63',
'AY585229':'human-CoV/OC43',
'AY597011.2':'human-CoV/HKU1',
'KY967357':'Human-CoV/229E/American-3/2015', 
'MG772934.1':'bats-SLCoV/ZXC21',
'MG772933.1':'bats-SLCoV/ZC45', #https://www.ncbi.nlm.nih.gov/nuccore/MG772933
'JX993987':'bats-SLCoV/Rp-Shaanxi2011',
'KC881005':'bats-SLCoV/RsSHC014', #RsSHC014 virus was found later in bat and has same SARS pathogeneity
'KC881006':'bats-SLCoV/Rs3367',   #Rs3367 virus was found later in bat and has same SARS pathogeneity
'KF367457':'bats-SLCoV/WIV1'
}

virusSeqFile ='keyCoVs.fasta'
# virusSeqFile ='POC2.fasta'
# virusSeqFile ='POC.fasta'

# idVirus,virusSeq = getDNASeqFastaX(virusSeqFile,0)
# idHost,hostSeq = getDNASeqFastaX(virusSeqFile,1)

'''
startVirus = 0
startHost = 0

# Align two sequences in the starting positions first

# Case 1.In alignment, viruse is longer than host seq. find if the virus is longer than host sequence using the first 100 bp host sequence to find the position of
# of the corresponding virus, example:
# Virus: TCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATC
# Host   --------------------------CACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATC	

maxDist = 10 #maximum 10 mistmaches are allowed

sHostSeq = hostSeq[0:100]
matchesInVirus = find_near_matches(sHostSeq, virusSeq, max_l_dist=maxDist)

print('matchesInVirus',matchesInVirus)
if len(matchesInVirus)>0: # case 1
    print(matchesInVirus[0].start)
    startVirus = matchesInVirus[0].start
sVirusSeq = virusSeq[0:100]
matchesInHost = find_near_matches(sVirusSeq, hostSeq, max_l_dist=maxDist) 
print('matchesInHost',matchesInHost)
if len(matchesInHost)>0: # case 1
    print(matchesInHost[0].start)
    startHost = matchesInHost[0].start

virusSeq = virusSeq[startVirus:len(virusSeq)]
hostSeq = hostSeq[startHost:len(virusSeq)]
'''
# Align two DNA sequences in the starting positions first
def alignDNASequences(virusSeq, hostSeq, win=100, maxDist = 10):
    
    startVirus = 0
    startHost = 0
    
    # Align two sequences in the starting positions first
    # Case 1.In alignment, viruse is longer than host seq. find if the virus is longer than host sequence using the first 100 bp host sequence to find the position of
    # of the corresponding virus, example:
    # Virus: TCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATC
    # Host   --------------------------CACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATC	
    
    sHostSeq = hostSeq[0:win]
    matchesInVirus = find_near_matches(sHostSeq, virusSeq, max_l_dist=maxDist)
    # print('matchesInVirus',matchesInVirus)
    if len(matchesInVirus)>0: # case 1
        print('Initial virus pos',matchesInVirus[0].start)
        startVirus = matchesInVirus[0].start
        
    sVirusSeq = virusSeq[0:win]
    matchesInHost = find_near_matches(sVirusSeq, hostSeq, max_l_dist=maxDist) 
    # print('matchesInHost',matchesInHost)
    if len(matchesInHost)>0: # case 1
        print('Initial host virus pos',matchesInHost[0].start)
        startHost = matchesInHost[0].start
    
    virusSeq = virusSeq[startVirus:len(virusSeq)]
    hostSeq = hostSeq[startHost:len(hostSeq)]

    return [virusSeq,hostSeq]

#*****************************************************************************
#[virusSeq,hostSeq] = alignDNASequences(virusSeq,hostSeq)

#virusSeqFile2='EPI_ISL_410721.fasta'
#idHost,hostSeq = getDNASeqFastaX(virusSeqFile2,0)

#lenS = len(virusSeq)
# data is obtained by maxDist =6 and win = 21
# i.e., I found these sequence segments of length 21 in 2019-nCoV, and these segments have at least 6 mismatches than the 
# RaTG13 strain (corresponding locations), so these segments are very different from these segments in RaTG13, these
# segments in 2019-nCoV are very mutated from RaTG13 or these segments may come from other species (recombination).

maxDist = 6
# If maxDist is zero, then only find perfect matches.
# If maxDist is 6, then 6 mistaches are allowed for fuzzy searching, the script return those mismatches that are larger than 6 in 21 mer, so these MISMATCHESs are our targets
# if maxDist is 3, then 3 mistaches are alloweed for fuzzy searching.  If the maxDist is smaller, the more records will be returned in this script. 
win = 21
#e = lenS-win+1

#Start positions of high variant segments ()
variants = []
nbrL = win
nbrR = win

'''
for i in range(win,e):
    t = i+win
    seqQ = virusSeq[i:t]
    hostSeqQ = hostSeq[i-nbrL:t+nbrR]
    
    matches = find_near_matches(seqQ, hostSeqQ,max_l_dist=maxDist)
    startList=[]
    endList=[]   
    distList=[]
    
    numMatches = len(matches)
    if numMatches>0:
    for m in matches:
        startList.append(m.start)
        endList.append(m.end)
        distList.append(m.dist)
    
    minDist = min(distList)
    idx = distList.index(minDist)
    
    start = startList[idx]
    end = endList[idx]
    
    if minDist== maxDist: # These MISMATCHESs are our targets
        variants.append(i)
        seqFound = hostSeqQ[start:end]
        print('2019-n',seqQ)
        print('Most Diversed',seqFound)
        
        [alignedSeqs,score]= alignLocal(seqQ,seqFound)
        print(alignedSeqs)
        print('score:',score)

print(variants)   
print('Completed checking', len(variants))

print('Virus1',idVirus)
print('Virus2',idHost)   
'''

#***************************************************Function**************************************************************
def get21M6Variants(virusSeq, hostSeq, win = 21, maxDist = 6, nbr = 250):
    
    # First get the initial positions aligned
    # [virusSeq,hostSeq] = alignDNASequences(virusSeq,hostSeq)
    
    lenS = len(virusSeq)
    # data is obtained by maxDist =6 and win = 21
    # i.e., I found these sequence segments of length 21 in 2019-nCoV, and these segments have at least 6 mismatches than the 
    # RaTG13 strain (corresponding locations), so these segments are very different from these segments in RaTG13, these
    # segments in 2019-nCoV are very mutated from RaTG13 or these segments may come from other species (recombination).
    
    # maxDist = 6
    # If maxDist is zero, then only find perfect matches.
    # If maxDist is 6, then 6 mistaches are allowed for fuzzy searching, the script return those mismatches that are larger than 6 in 21 mer, so these MISMATCHESs are our targets
    # if maxDist is 3, then 3 mistaches are alloweed for fuzzy searching.  If the maxDist is smaller, the more records will be returned in this function. 
    
    # win = 21
    e = lenS - win + 1
    
    #Start positions of high variant segments ()
    variants=[]
    # nbr = 250
    identities = []
    for i in range(0,e):
        t = i + win
        seqQ = virusSeq[i:t]
        
        if i-nbr > 0:
            hostSeqQ = hostSeq[i-nbr:t+nbr]
        else:
            hostSeqQ = hostSeq[0:t+nbr]
            
        matches = find_near_matches(seqQ, hostSeqQ, max_l_dist = maxDist)
        startList = []
        endList = []   
        distList = []
        
        numMatches = len(matches)
        if numMatches > 0:
            for m in matches:
                startList.append(m.start)
                endList.append(m.end)
                distList.append(m.dist)
        
        minDist = min(distList)
        
        identity = round((win-minDist)/win,4)
        identities.append(identity)
        
        idx = distList.index(minDist)
        
        start = startList[idx]
        end = endList[idx]
        
        if minDist == maxDist: # These MISMATCHESs are our targets
            print('Found diversed mismatches',matches)
            variants.append(i)
            seqFound = hostSeqQ[start:end]
            print('2019-n',seqQ)
            print('Most Diversed',seqFound)
            
            [alignedSeqs,score] = alignLocal(seqQ,seqFound)
            print(alignedSeqs)
            print('score:',score)
        
    # print(variants)   
    print('Completed checking', len(variants))
    
    # print('Virus1',idVirus)
    # print('Virus2',idHost)   
    return [variants,identities]


#-----------------------------------------------------TEST------------------------------------------------------
virusSeqFile = 'SARS-CoVs.fasta'
colorLabels = []
positionData = []

#Key 9 SARS viruses for transform study+Pangolin
genes = {
'MN908947.3':'2019-nCoV',          
'MN996532.1':'SLCoV/RaTG13',
'MG772934.1':'SLCoV/ZXC21',
'MG772933.1':'SLCoV/ZC45', #https://www.ncbi.nlm.nih.gov/nuccore/MG772933
'AY274119':'SCoV/Tor2',
'KC881005':'SLCoV/RsSHC014', #RsSHC014 virus was found later in bat and has same SARS pathogeneity
'KC881006':'SLCoV/Rs3367',   #Rs3367 virus was found later in bat and has same SARS pathogeneity
'NC_019843.3':'MERS-CoV',
'JX993987':'SLCoV/Shaanxi2011'
#'KF367457':'SLCoV/WIV1'
}

#idVirus,virusSeq = getDNASeqFastaX(virusSeqFile,0)
#idHost,hostSeq = getDNASeqFastaX(virusSeqFile,4)
#variants= get21M6Variants(virusSeq,hostSeq, win =21, maxDist=6)

#positionData.append(variants)
#gg.eventPlot(idVirus,colorLabels,positionData) 

#***************************************************************************************************
# This method does not depend on MSA data, the result is almost identical with the result of MSA data
'''
idx = 0
idVirus,virusSeq = getDNASeqFastaX(virusSeqFile,idx)
for i in range(0,11):
    if i != idx:
        idHost,hostSeq = getDNASeqFastaX(virusSeqFile,i)
        print('Name',idHost)
        colorLabels.append(idHost)
        #[virusSeq,hostSeq] = alignDNASequences(virusSeq,hostSeq)
        [variants,identities]= get21M6Variants(virusSeq,hostSeq, win =21, maxDist=6)
        positionData.append(variants)

gg.eventPlot(idVirus,colorLabels,positionData) 
'''

#**************************************************************************************************
#This method depends on MSA data (clusteral W)
from Bio import AlignIO
#align = AlignIO.read("SARS-CoVs_ClustalW.txt", "clustal")
'''
alignX = align[1][0:120]

#print(alignX.id)
seqX = alignX.seq
print(seqX)

alignY = align[2][0:120]
#print(alignX.id)
seqY = alignY.seq
print(seqY)
'''

'''
maxDist = 10
matches = find_near_matches(seqQ, hostSeqQ,max_l_dist=maxDist)
print(matches)

'''
def getMismatch(seq,hostSeq):
    cntMatch = 0
    cntTotal = 0
    SNPs = {}
    
    for i in range(0,len(seq)):
        a = seq[i]
        b = hostSeq[i]
        # For exact matching only consider non-gaps alignments
        if a != '-' and b != '-':
            if a==b:
                cntMatch = cntMatch+1
            cntTotal = cntTotal+1 #do not count the both are deleted in the MSA file   
    
        # For SNPs, consider both non-gaps and gaps alignments
        if a != b and (a != '-' and b != '-' and a != 'N' and b != 'N'):
            SNPs[i] = a+':'+b
            
    dist = cntTotal-cntMatch
    # print(cntTotal,cntMatch,dist)
    # identity = 1- round(cntMatch/cntTotal,4)
    return [dist, SNPs]#,identity]

'''
[dist,variants] = getMismatch(seqY,seqX)
print(variants)

keyX=list(variants.keys())
print(keyX)
'''

# Totally based on MSA results
def get21M6VariantsClustalW(virusSeq,hostSeq,win = 21):
    variants =[]
    cntInDel = 0
    lenS = len(virusSeq)
    e = lenS-win+1
    
    for i in range(win,e):
        t = i+win
        seqQ = virusSeq[i:t]
        hostSeqT = hostSeq[i:t]     
            
        [dist, SNPs] = getMismatch(seqQ,hostSeqT)
        # print(identity)
        
        if dist>=6:
            posVirus = i-cntInDel
            # print(posVirus)
            variants.append(posVirus)
        
        if virusSeq[i]=='-':
            cntInDel=cntInDel+1
            
    return variants
# print(round(cntMatch/cntTotal),4)

# print(alignX.id)
# virusSeq = align[1].seq
# hostSeq = align[2].seq
# variants = get21M6VariantsClustalW(virusSeq,hostSeq)
# print(variants)


#=======================================================TEST2===================================================
'''
idx = 0
virusSeq = align[idx].seq
idVirus = align[idx].id
positionData=[]
colorLabels=[]

for i in range(0,9):
    if i != idx:
        hostSeq = align[i].seq
        variants = get21M6VariantsClustalW(virusSeq,hostSeq)
        print(align[i].id)
        colorLabels.append(align[i].id)
        positionData.append(variants)

gg.eventPlot(idVirus,colorLabels,positionData) 
'''

#colorsStrains,labelStrains  = getStrainColors(strainNames)

# Examples headers of the genome fasta records
strainNames = ['USA|CA1|ELS1234|2020-02-12','USA|WA1|ELS1234|2020-0223','USA|IL2|ELS1234|2020-0223','USA|CA2|ELS1234|2020-0223','USA|IL1|ELS1234|2020-0223']

def getCityName(strainName):
    # strainName = 'USA|CA1|ELS1234|2020-0223'
    strainNameS=strainName.split('|')
    # print(strainNameS)
    
    labelName = strainNameS[1]+ '|'+strainNameS[3]
    # print(labelName)
    
    cityName = strainNameS[1][0:2] # The first two characters of cities are the same
    # print(cityName)
    return cityName,labelName

def getStrainColors(strainNames):
    cityNames = []
    colorsCity ={}
    cnt = 0
    for strainName in strainNames:
        cityName,labelName = getCityName(strainName)
        # print(cityName,labelName)
        
        if cityName not in cityNames:
            cnt = cnt + 1
            cityNames.append(cityName)
            c= 'C{}'.format(cnt)
            colorsCity[cityName] = c

    # label each strain by a color
    colorsStrains = []
    labelStrains = []
    for strainName in strainNames:
        cityName,labelName = getCityName(strainName)
        colorsStrains.append(colorsCity[cityName])
        labelStrains.append(labelName)
        
    # print(colorsStrains)
    return colorsStrains,labelStrains

#colorsStrains,labelStrains  = getStrainColors(strainNames)
#print(colorsStrains)
#print(labelStrains)

from Bio import AlignIO

#align = AlignIO.read("GenomesGISAID_MSA_USA.txt", "clustal")
#align = AlignIO.read("GeneomeMSA_SNP_02282020_Japan.txt", "clustal")
#align = AlignIO.read("GeneomeMSA_SNP_02282020_China.txt", "clustal")
#align = AlignIO.read("GeneomeMSA_SNP_02282020_Korea.txt", "clustal")
'''
idx = 47 #China
idx = 17
seqRef = align[idx].seq

refVirus = align[idx].id
print(refVirus)
positionData=[]
colorLabels=[]
viruses ={}

for i in range(0,10):
#for i in range(10,19):
    if i != idx:
        seqVirus = align[i].seq
        virusName = align[i].id
        [dist,SNPs] = getMismatch(seqRef,seqVirus)
        colorLabels.append(align[i].id)
        print('positions',list(SNPs.keys()))

        positions = list(SNPs.keys())
        positionData.append(positions)
    
    if len(positions)>0:
        viruses[virusName]=list(SNPs.keys())
refVirus='2019-nCoV'
gg.eventPlot(refVirus,colorLabels,positionData,lineWidth=4) 

#==============================================================================

idx = 17
seqRef = align[idx].seq

refVirus = align[idx].id
print(refVirus)
positionData=[]
colorLabels=[]
viruses ={}

strainNamesX =[]
for i in range(0,19):

    if i != idx:
        seqVirus = align[i].seq
        virusName = align[i].id
        [dist,SNPs] = getMismatch(seqRef,seqVirus)
        colorLabels.append(align[i].id)
        print('positions',list(SNPs.keys()))

        positions = list(SNPs.keys())
        positionData.append(positions)
        strainNamesX.append(virusName)
    
    if len(positions)>0:
        viruses[virusName]=list(SNPs.keys())

#strainNames = list(viruses.keys())
colorsStrains,labelStrains  = getStrainColors(strainNamesX)
refVirus=''
gg.eventPlot2(refVirus,colorLabels,positionData, colorsStrains,lineWidth=4)

#-----------------------------------------------------------------------------
#Note: virus is a dictionary, the key is the strain name, the value is the position list of mutations (SNPs)
print(viruses)

# MDS Plotting
import numpy as np
from matplotlib import pyplot as plt
from sklearn import manifold

strainNames = list(viruses.keys())
n = len(strainNames)
distM=np.zeros([n, n]) 
distV=[]

for i in range(0,len(strainNames)):
    nameA = strainNames[i]
    posA = viruses[nameA]
    for j in range (0,n):
        nameB = strainNames[j]
        posB = viruses[nameB]
        dist =  distJaccard(posA ,posB)
        distV.append(dist)
        distM[i,j]=dist
#-----------------------------------------------------------------------------
import random
random.seed(9000)

mds = manifold.MDS(n_components=2, max_iter=500, dissimilarity="precomputed", n_jobs = 1)
pos = mds.fit(distM).embedding_

#model = manifold.TSNE(metric='precomputed')
#pos = model.fit_transform(distM) 
#not good as MDS

#import umap
#U = umap.UMAP(metric='precomputed')
#pos= U.fit_transform(distM)
#NOT as good as MDS

shortProteinNames=[]
# Plot the points
fig = plt.figure(figsize = (10,10))
m=len(strainNames)
colors = ['C{}'.format(i%10) for i in range(m)]

#markerSizes=m*[2] #Illustration for realative lengths of amino acid sequences
#v=[0.02*j for j in markerSizes]

colorsStrains,labelStrains  = getStrainColors(strainNames)

for i in range(len(pos)):
    plt.plot(pos[i, 0], pos[i, 1],marker='o',c=colorsStrains[i],markersize=20)
    plt.text(pos[i, 0]+0.02, pos[i, 1]+0.02, labelStrains[i], fontsize=10,fontweight = 'bold') 

plt.xticks(fontsize=14, fontweight = 'bold')
plt.yticks(fontsize=14, fontweight = 'bold')
plt.xlabel('Similitude longitude', labelpad=5,fontsize=12,fontweight = 'bold')
plt.ylabel('Similitude latitude',labelpad=5,fontsize=12,fontweight = 'bold')
plt.show()
print('Completed')

'''
'''
positionData=[]
colorLabels=[]
for i in range(9,19):
    if i != idx:
        seqVirus = align[i].seq
        [dist,SNPs] = getMismatch(seqRef,seqVirus)
        #print(SNPs)
        # print(align[i].id)
        colorLabels.append(align[i].id)
        positionData.append(list(SNPs.keys()))

gg.eventPlot(refVirus,colorLabels,positionData,lineWidth=4) 
'''

#https://www.biostars.org/p/382859/
#https://biopython.org/DIST/docs/api/Bio.Align.Applications._ClustalOmega.ClustalOmegaCommandline-class.html
#===============================================================
#Get the reference genome index and total number of records in the MSA file
def getReferenceMSA(msaName):
    # align = AlignIO.read("GenomesGISAID_MSA_Japan.txt", "clustal")
    align = AlignIO.read(msaName, "clustal")
    # print(len(align))
    N = len(align) #Total sequences in the alignment files
    
    idx = -1
    for record in align:
        recordId = (record.id)
        idx = idx+1
        if '2019-nCoV' in recordId:
            print(recordId,idx)
            break
    refIdx = idx
    # print('Reference',align[idx],idx)
    return [N,refIdx]
'''
msaName = 'GenomesGISAID_MSA_Japan.txt'
[N,refIdx] = getReferenceMSA(msaName)
print(N,refIdx)

#==========================================================================
    
import Gene2020 as gene
seq='GTTGCAACTGCAGAAGCTGAA' #Wuhan GenomeAnalysis_02102020.py
seq='GTTGCTACAGCTCACAGCGAG' #SARS

[proteinSeq0,proteinSeq1,proteinSeq2] = gene.translateDNA3FramesX(seq)
print(proteinSeq0,proteinSeq1,proteinSeq2)

mutations={8782: 'C->T', 28144: 'T->C'}

refFasta='NC_045512.fasta'
idx = 0
idd,seq = getDNASeqFastaX(refFasta,idx)
print(idd)

print('Length genome',len(seq))
pos = 8782-1 #C-> Nonsense mutation ;https://www.ncbi.nlm.nih.gov/protein/1796318597?report=fasta
pos = 28144-1 # T->C #what is this?? https://www.uniprot.org/proteomes/?query=taxonomy:694009
#https://www.uniprot.org/uniprot/?query=proteome:UP000000354+AND+proteomecomponent:%22Genome%22&sort=score
#pos = 26144-1 #(G->T)
#NS protein
#https://www.ncbi.nlm.nih.gov/protein/QHR63301.1?report=genbank&log$=prottop&blast_rank=1&RID=64W20U8J016

pos = 23403-1
pos = 14408-1 #C-T
pos = 3037-1 #C->T
#241: 'C->T'
pos = 241-1

# {11083: 'G->T'}
'''
#----------------------------------------------------------------------------------
codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }


def translateDNA3Frames(seq):  
    proteinSeq0 = ''
    proteinSeq1 = ''
    proteinSeq2 = ''
    #seq = 'ATAACAAACAGCTTA'
    # start = seq.find('ATG')
    # seqStart = seq[int(start):]
    # stop = seqStart.find('TAA')
    # cds = str(seqStart[:int(stop)+3])
    
    for i in range(0,len(seq),3):
        codon0 = seq[i:i+3]
        codon1 = seq[i+1:i+4]
        codon2 = seq[i+2:i+5]
        
        if len(codon0)==3:
            proteinSeq0 += codontable[codon0]
            
        if len(codon1)==3:
            proteinSeq1 += codontable[codon1]
        
        if len(codon2)==3:
            proteinSeq2 += codontable[codon2]
        
    proteinSeq = ''    
    if '_' not in proteinSeq0: 
        proteinSeq = proteinSeq0
    
    if '_' not in proteinSeq1: 
        proteinSeq = proteinSeq1
        
    if '_' not in proteinSeq2: 
        proteinSeq = proteinSeq2  
    # print(proteinSeq0)
    # print(proteinSeq1)
    # print(proteinSeq2)
    return proteinSeq


# Removed these ORFs that have "stop codon"
def goodORFs(ORFs):
    #ORFs=['FYENAFLPFA', 'FMKMPFYLL', 'L_KCLFTFC']
    goodORFs = []
    for orf in ORFs:
        if '_' not in orf:
            goodORFs.append(orf)
    return goodORFs

def translateDNA3FramesX(seq):  
    proteinSeq0 = ''
    proteinSeq1 = ''
    proteinSeq2 = ''
    # seq = 'ATAACAAACAGCTTA'
    # start = seq.find('ATG')
    # seqStart = seq[int(start):]
    # stop = seqStart.find('TAA')
    # cds = str(seqStart[:int(stop)+3])
    codon0s = ''
    codon1s = ''
    codon2s = ''
    for i in range(0,len(seq),3):
        codon0 = seq[i:i+3]
        codon0s = codon0s + codon0+ '|'
        codon1 = seq[i+1:i+4]
        codon1s = codon1s + codon1+ '|'
        codon2 = seq[i+2:i+5]
        codon2s = codon2s + codon2+ '|'
        
        if len(codon0) == 3:
            proteinSeq0 += codontable[codon0]
            
        if len(codon1) == 3:
            proteinSeq1 += codontable[codon1]
        
        if len(codon2) == 3:
            proteinSeq2 += codontable[codon2]
        
    ORFs =  [proteinSeq0, proteinSeq1, proteinSeq2] 
    ORFs = goodORFs(ORFs)
    return ORFs  


# seq='GTTGCAACTGCAGAAGCTGAA' # Wuhan GenomeAnalysis_02102020.py
# seq='GTTGCTACAGCTCACAGCGAG' # SARS

# [proteinSeq0,proteinSeq1,proteinSeq2] = translateDNA3FramesX(seq)
# print(proteinSeq0,proteinSeq1,proteinSeq2)


# get three ORFs starting from snpPos in the genome sequence seq
def getORFs(snpPos, mutatedTo, seq):
    # snpPos = 11083 #THis is actual position in genome order (starts from 1)
    pos = snpPos - 1 #This is the position in Python list in programming
    # mutatedTo = 'T'  #For example, from C->T mutation in SNP
    
    # print('pos0:',seq[pos])
    # print('pos_1:',seq[pos-1])
    # print('pos_2:',seq[pos-2])
    seqT = seq[pos - 2 : pos + 28]
    # print('Before mutation:',seqT)
    
    ORFs = translateDNA3FramesX(seqT)
    print('Amino acid before:',ORFs)
    
    seqM = seqT[0:2]+mutatedTo+seqT[3:len(seqT)]
    # print('After mutation:',seqM)
    mORFs = translateDNA3FramesX(seqM)
    
    print('Amino acid after:',mORFs)
    return [ORFs,mORFs]

#=============================================================================================
# Passed two ORFS, Second one contains mutations, 
# the first one is from reference geneom (no mutation)    
def alignORFs2RefProteom(ORFs, mORFs, refProFasta):
    maxDist = 1
    proteinMutation = {}
    m = len(mORFs)
    print('ORFs:',ORFs)
    print('mORFs:',mORFs)
    print(m)
    
    for i in range(0,m):    # using index mORFs:
        mProtein = mORFs[i]  #muated protein derived from SNP using reference genome
        protein = ORFs[i]    #no mutation in reference genome
        
        for record in SeqIO.parse(refProFasta, "fasta"):
            header = record.id
            refProSeq = str(record.seq)
            
            matches = find_near_matches(mProtein,refProSeq,max_l_dist = maxDist)
            if len(matches)==0:
                # maxDist=3
                # print('what?',matches,mProtein,refProSeq)
                matches = find_near_matches(mProtein,refProSeq,max_l_dist = 2)
                # print('what?',matches,mProtein)
                # print('matches?',matches)
            dists = []
            
            if len(matches)>0:
                # print('2019-nCoV mutation protein:',mProtein)
                # print('SARS-CoV protein name:',header)
                print(matches)
            
                for i in range(0,len(matches)):
                    dists.append(matches[i].dist)
                    seqMatched = refProSeq[matches[i].start-1:matches[i].end+1] 
                    # Move left 1 because python start one less than actual genome position, move right becaue end is one less in python
                    [alignedSeqs,score] = alignLocal(mProtein,seqMatched)
                    print(alignedSeqs)#,'score',score)
                    
                    aminoAcidPos = matches[i].start # in actual protein position
                    print('Original',aminoAcidPos)
                    
                    aminoAcid_original = protein[0]
                    aminoAcid_mutated = mProtein[0]
                    mutationAA = aminoAcid_original+':'+aminoAcid_mutated
                    print('mutated',mutationAA)
                    
                    proteinMutation['protein'] = header
                    proteinMutation['position'] = aminoAcidPos
                    proteinMutation['mutation'] = mutationAA
    return proteinMutation
#============================================================
# https://www.uniprot.org/proteomes/UP000000354
refProFasta = 'SARS-CoV-2_proteome.fasta' #THis is protem from SARS-CoV
# I need to make proteins from 2019-nCoV genome 'NC_045512.fasta' using reference proteins from SARS-CoV!!
# refFasta='NC_045512.fasta'
# idx = 0
# idd,seq = getDNASeqFastaX(refFasta,idx)
# print(idd)

# snpPos = 11083
# This is actual position in genome order (starts from 1)  
# {11083: 'G->T'}
# mutatedTo = 'T' 
# note: 
# pos = 23403
# mutatedTo='G'
# [ORFs,mORFs] = getORFs(snpPos,mutatedTo,seq) 

# proteinMutation =alignORFs2RefProteom(ORFs,mORFs,refProFasta)
# print('YES',snpPos,mutatedTo,proteinMutation)

#=============================================================================
import csv
from csv import DictReader
# virusSNPFile = open('SARS-CoV-2_SNPs.csv', 'w')
# virusDict = {'US|Key|2020-01-02': [1,10], 'US|Key2|2020-01-03':[20,30,100]}
def writeCSV_virusDict(virusDict):
    # csv_columns = ['No','Name']
    virusSNPFile = 'SARS-CoV-2_SNPs.csv'
    
    with open(virusSNPFile, 'w') as csvfile:
        writer = csv.writer(csvfile)
        # writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        # writer.writeheader()
        
        # virusDict = {'US|Key|2020-01-02': [1,10], 'US|Key2|2020-01-03':[20,30,100]}
        for seqId, snp in virusDict.items():
            writer.writerow([seqId, snp])

# virusDict = {'US|Key|2020-01-02': [1,10], 'US|Key2|2020-01-03':[20,30,100]}
# writeCSV_virusDict(virusDict)

# The following is a list of sets, not a list of dictionary
# virusSNPlist=[{241: 'C->T'},{3037: 'C->T'},{8782: 'C->T'},{11083: 'G->T'},{14408: 'C->T'},{23403: 'A->G'},{28144: 'T->C'}]
# virusSNPlist.append({26144,'G->T'})


# Note a list of dictionaries, the dictionary element must have same keys
# Read csv into list of dictionaries using python

def snpRecordsFromCSV(csvSNPName):
    with open(csvSNPName, 'r') as recordsObj:
        dict_reader = DictReader(recordsObj)
        records = list(dict_reader)
        return records  
'''   
csvSNPName ='SARS-CoV-2_SNPs.csv'
records = snpRecordsFromCSV(csvSNPName)
print(records)
'''
#virusSNPs is a list of SNP dictionary
def proteinSNPs(virusSNPs):
    # snp0={'pos':241, 'SNP':'C->T'}
    # snp1={'pos':11083, 'SNP':'G->T'}
    refProFasta ='SARS-CoV-2_proteome.fasta' #THis is protem from SARS-CoV
    
    # virusSNPs=[]
    # virusSNPs.append(snp0)
    # virusSNPs.append(snp1)
    refFasta = 'NC_045512.fasta'
    idx = 0
    idd,refSeq = getDNASeqFastaX(refFasta,idx)
    # print(seq)
    # print(idd)
    
    proteinInfos=[]
    
    for virusDict in virusSNPs:
        # print(virusDict['pos'],virusDict['SNP'])
        snpMutation = {}
        print('SNP records:',virusDict)
        # print('seqId:',virusDict['seqId'])
        snpPos = virusDict['pos']
        mutateInfo = virusDict['SNP']
        mutatedTo = virusDict['SNP'][3]
        snpMutation['pos'] = snpPos
        snpMutation['mutation']= mutateInfo
        
        [ORFs,mORFs] = getORFs(snpPos,mutatedTo,refSeq) 
        print('Initial',ORFs,mORFs)
        proteinMutation = alignORFs2RefProteom(ORFs,mORFs,refProFasta)
        
        # print(snpPos,mutatedTo, proteinMutation)
        # proteinInfo = '{\'SNP\':\''+str(mutateInfo) +'\'},'+ str(proteinMutation)
        proteinInfo = str(snpMutation)+','+str(proteinMutation)
        proteinInfos.append(proteinInfo)
        print(proteinInfo)
    return proteinInfos

# pos is the actual position on genome, need to subtract 1 for python list
virusSNPs = [{'pos': 28144, 'SNP': 'T->A'}]
# proteinInfos = proteinSNPs(virusSNPs)

#virusSNPs is a list of SNP dictionary
def mapGenomeSNPs2Proteom(virusSNPs):
    # snp0={'pos':241, 'SNP':'C->T'}
    # snp1={'pos':11083, 'SNP':'G->T'}
    refProFasta = 'SARS-CoV-2_proteome.fasta' #THis is protem from SARS-CoV
    
    # virusSNPs = []
    # virusSNPs.append(snp0)
    # virusSNPs.append(snp1)
    refFasta = 'NC_045512.fasta'
    idx = 0
    idd,refSeq = getDNASeqFastaX(refFasta,idx)
    # print(seq)
    # print(idd)
    
    proteinInfos=[]
    
    for virusDict in virusSNPs:
        # print(virusDict['pos'],virusDict['SNP'])
        snpMutation = {}
        # print('SNP records:',virusDict)
        # print('seqId:',virusDict['seqId'])
        snpPos = virusDict['pos']
        mutateInfo = virusDict['SNP']
        mutatedTo = virusDict['SNP'][3]
        snpMutation['pos'] = snpPos
        snpMutation['mutation'] = mutateInfo
        
        [ORFs,mORFs] = getORFs(snpPos,mutatedTo,refSeq) 
        # print('Initial',ORFs,mORFs)
        proteinMutation = alignORFs2RefProteom(ORFs,mORFs,refProFasta)
        
        # print(snpPos,mutatedTo, proteinMutation)
        # proteinInfo = '{\'SNP\':\''+str(mutateInfo) +'\'},'+ str(proteinMutation)
        proteinInfo = proteinMutation
        proteinInfos.append(proteinInfo)
        
    return proteinInfos


proteinInfos = mapGenomeSNPs2Proteom(virusSNPs)
print('PPI',proteinInfos)

'''
refProFasta ='SARS-CoV-2_proteome.fasta' #THis is protem from SARS-CoV
# return the translated ORFs after comparing with the reference proteins
ORFs = ['LPFTINCQE', 'YLLQLIARN']
maxDist = 0
for orf in ORFs:
    for record in SeqIO.parse(refProFasta, "fasta"):
        header = record.id
        refProSeq = str(record.seq)
        matches = find_near_matches(orf,refProSeq,max_l_dist = maxDist)
        if len(matches)>0:
            print(matches)
            print(orf,header)
'''    


'''
TEST {'pos': 28144, 'SNP': 'T->A'}
Amino acid before: ['LPFTINCQE', 'YLLQLIARN']
Amino acid after: ['NLLQLIARN']
ORFs: ['LPFTINCQE', 'YLLQLIARN']
mORFs: ['NLLQLIARN']
{'pos': 28144, 'mutation': 'T->A'},{}
'''


#****************************************************************************
'''
virusSNPs = [{'pos':3037,'SNP':'C->T'},{'pos':8782, 'SNP':'C->T'},{'pos':11083, 'SNP':'G->T'},{'pos':14408, 'SNP':'C->T'}]
snp2 = {'pos':23403, 'SNP':'A->G'}
snp0 = {'pos':26144, 'SNP':'G->T'}
snp1 = {'pos':28144, 'SNP':'T->C'}
snp3 = {'pos':18060, 'SNP':'A->G'}#?
snp4 = {'pos':17858, 'SNP':'A->T'}#?
snp5 = {'pos':17747, 'SNP':'T->G'}#?

virusSNPs.append(snp0)
virusSNPs.append(snp1)
virusSNPs.append(snp2)
virusSNPs.append(snp3)
virusSNPs.append(snp4)
virusSNPs.append(snp5)

proteinInfos = proteinMutation(virusSNPs)
SNPFile = 'SNPs.txt' 
with open(SNPFile, 'w') as snpFile:
    for proteinInfo in proteinInfos:
        snpFile.write(proteinInfo+'\n')

seq='AAAAABBBBBCCCCCDDDDDEEEEEFFFF'
win=5
step=1
seqs=getSlidingWindows(seq,win,step)
print(seqs)

refProFasta ='SL-CoV_RaTg13.fasta' #THis is protem from SARS-CoV
proSeqs=[]
for record in SeqIO.parse(refProFasta, "fasta"):
        header = record.id
        refProSeq = str(record.seq)
        proSeqs.append(refProSeq)

spike =proSeqs[2]
print(spike)
'''

def mapSNPs2Proteins(virusSNPs, numNexts):
    # virusSNPs=[{'pos': 241, 'SNP': 'C->T'}, {'pos': 3037, 'SNP': 'C->T'}, {'pos': 23403, 'SNP': 'A->G'}, {'pos': 14408, 'SNP': 'C->T'}, {'pos': 28881, 'SNP': 'G->A'}, {'pos': 28882, 'SNP': 'G->A'}, {'pos': 28883, 'SNP': 'G->C'}, {'pos': 8782, 'SNP': 'C->T'}, {'pos': 28144, 'SNP': 'T->C'}, {'pos': 26144, 'SNP': 'G->T'}, {'pos': 11083, 'SNP': 'G->T'}, {'pos': 17747, 'SNP': 'C->T'}, {'pos': 17858, 'SNP': 'A->G'}, {'pos': 18060, 'SNP': 'C->T'}]
    proteinInfos = proteinSNPs(virusSNPs)
    SNPFile = './snpRecords/Protein_SNPs_%s2020_numNexts_%d.txt'%(SubDate,numN) 
    with open(SNPFile, 'w') as snpFile:
        for proteinInfo in proteinInfos:
            snpFile.write(proteinInfo+'\n')

#-----------------------------------------------------------------------------
'''
msaName GenomesGISAID_SARS-CoV-2_US_MSA.txt
SNPs[8781, 17747, 17858, 18060, 28144] on reference genome: 2019-nCoV|WH01|NC_045512|2020-01-05 {8782: 'C', 17747: 'C', 17858: 'A', 18060: 'C', 28144: 'T'}
SNPs[8781, 17747, 17858, 18060, 28144] on virus genome: US|WA-UW33|EPI_ISL_414620|2020-03-08 {8524: 'T', 17489: 'T', 17600: 'G', 17802: 'T', 27886: 'C'}
'''

'''
proteinFasta ='helicase.fasta'          
#For mutations
refFasta='NC_045512.fasta'

idx = 0
idd,protein = getDNASeqFastaX(proteinFasta,idx)
print(protein)

idd,genome = getDNASeqFastaX(refFasta,idx)
#seq = genome[16236:18039]
#print(seq)

nt=genome[8781]
nt=genome[18060]
print(nt)
'''

'''
ORFs = gene.translateDNA3FramesX(seq)
print('Amino acid before:',ORFs[0])

#{'pos': 17859, 'mutation': 'A->G'}
snpPos = 17859
pos = snpPos-1 #This is the position in Python list in programming

print('nucleotide?',genome[pos:pos+20])
seqO = genome[pos-2:pos+28]
print(seqO)


#----------------------------------------------------------------------------
ORFs = gene.translateDNA3FramesX(seqO)
print('Amino acid sequence original:',ORFs)#[0])
orf0=ORFs[0]
orf1=ORFs[1]

maxDist=1
matches = find_near_matches(orf0,protein,max_l_dist = maxDist)
print('matches',matches)         
seqMatched = protein[matches[0].start:matches[0].end]
print(seqMatched)

[alignedSeqs,score] = alignLocal(orf0,seqMatched)
print(alignedSeqs)#,'score',score)
'''