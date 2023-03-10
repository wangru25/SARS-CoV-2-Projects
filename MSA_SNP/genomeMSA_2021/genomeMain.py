# -*- coding: utf-8 -*-
'''
@Author: Rui Wang
@Date: 2020-04-30 15:45:48
@LastModifiedBy: Rui Wang
@LastEditTime: 2020-05-13 09:43:54
@Email: wangru25@msu.edu
@FilePath: /38_Influ/genomeMSA/genomeMain.py
@Description: Based on Dr. Changchuan Yin's code. Contact: yin1@uic.com
'''
import os
import sys
import csv
import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from os import listdir
from os.path import isfile, join


def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
            print('folder created')
        else:
            print('folder exists ')
    except OSError:
        print ('Error: Creating directory. ' +  directory)
        
class jaccardDis:
    def __init__(self, listA, listB):
        self.listA = listA
        self.listB = listB

    def Jaccard(self):
        unionAB = list(set(self.listA).union(self.listB))
        interAB = list(set(self.listA).intersection(self.listB))
        dist = 1 - len(interAB)/len(unionAB)
        return dist

    def DrJaccard(self):
        dist = self.Jaccard()
        if len(self.listA) > len(self.listB):
            dist = dist*(-1)
        return dist

    def isInSet(self):
        A = set(self.listA)
        B = set(self.listB)
        if len(A) < len(B):
            return A.issubset(B)
        else:
            return False

    def JaccardNextGens(self):
        isNext = False
        dist = self.DrJaccard()
        if self.isInSet():
            isNext = True
        return dist, isNext

    def JaccardRelatives(self):
        isRelative = False
        dist = self.DrJaccard()
        isRelative = True
        return dist, isRelative


class generateSNP:
    def __init__(self, inFile, outFile, indexFile):
        self.inFile = inFile      # inFile = GenomesGISAID_SARS-CoV-2_04162020_MSA_1.txt
        self.outFile = outFile    # outFile = snpRecords_04162020.csv  './snpRecords/snpRecords_%s2020_%d.csv'%(SubDate,count)
        self.align = AlignIO.read(self.inFile, "clustal")
        self.refSeqId = '2019-nCoV|WH01|NC_045512|2020-01-05' #For Clustal Omega: '2019-nCoV|WH01|NC_045512|2020-01-05'
        # self.refSeqId = '2017-Bat|Bat|NC_014470|2008-02-01' #For Clustal Omega: '2019-nCoV|WH01|NC_045512|2020-01-05'
        # self.refSeqId = '2017-Bat|Bat|MG772933B|2017-02-01'
        # self.refSeqId = '2003-nCoV|WH01|NC_045512|2003-01-05'
        self.indexFile = indexFile

    def getRefSeqId(self):
        f = pd.read_csv(self.indexFile)
        fullName = f['Full_Name'].tolist()
        mapName = f['Map_Name'].tolist()
        if self.refSeqId in fullName:
            idxT = fullName.index(self.refSeqId)
            refSeqId = mapName[idxT]
        return refSeqId
        
    def getReferenceMSA(self):
        N = len(self.align)  
        idx = -1
        f = pd.read_csv(self.indexFile)
        fullName = f['Full_Name'].tolist()
        mapName = f['Map_Name'].tolist()
        if self.refSeqId in fullName:
            idxT = fullName.index(self.refSeqId)
            keyword = str(mapName[idxT])
        print('keyword',keyword)
        for record in self.align:
            recordId = (record.id)
            # print('recordId',recordId)
            idx = idx+1
            if keyword == recordId:
                print(recordId, idx)
                break
        refIdx = idx
        # print("?????????",refIdx)
        return [N,refIdx]

    def getMismatch(self, seq, hostSeq):
        cntMatch = 0
        cntTotal = 0
        SNPs = {}
        for i in range(0,len(seq)):
            a = seq[i].upper()
            if a == 'U':
                a = 'T'
            b = hostSeq[i].upper()
            if b == 'U':
                b = 'T'
            # For exact matching only consider non-gaps alignments
            if a != '-' and b != '-':
                if a == b:
                    cntMatch = cntMatch + 1
                cntTotal = cntTotal + 1 #do not count the both are deleted in the MSA file   
            # For SNPs, consider both non-gaps and gaps alignments
            if a != b and (a != '-' and b != '-' and a != 'N' and b != 'N'):
                SNPs[i] = a+':'+b
        dist = cntTotal - cntMatch
        # print(cntTotal,cntMatch,dist)
        # identity = 1- round(cntMatch/cntTotal,4)
        return [dist, SNPs]#,identity]

    def getVirusSNPs(self): 
        [N,refIdx] = self.getReferenceMSA()
        seqRef = self.align[refIdx].seq.upper()
        refVirus = self.align[refIdx].id
        print('RefVirus:', refVirus)
        positionData = []
        virusDict = {}
        for i in range(0, N):
            if i != refIdx:
                seqVirus = self.align[i].seq.upper()
                virusName = self.align[i].id
                # virusNameT = virusName.split('|')
                # virusName = virusNameT[0]+'|'+virusNameT[1]+'|'+virusNameT[3]
                [dist, SNPs] = self.getMismatch(seqRef,seqVirus)
                print('positions', list(SNPs.keys()))
                positions = list(SNPs.keys())  # This is the python position
                positionData.append(positions)
                if len(positions) > 0:
                    virusDict[virusName] = list(SNPs.keys())
        refVirus = ''
        return virusDict

    def mapMSA2RefGenome(self, posMSA):
        refSeq = ''
        idd = ''
        refSeqId = self.getRefSeqId()
        for i in range (0, len(self.align)):
            idd = self.align[i].id
            if (idd == refSeqId):
                refSeq = self.align[i].seq.upper()
                break
        preSeq = refSeq[0: posMSA + 1]  #The sequence before the mismatching position in MSA entries. There are dashed in the sequences (gaps)
        # print(preSeq)
        lastDash = 0
        for ch in preSeq:
            if ch == '-':
                lastDash = lastDash + 1
        ntRefGenome = refSeq[posMSA]     #-lastDash] #Python list start from 0
        posRefGenome = posMSA - lastDash + 1 #actual position in a genome start from 1, python start from 0
        return [posRefGenome, ntRefGenome]


    def mapMSA2VirusGenome(self, seqId, posMSA):
        seqId = seqId.strip()
        seq = ''
        idd = ''   
        for i in range (0, len(self.align)):
            idd = self.align[i].id
            if (idd == seqId):
                seq = self.align[i].seq.upper()
                break         
        preSeq = seq[0 : posMSA + 1]  #The sequence before the mismatching position in MSA entries. There are dashed in the sequences (gaps)
        lastDash = 0
        for ch in preSeq:
            if ch == '-':
                lastDash = lastDash + 1
        posGenome = posMSA - lastDash + 1 #actual position in a genome start from 1, python start from 0
        ntGenome = seq[posMSA] 
        return [posGenome,ntGenome]

    def getGenomePositions(self, seqId, snps):
        genomeRef_PosNTs = {}
        genomeVirus_PosNTs = {}
        genome_mutations = {}
        snpsRef = []
        for posMSA in snps:
            [posGenomeRef,ntRef] = self.mapMSA2RefGenome(posMSA) #only find ref in the MSA file one #YES check it
            snpsRef.append(posGenomeRef)
            # print(seqId,posMSA,'Reference genome position:'+str(posGenomeRef),'NT:'+ntRef)
            genomeRef_PosNTs[posGenomeRef] = ntRef    
            # seqId = strainIds[0] #Assume all the mutation in the same position are the same mutation??? Maynot be true.
            # print('posMSA?',posMSA)
            [posGenome,nt] = self.mapMSA2VirusGenome(seqId,posMSA)
            # print(seqId,posMSA,'Virus genome position:'+str(posGenome),'NT:'+nt)
            if nt == 'U':
                nt = 'T'
            genomeVirus_PosNTs[posGenome]= nt
            ntMutation = ntRef+'->'+nt
            genome_mutations[posGenomeRef] = ntMutation
        return genomeRef_PosNTs, genomeVirus_PosNTs, genome_mutations, snpsRef

    def mapMutations(self, genomeRef_PosNTs, genomeVirus_PosNTs):
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
            if NT_v == 'U':
                NT_v = 'T'
            NTInfo = NT_r+'->'+NT_v
            # print(NTInfo)    
            mutationDict[r] = NTInfo
        # print(mutaionDict)
        return mutationDict
        
    def getAllNextsFromRecords(self, records):
        for recordQ in records:
            snpQ = recordQ['snps']
            cnt = 0
            for recordP in records:
                snp = recordP['snps']
                # print('jaccard snp', snp)
                # print('jaccard snpQ', snpQ)
                self.jaccard = jaccardDis(snp, snpQ)
                dist, isNext = self.jaccard.JaccardNextGens()
                if isNext:
                    cnt = cnt + 1
            recordQ['numNexts'] = cnt
        return records

    def msaSNP2Genome(self):
        snpRecords = []
        virusDict = self.getVirusSNPs()
        refSeqId = self.getRefSeqId()
        f = pd.read_csv(self.indexFile)
        fullName = f['Full_Name'].tolist()
        mapName = f['Map_Name'].tolist()
        for seqId, snps in virusDict.items():
            genomeRef_PosNTs, genomeVirus_PosNTs, genome_mutations, snpsRef = self.getGenomePositions(seqId, snps)
            print('msaName', self.inFile)
            print('SNPs'+str(snps)+' on reference genome:', refSeqId, genomeRef_PosNTs)
            print('SNPs'+str(snps)+' on virus genome:', seqId, genomeVirus_PosNTs)   
            record = {}
            if seqId in mapName:
                idxT = mapName.index(seqId)
                seqId = fullName[idxT]
            record['seqId'] = seqId
            record['snps'] = snps # This is the positions in align entries in the clustal file
            record['snpsRef'] = snpsRef
            record['refNTs'] = genomeRef_PosNTs # This is the reference locations from the aligned positions
            record['mutatedNTs'] = genomeVirus_PosNTs
            mutationMap = self.mapMutations(genomeRef_PosNTs, genomeVirus_PosNTs)
            record['mapNTs'] = mutationMap
            snpRecords.append(record) # A list of record dictionary
        print('Total genomes:',len(virusDict))
        snpRecords = self.getAllNextsFromRecords(snpRecords) #ad numNexts (number of desenstants)
        csv_columns = ['seqId','snps','snpsRef','refNTs','mutatedNTs','mapNTs','numNexts']
        try:
            with open(self.outFile, 'w') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
                writer.writeheader()
                for record in snpRecords:
                    writer.writerow(record)
        except IOError:
            print("I/O error")  

def main(inFile, outFile, indexFile):
    GenerateSNP = generateSNP(inFile,outFile,indexFile)
    GenerateSNP.msaSNP2Genome()

if __name__ == "__main__":
    '''
    SubDates = sys.argv[1].split(',')
    
    oldDate = SubDates[0]
    for idx, SubDate in enumerate(SubDates):
        if idx > 0:
            newDate = SubDates[idx]
            mypath = './clustalW/%s'%newDate
            onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
            createFolder('./snpRecords/%s/'%newDate)
            prefix_max = len(onlyfiles)
            print("processing " + newDate + ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            for i in range(1, prefix_max+1):
                prefix = str(i)
                inFile = './clustalW/%s/GenomesGISAID_SARS-CoV-2_%s2021_MSA_%s.txt'%(newDate,newDate, prefix)
                
                outFile = './snpRecords/%s/snpRecords_%s2021_%s.csv'%(newDate, newDate, prefix)
                if not os.path.exists(outFile):
                    indexFile = '../data/GISAID_2021/%s/%s_mapIndex_%s.csv'%(newDate, prefix, newDate)
                    main(inFile, outFile, indexFile)
            print("Finished processing " + newDate + ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print("Merging " + oldDate + " and " + newDate + ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            arg = ' '.join([oldDate, newDate, str(prefix_max)])
            os.system("python snpRecords/mergeFiles.py %s"%arg)
            print("Finished Merging " + oldDate + " and " + newDate + ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            oldDate=newDate
    '''
    inpath_mast = "clustalW/"
    folder = [x[0] for x in os.walk(inpath_mast)]
    SubDates = []
    outpath_master = "snpRecords/"
    for f in folder:
        SubDates.append(f.split('/')[1])
    
    for date in SubDates:
        outpath = outpath_master + date + "/"
        filename_master = "snpRecords_%s2021_"%date
        if os.path.exists(outpath + filename_master + "new.csv"):
            print(filename_master + "new.csv", "exists")
        else:
            inpath = inpath_mast + date + "/"
            onlyfiles = [f for f in listdir(inpath) if isfile(join(inpath, f))]
            createFolder('./snpRecords/%s/'%date)
            prefix_max = len(onlyfiles)
            print("processing " + date + ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            for i in range(1, prefix_max+1):
                prefix = str(i)
                inFile = './clustalW/%s/GenomesGISAID_SARS-CoV-2_%s2021_MSA_%s.txt'%(date,date, prefix)
                outFile = './snpRecords/%s/snpRecords_%s2021_%s.csv'%(date, date, prefix)
                if not os.path.exists(outFile):
                    indexFile = '../data/GISAID_2021/%s/%s_mapIndex_%s.csv'%(date, prefix, date)
                    main(inFile, outFile, indexFile)
            print("Finished processing " + date + ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    '''
    SubDate = "1210"
    prefix_min = 45
    prefix_max = 54
    
    for i in range(prefix_min, prefix_max + 1):
        
        prefix = str(i)
        inFile = './clustalW/GenomesGISAID_SARS-CoV-2_%s2020_MSA_%s.txt'%(SubDate, prefix)
        outFile = './snpRecords/snpRecords_%s2020_%s.csv'%(SubDate, prefix)
        indexFile = '../data/GISAID/%s/%s_mapIndex_%s.csv'%(SubDate, prefix, SubDate)
        main(inFile, outFile, indexFile)
    '''
    
    '''
    SubDate = "1019"
    prefix = '22'
    inFile = './clustalW/GenomesGISAID_SARS-CoV-2_%s2020_MSA_%s.txt'%(SubDate, prefix)
    outFile = './snpRecords/snpRecords_%s2020_%s.csv'%(SubDate, prefix)
    indexFile = '../data/GISAID/%s/%s_mapIndex_%s.csv'%(SubDate, prefix, SubDate)
    main(inFile, outFile, indexFile)
    '''

# def reverseComplement(seq):
#     seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
#     return "".join([seq_dict[base] for base in reversed(seq)])

# # print(reverseComplement('ATATTGCAGCAGTACGCACACA'))

# from Bio import SeqIO
# from Bio import Entrez
# from Bio import pairwise2
# from Bio.pairwise2 import format_alignment
# def getDNASequencesFasta(fastaFileName,idd):
#     '''
#     fastaFileName='./Protein/SARS-CoV-2_proteome.fasta'
#     idd = 'nsp3'
#     proteinSeq = getDNASequencesFasta(fastaFileName,idd)
#     print(proteinSeq[0])
#     '''
#     from Bio.Seq import Seq
#     from Bio import SeqIO
#     seqR =''
#     for record in SeqIO.parse(fastaFileName, "fasta"):
#         #print(record.id) #header
#         idx = record.id
#         seq = record.seq
#         if idx == idd:
#             seqR = seq
#     return seqR

#     # Usage example:
#     '''
#     fastaFileName='SARS-CoV-2_proteome.fasta'
#     idd = 'YP_009725300.1|nsp4|8555:10054'
#     proteinSeq = getDNASequencesFasta(fastaFileName,idd)

#     actualPos= 75
#     print(proteinSeq[actualPos-1])
#     '''

# # fastaFileName = 'PMRP.fasta'
# # idd = 'NM_001104546.2'
# fastaFileName = 'NC_045512.fasta'
# idd = '2019-nCoV|WH01|NC_045512|2020-01-05'
# proteinSeq = getDNASequencesFasta(fastaFileName,idd)
# # seq = 'TTGCTGCTGCTTGACAGATT'
# seq = reverseComplement('GACCCCAAAATCAGCGAAAT')
# print(seq)
# nPos = proteinSeq.find(seq)
# print(nPos)
# a = '%d:%d'%(nPos+1, nPos+len(seq))
# print('actual postion:', a)
