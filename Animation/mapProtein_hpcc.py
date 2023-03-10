# -*- coding: utf-8 -*-
'''
@Author: Rui Wang
@Date: 2020-05-03 08:30:58
@LastModifiedBy: Rui Wang
LastEditTime: 2021-01-25 16:24:42
@Email: wangru25@msu.edu
FilePath: /38_Influ/analysis/Animation/mapProtein.py
@Description: Run mapProtein.py 01012020 first, and then proteinMain.py, and then animationFreq_Pro.py
'''
import pandas as pd 
import numpy as np
import os
import sys
import ast
from Bio import AlignIO
from Bio import SeqIO
from Bio import Entrez
from Bio import pairwise2
from fuzzysearch import find_near_matches

# Convert long records to virus dictionary
def getVirusDictFromRecords(records):
    virusDict = {}
    for i in range(records.shape[0]):
        seqId = records['seqId'][i]
        snp = records['snpsRef'][i]
        virusDict[seqId]=snp
    return virusDict

def createSNP(mapNTs):
    snps = []
    for pos,mutation in mapNTs.items():
        snp = {}
        snp['pos'] = pos
        snp['SNP'] = mutation
        snps.append(snp)
    print(snps)
    return snps

def uniqueSNPfromRecords(records,numN):
    # records = pd.read_csv('snpRecords_04052020.csv')
    snpss = []
    for i in range(records.shape[0]):
        if records['numNexts'][i] >= numN:       
            # print('High',records.loc[[i]])
            mapNTs = records['mapNTs'][i]
            mapNTs = ast.literal_eval(mapNTs)
            snps = createSNP(mapNTs)
            snpss = snpss+snps
    # print('snpss',snpss)
    virusSNPs = [i for n,i in enumerate(snpss) if i not in snpss[n+1:]]
    # virusSNPs = list({v['pos']:v for v in snpss}.values()) #remove duplicated dictionaries.   
    return  virusSNPs    


def mapSNPs2Proteins(virusSNPs, numNexts, Date, OldDate,nums):
    # virusSNPs=[{'pos': 241, 'SNP': 'C->T'}, {'pos': 3037, 'SNP': 'C->T'}, {'pos': 23403, 'SNP': 'A->G'}, {'pos': 14408, 'SNP': 'C->T'}, {'pos': 28881, 'SNP': 'G->A'}, {'pos': 28882, 'SNP': 'G->A'}, {'pos': 28883, 'SNP': 'G->C'}, {'pos': 8782, 'SNP': 'C->T'}, {'pos': 28144, 'SNP': 'T->C'}, {'pos': 26144, 'SNP': 'G->T'}, {'pos': 11083, 'SNP': 'G->T'}, {'pos': 17747, 'SNP': 'C->T'}, {'pos': 17858, 'SNP': 'A->G'}, {'pos': 18060, 'SNP': 'C->T'}]
    proteinInfos = proteinSNPs(virusSNPs)
    SNPFile = '/mnt/home/wangru25/Rui_Wang/1_Training/38_Influ/analysis/proteinPlot/Protein/Protein_SNPs_%s_new_%d.txt'%(Date,nums)
    with open(SNPFile, 'w') as snpFile:
        for proteinInfo in proteinInfos:
            snpFile.write(proteinInfo+'\n')
    os.chdir('/mnt/home/wangru25/Rui_Wang/1_Training/38_Influ/analysis/proteinPlot/Protein/')
    os.system('cat Protein_SNPs_%s.txt Protein_SNPs_%s_new_%d.txt | sort |uniq > Protein_SNPs_%s.txt'%(OldDate, Date, nums,Date))

def proteinSNPs(virusSNPs):
    # snp0={'pos':241, 'SNP':'C->T'}
    # snp1={'pos':11083, 'SNP':'G->T'}
    refProFasta ='SARS-CoV-2_proteome.fasta' #This is protem from SARS-CoV-2
    # print(virusSNPs)
    
    # virusSNPs=[]
    # virusSNPs.append(snp0)
    # virusSNPs.append(snp1)
    refFasta = 'NC_045512.fasta'
    idx = 0
    idd,refSeq = getDNASeqFastaX(refFasta,idx)
    # print('seq',refSeq)
    # print('idd',idd)
    
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
        # print('snpPos', snpPos)
        # print('mutateInfo', mutateInfo)
        # print('mutatedTo', mutatedTo)

        [ORFs,mORFs] = getORFs(snpPos,mutatedTo,refSeq) 
        print('Initial',ORFs,mORFs)
        proteinMutation = alignORFs2RefProteom(ORFs,mORFs,refProFasta)
        
        # print(snpPos,mutatedTo, proteinMutation)
        # proteinInfo = '{\'SNP\':\''+str(mutateInfo) +'\'},'+ str(proteinMutation)
        proteinInfo = str(snpMutation)+','+str(proteinMutation)
        proteinInfos.append(proteinInfo)
        print(proteinInfo)
    return proteinInfos

def getDNASeqFastaX(seqFasta,idx):
    seqs = []
    ids = []
    for record in SeqIO.parse(seqFasta, "fasta"):
        print(record.id)  
        idd = record.id
        idd = idd.replace('|','_')
        ids.append(idd)
        seq = record.seq
        seqs.append(seq)
    seq = removeNonATCG(seq) 
    seq = seqs[idx]
    idd = ids[idx]
    return idd,seq

def removeNonATCG(seq):
    seq = str(seq)
    for ch in seq:
        if ch not in ['A','T','C','G']:
            seq = seq.replace(ch,'')
    return seq

def getORFs(snpPos, mutatedTo, seq):
    # Four parameters (2,-1,-1,-1) in alignment:Identical characters are given 2 points, 1 point is deducted for each non-identical character.
    # 1 point is deducted when opening a gap, and 1 point is deducted when extending it.
    
    # snpPos = 11083 #THis is actual position in genome order (starts from 1)
    pos = snpPos - 1 #This is the position in Python list in programming
    # mutatedTo = 'T'  #For example, from C->T mutation in SNP
    
    print('pos_0:',seq[pos])
    print('pos_1:',seq[pos-1])
    print('pos_2:',seq[pos-2])
    seqT = seq[pos - 2 : pos + 28]
    print('Before mutation:',seqT)
    
    ORFs = translateDNA3FramesX(seqT)
    print('Amino acid before:',ORFs)
    
    seqM = seqT[0:2]+mutatedTo+seqT[3:len(seqT)]
    print('After mutation:',seqM)
    mORFs = translateDNA3FramesX(seqM)
    
    print('Amino acid after:', mORFs)
    return [ORFs,mORFs]

def alignLocal(seqQuery,seqFound):
    # print('Seq1',seqQuery)
    # print('Seq2',seqFound)
    seqQuery = str(seqQuery) #It can take Biopython Seq object
    seqFound = str(seqFound)
    alignments = pairwise2.align.localms(seqQuery,seqFound, 2, -1, -1, -1)

    print("??????????????", alignments)
    if len(alignments) == 0:
        alignments = ''
        score = 0
    else:
        match = []
        score=0
        print('??????????????????',alignments[0][0])
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

        print('???????======================', [alignedSeqs,score])
        
        return [alignedSeqs,score]

# Removed these ORFs that have "stop codon"
def goodORFs(ORFs):
    #ORFs=['FYENAFLPFA', 'FMKMPFYLL', 'L_KCLFTFC']
    goodORFs = []
    for orf in ORFs:
        goodORFs.append(orf)

    # for orf in ORFs:
    #     if '_' not in orf:
    #         goodORFs.append(orf)
    return goodORFs

def isValid(codon):
    #ORFs=['FYENAFLPFA', 'FMKMPFYLL', 'L_KCLFTFC']
    notPossibleNeu = 'BDEFHIJKLMNOPQRSUVWXYZabcdefghijklmnopqrstuvwxyz_'
    saveIdx =[]
    for a in codon:
        if a not in notPossibleNeu:
            idx = 1
            saveIdx.append(idx)
        else:
            idx = 0
            saveIdx.append(idx)
    if len(saveIdx) != np.sum(saveIdx):
        codon = 'TAA'
    else:
        codon = codon
    return codon

def geneticCodon(codon):
    aa =''
    isCodon = True
    for ch in codon:
       if ch not in ['A','T','C','G']:
          isCodon = False
          break
    if isCodon:
        aa = codontable[codon]
    else:
        aa = '#'
    return aa

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
        print('codon0s', codon0s)
        codon1 = seq[i+1:i+4]
        codon1s = codon1s + codon1+ '|'
        print('codon1s', codon1s)
        codon2 = seq[i+2:i+5]
        codon2s = codon2s + codon2+ '|'
        print('codon2s', codon2s)
        
        if len(codon0) == 3:
            # proteinSeq0 += codontable[isValid(codon0)]
            proteinSeq0 += geneticCodon(codon0) #codontable[codon0]
            
        if len(codon1) == 3:
            # proteinSeq1 += codontable[isValid(codon1)]
            proteinSeq1 += geneticCodon(codon1) #codontable[codon1]
        
        if len(codon2) == 3:
            # proteinSeq2 += codontable[isValid(codon2)]
            proteinSeq2 += geneticCodon(codon2) #codontable[codon2]
        
    ORFs =  [proteinSeq0, proteinSeq1, proteinSeq2] 
    ORFs = goodORFs(ORFs)
    return ORFs

def alignORFs2RefProteom(ORFs, mORFs, refProFasta,group=False):
    maxDist = 1
    proteinMutation = {}
    m = len(mORFs)
    # print('ORFs:',ORFs)
    # print('mORFs:',mORFs)
    print(m)
    
    if len(ORFs) != len(mORFs):
        pass
    else:
        for i in range(0,m):    # using index mORFs:
            if mORFs[i] == '':
                proteinMutation = proteinMutation
            else:
                mProtein = mORFs[i]  #muated protein derived from SNP using reference genome
                protein = ORFs[i]    #no mutation in reference genome
                print('mProtein===================================================', mProtein)
                print('protein****************************************************', protein)
                
                for record in SeqIO.parse(refProFasta, "fasta"):
                    header = record.id
                    refProSeq = str(record.seq)
                    # print('=====================================%s==================================\n'%header,refProSeq)
                    # print('mProtein \n',mProtein)
                    # print('refProSeq \n',refProSeq)
                    
                    matches = find_near_matches(mProtein,refProSeq, max_l_dist = maxDist)
                    # print('matches',matches)
                    # print('len', len(matches))
                    if len(matches)==0:
                        # maxDist=2
                        # print('what?',matches,mProtein,refProSeq)
                        matches = find_near_matches(mProtein,refProSeq, max_l_dist = 1)
                        print('what?',matches,mProtein)
                        # print('matches2?',matches)
                    dists = []
                    
                    if len(matches)>0:
                        # print('2019-nCoV mutation protein:',mProtein)
                        # print('SARS-CoV protein name:',header)
                        # print('matches',matches)
                        # print(len(matches))
                    
                        for i in range(0,len(matches)):
                            dists.append(matches[i].dist)
                            # seqMatched = refProSeq[matches[i].start-1 : matches[i].end + 1] 
                            seqMatched = refProSeq[matches[i].start : matches[i].end] 
                            print('start',matches[i].start)
                            print('end',matches[i].end)

                            print('seqMatched',seqMatched)
                            # Move left 1 because python start one less than actual genome position, move right becaue end is one less in python
                            # [alignedSeqs,score] = alignLocal(mProtein,seqMatched)
                            # print('alignedSeqs',alignedSeqs)#,'score',score)
                            
                            aminoAcidPos = matches[i].start # in actual protein position
                            print('Original', aminoAcidPos)
                            
                            aminoAcid_original = protein[0]
                            aminoAcid_mutated = mProtein[0]
                            mutationAA = aminoAcid_original+':'+aminoAcid_mutated
                            print('mutated',mutationAA)

                            if aminoAcid_mutated == aminoAcid_original:
                                realPosition = aminoAcidPos + 1
                            else: 
                                realPosition = aminoAcidPos
                            
                            proteinMutation['protein'] = header
                            proteinMutation['position'] = realPosition
                            proteinMutation['mutation'] = mutationAA
    return proteinMutation

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

mapProteins = True
if mapProteins:
    Date = sys.argv[1] # 09252020
    OldDate = sys.argv[2] #01202021
    nums = int(sys.argv[3]) #int 0, 1, 2
    # virusSNPs={8782: 'C->T', 17747: 'C->T', 17858: 'A->G', 18060: 'C->T', 28144: 'T->C'}
    # virusSNPs = list({v['pos']:v for v in snpss}.values()) #remove duplicated dictionaries 
    # records = pickle.load(open("./snpRecords/snpRecords_%s2020_%d.pkl"%(SubDate, count), "rb" ))
    records = pd.read_csv('snpRecords_%s_new_%d.csv'%(Date,nums))
    virusDict = getVirusDictFromRecords(records)
    print(virusDict)
    virusSNPs = uniqueSNPfromRecords(records,0)
    print(virusSNPs)
    mapSNPs2Proteins(virusSNPs, 0, Date, OldDate, nums) #produce 'Protein_SNPs_0417.txt'
    print('Completed protein mapping')  
