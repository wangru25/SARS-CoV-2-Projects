# -*- coding: utf-8 -*- 
""" 
* @Author: Changchuan Yin
* @Date: 2020-04-06 17:14:26 
* @Last Modified by:   Rui Wang 
* @Last Modified time: 2020-04-06 17:14:26  
* @Contact: cyin1@uic.edu
* @Contact: wangru25@msu.edu 
"""

# DO NOT remove the usage examples
import os
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
# https://pypi.python.org/pypi/fuzzysearch/0.3.0#downloads
from fuzzysearch import find_near_matches
import genomeGraph as gg
import genomeUtility as gh

SubDate = sys.argv[1]   # 0405: Represent for submission date: April 05
count = int(sys.argv[2])
numN = int(sys.argv[3])  # 10: Extract sequences with at least 10 descendants. A nonnegative integer.

#==============================================================================
# https://en.wikipedia.org/wiki/Jaccard_index
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
        seq = record.seq
        seq = str(seq) 
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
    # handle = open(fileName, 'w')
    handle.write('>'+header+'\n')
    width = 70
    for i in range(0, len(seq), width):
        seqT = seq[i:i+width]
        handle.write(seqT + '\n')

#==============================================================================
# To write a fasta file of a list of headers and sequences:
# Inputs: list headers, list sequences, and file name
# OutputS: fasta file containing the headers and sequences
#==============================================================================
def writeFastaFile(headers,seqs,handle):
    # handle = open(fileName, 'w')
    for j in range(0,len(seqs)):
        header = headers[j]
        handle.write('>'+header+'\n')
        seq = seqs[j]
        width = 70
        for i in range(0, len(seq), width):
            seqT = seq[i:i+width]
            handle.write(seqT + '\n')

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
    # seq=Seq(seq)
    complementSeq = seq.reverse_complement()
    # complementSeq = str(complementSeq)
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
    n = len(seq)
    lastPos = n - win + step
    seqs=[]
    for i in range(0,lastPos,step):
        startPos = i
        endPos = i + win
        seqW = seq[startPos:endPos]
        seqs.append(seqW)
        # print('Seq:',seqW)
    # print('Last seq',seqW)
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
# Search forward subsequence in seqT that matches query sequence seqQ
def searchSequence(seqQ, seqT, maxDist):
    seqQ = seqQ.upper()
    seqT = seqT.upper()
    
    isAlmost = False
    seqFound = ''
    distList = []
    startList = []
    endList = []
    minDist = 100
    start = 0
    end = 0
    
    matches = find_near_matches(seqQ, seqT,max_l_dist=maxDist)
    # print('Matches',matches)
    
    numMatches = len(matches)
    if numMatches > 0:
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

# --------------------------------------------------------------------------------
def getDNASeqFastaX(seqFasta, k):
    seqs = []
    ids = []
    for record in SeqIO.parse(seqFasta, "fasta"):
        # print(record.id)  
        idx = record.id
        idx = idx.replace('|','_')
        ids.append(idx)
        # print(record.seq)
        seq = record.seq
        # print('Len',len(seq))
        seqs.append(seq)
    seq = removeNonATCG(seq) 
    seq = seqs[k]
    _id = ids[k]
    # print(_id)
    return _id,seq


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
    handle = Entrez.efetch(db = 'nuccore',  rettype = "gb", id = genBankId, retmode = 'text')
    for seqRecord in SeqIO.parse(handle, 'genbank'):
        genBankId = seqRecord.id
        desc = seqRecord.description
        seq = str(seqRecord.seq)
    return [desc,seq]  
    # rettype = "gb"
    # usage example

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
        cds = seq[int(startPos)-1:int(endPos)]
        
        if strand == -1:#=='complement':
            cds = cds.reverse_complement()

        aa = cds.translate()
        print('proteinSeq translated from genBankId:\n',aa)
    cds = str(cds)
    return cds 
    

#==============================================================================
# Function retrieves gene and protein sequence for a given gene name
# Inputs: a genome accessID, and gene name
# OutputS: proteinId, gene sequence and protein sequence of the gene name
#==============================================================================
def getGeneSequence(genBankId,gene):
    # gene = 'vanRM'
    # genBankId = 'FJ349556'
    startPos = 0
    endPos = 0
    strand = 0
    geneSeq = ''
    AASeq = ''
    # AASeq2 = ''
    proteinId = ''
    # printeinSeq =''
    # key ='translation'
    keys= {'gene', 'translation'}# <= set(some_dict)
    Entrez.email = 'cyin1@uic.edu'
    handle = None

    handle = Entrez.efetch(db = 'nuccore', rettype = 'gbwithparts', id = genBankId, retmode = 'text')
    for seqRecord in SeqIO.parse(handle, 'genbank'):
        genBankId=seqRecord.id
        # desc = seqRecord.description
        seq = seqRecord.seq
        
        # print('seqRecord:\n',seqRecord)
        # print('Id:', genBankId)
        # print('desc:',desc) 
        # print('Seq',seq)
        
        for seqFeature in seqRecord.features :
            #print('sequence feature:\n',seqFeature)
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
                    geneSeq = seq[startPos:endPos]
                    # AASeq2=Seq.translate(geneSeq)
                    # print('AA_TranslatedGene',AASeq2)
                else:
                    geneSeq = seq[startPos:endPos]
                    # print('geneSeq:\n',geneSeq)
                    geneSeq = geneSeq.reverse_complement()
                    # tb = seqFeature.qualifiers['transl_table']
                    # AASeq2 = Seq.translate(geneSeq)
                    # print('AA_TranslatedGene',AASeq2)
                break
            
    return [proteinId,str(geneSeq),str(AASeq)]    

# Usage example 1 
# https://www.ncbi.nlm.nih.gov/nuccore/FJ349556
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
genomeGenBankId = 'KJ660346'
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
    s.append(alignments[0][0]+ '\n')
    s.append(m + '\n')
    s.append(alignments[0][1])
    # s = str(s)
    alignedSeqs = "".join(s)
    
    return [alignedSeqs,score]

#=====================================================TEST====================================================================
# Align two DNA sequences in the starting positions first
def alignDNASequences(virusSeq,hostSeq,win=100,maxDist = 10):
    
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

# virusSeqFile2 = 'EPI_ISL_410721.fasta'
# idHost,hostSeq = getDNASeqFastaX(virusSeqFile2,0)

# lenS = len(virusSeq)
# data is obtained by maxDist = 6 and win = 21
# i.e., I found these sequence segments of length 21 in 2019-nCoV, and these segments have at least 6 mismatches than the 
# RaTG13 strain (corresponding locations), so these segments are very different from these segments in RaTG13, these
# segments in 2019-nCoV are very mutated from RaTG13 or these segments may come from other species (recombination).


def getMismatch(seq,hostSeq):
    cntMatch = 0
    cntTotal = 0
    SNPs={}
    for i in range(0,len(seq)):
        a = seq[i]
        b = hostSeq[i]
        # For exact matching only consider non-gaps alignments
        if a != '-' and b != '-':
            if a == b:
                cntMatch = cntMatch+1
            cntTotal = cntTotal+1 #do not count the both are deleted in the MSA file   
    
        # For SNPs, consider both non-gaps and gaps alignments
        if a != b and (a != '-' and b != '-' and a != 'N' and b != 'N'):
            SNPs[i] = a+':'+b
            
    dist = cntTotal-cntMatch
    # print(cntTotal,cntMatch,dist)
    # identity = 1- round(cntMatch/cntTotal,4)
    return [dist, SNPs]#,identity]

def getGaps(seq,hostSeq):

    gapsDict={}
    for i in range(0,len(seq)):
        a = seq[i]
        b = hostSeq[i]
        # only check gaps not in the two ends alignments
        if (i > 250 and i < 28000):
            # For SNPs, consider both non-gaps and gaps alignments
            if (a == '-' or b == '-'):
                gapsDict[i]=a + ':' + b

    return gapsDict
'''
[dist,variants] = getMismatch(seqY,seqX)
print(variants)

keyX = list(variants.keys())
print(keyX)


gaps = getGaps(seqY,seqX)
print('Gaps',gaps)
'''

# Totally based on MSA results
def get21M6VariantsClustalW(virusSeq, hostSeq, win = 21):
    variants = []
    cntInDel = 0
    lenS = len(virusSeq)
    e = lenS-win+1
    
    for i in range(win,e):
        t = i + win
        seqQ = virusSeq[i:t]
        hostSeqT = hostSeq[i:t]     
            
        [dist, SNPs] = getMismatch(seqQ,hostSeqT)
        # print(identity)
        
        if dist >= 6:
            posVirus = i - cntInDel
            # print(posVirus)
            variants.append(posVirus)
        
        if virusSeq[i] == '-':
            cntInDel = cntInDel + 1
            
    return variants

# colorsStrains,labelStrains  = getStrainColors(strainNames)
# Examples headers of the genome fasta records
strainNames = ['USA|CA1|ELS1234|2020-02-12','USA|WA1|ELS1234|2020-0223','USA|IL2|ELS1234|2020-0223','USA|CA2|ELS1234|2020-0223','USA|IL1|ELS1234|2020-0223']
def getCityName(strainName):
    # strainName = 'USA|CA1|ELS1234|2020-0223'
    strainNameS = strainName.split('|')
    # print(strainName)
    
    labelName = strainNameS[0]+ '|'+ strainNameS[1]#strainNameS[3]
    # print(labelName)
    
    cityName = strainNameS[1][0:2] # The first two characters of cities are the same
    # print(cityName)
    return cityName, labelName

# Color a strain sequence in the same color if the strain is in the same city.
def getStrainColors(strainNames):
    cityNames = []
    colorsCity = {}
    cnt = 0
    for strainName in strainNames:
        cityName, labelName = getCityName(strainName)
        # print(cityName, labelName)
        
        if cityName not in cityNames:
            cnt = cnt + 1
            cityNames.append(cityName)
            c = 'C{}'.format(cnt%10) #only 10 colors are supported
            colorsCity[cityName] = c

    # label each strain by a color
    colorsStrains = []
    labelStrains =[]
    for strainName in strainNames:
        cityName, labelName = getCityName(strainName)
        colorsStrains.append(colorsCity[cityName])
        labelStrains.append(labelName)
    # print(colorsStrains)
    return colorsStrains,labelStrains

# colorsStrains,labelStrains  = getStrainColors(strainNames)
# print(colorsStrains)
# print(labelStrains)

from Bio import AlignIO
#-----------------------------------------------------------------------
def plotSNPs(msaName):
    align = AlignIO.read(msaName, "clustal")
    
    [N,refIdx] = gh.getReferenceMSA(msaName)
    
    seqRef = align[refIdx].seq
    refVirus = align[refIdx].id
    print(refVirus)
    
    positionData = []
    colorLabels = []
    virusDict = {}
    strainNamesX = []
    
    for i in range(0,N):
        if i != refIdx:
            seqVirus = align[i].seq
            virusName = align[i].id
            virusNameT = virusName.split('|')
            virusName = virusNameT[0]+'|'+virusNameT[1]+'|'+virusNameT[3]
            [dist,SNPs] = getMismatch(seqRef,seqVirus)
            # colorLabels.append(align[i].id)
            colorLabels.append(virusName)
            print('positions',list(SNPs.keys()))
        
            positions = list(SNPs.keys())
            positionData.append(positions)
            strainNamesX.append(virusName)
            
            if len(positions)>0:
                virusDict[virusName] = list(SNPs.keys())
            
    # strainNames = list(viruses.keys())
    colorsStrains, labelStrains = getStrainColors(strainNamesX)
    refVirus = ''
    gg.eventPlot2(refVirus, colorLabels, positionData, colorsStrains,lineWidth=4)
    
    return virusDict

def plotSNPs_subset(msaName,stateName):
    align = AlignIO.read(msaName, "clustal")
    
    [N,refIdx] = gh.getReferenceMSA(msaName)
    # N is the number of sequences in the MSA file
    
    seqRef = align[refIdx].seq
    refVirus = align[refIdx].id
    print(refVirus)
    
    positionData=[]
    colorLabels=[]
    virusDict ={}
    strainNamesX =[]

    for i in range(0,N):
        if i != refIdx:
            seqVirus = align[i].seq
            virusName = align[i].id
            
            # if 'CruiseA' in virusName: get those having CruiseA
            if 'CruiseA' not in virusName:
                virusNameT = virusName.split('|')
                virusName = virusNameT[0]+'|'+virusNameT[1]+'|'+virusNameT[3]
                [dist,SNPs] = getMismatch(seqRef,seqVirus)
                # colorLabels.append(align[i].id)
                colorLabels.append(virusName)
                print('positions',list(SNPs.keys()))
            
                positions = list(SNPs.keys())
                positionData.append(positions)
                strainNamesX.append(virusName)
            
            if len(positions) > 0:
                virusDict[virusName] = list(SNPs.keys())
    
    # strainNames = list(viruses.keys())
    colorsStrains,labelStrains = getStrainColors(strainNamesX)
    refVirus = ''
    gg.eventPlot2(refVirus, colorLabels, positionData, colorsStrains, lineWidth = 4)
    
    return virusDict

def plotSNPs_total(msaName):
    align = AlignIO.read(msaName, "clustal")
    
    [N,refIdx] = gh.getReferenceMSA(msaName)
    
    seqRef = align[refIdx].seq
    refVirus = align[refIdx].id
    print(refVirus)
    
    positionData = []
    colorLabels = []
    virusDict = {}
    strainNamesX = []
    
    for i in range(0,N):
        if i != refIdx:
            seqVirus = align[i].seq
            virusName = align[i].id
            virusNameT = virusName.split('|')
            virusName = virusNameT[0]+'|'+virusNameT[1]+'|'+virusNameT[3]
            [dist,SNPs] = getMismatch(seqRef,seqVirus)
            # colorLabels.append(align[i].id)
            colorLabels.append(virusName)
            print('positions',list(SNPs.keys()))
        
            positions = list(SNPs.keys())
            positionData.append(positions)
            strainNamesX.append(virusName)
            
            if len(positions) > 0:
                virusDict[virusName] = list(SNPs.keys())
    
    # strainNames = list(viruses.keys())
    colorsStrains,labelStrains  = getStrainColors(strainNamesX)
    refVirus = ''
    gg.eventPlotT(refVirus,colorLabels,positionData, colorsStrains,lineWidth=4)
    
    return virusDict


def plotSNPs_noID(virusDict):
    positionData = []
    colorLabels = []
    # virusDict = {}
    strainNamesX =[]
    
    for seqId,snp in virusDict.items():
        # if i != refIdx:
            # seqVirus = align[i].seq
            virusName = seqId#align[i].id
            virusNameT = virusName.split('|')
            virusName = virusNameT[0]+'|'+virusNameT[1]+'|'+virusNameT[3]
            # [dist,SNPs] = getMismatch(seqRef,seqVirus)
            # colorLabels.append(align[i].id)
            colorLabels.append(virusName)
        # print('positions',list(SNPs.keys()))
        
            positions = snp#list(SNPs.keys())
            positionData.append(positions)
            strainNamesX.append(virusName)
            
            # if len(positions)>0:
            #     virusDict[virusName] = list(SNPs.keys())
    
    # strainNames = list(viruses.keys())
    colorsStrains,labelStrains  = getStrainColors(strainNamesX)
    refVirus = ''
    gg.eventPlotT(refVirus,colorLabels,positionData, colorsStrains,lineWidth=4)
    
    # return virusDict

# Plot with country code
def plotSNPs_Y(virusDict,yLabel):

    positionData = []
    colorLabels = []
    # virusDict = {}
    strainNamesX = []
    
    for seqId, snp in virusDict.items():
        # if i != refIdx:
            #seqVirus = align[i].seq
            virusName = seqId#align[i].id
            virusNameT = virusName.split('|')
            virusName = virusNameT[0]+'|'+virusNameT[1]+'|'+virusNameT[3]
            # [dist,SNPs] = getMismatch(seqRef,seqVirus)
            # colorLabels.append(align[i].id)
            colorLabels.append(virusName)
        # print('positions',list(SNPs.keys()))
        
            positions = snp#list(SNPs.keys())
            positionData.append(positions)
            strainNamesX.append(virusName)
            
            # if len(positions)>0:
            #     virusDict[virusName]=list(SNPs.keys())
    
    # strainNames = list(viruses.keys())
    colorsStrains,labelStrains  = getStrainColors(strainNamesX)
    refVirus=''
    gg.eventPlotY(yLabel,refVirus,colorLabels,positionData, colorsStrains,lineWidth=4)


def plotSNPs_ID(virusDict):
    # align = AlignIO.read(msaName, "clustal")
    
    # [N,refIdx] = gm.getReferenceMSA(msaName)
    
    # seqRef = align[refIdx].seq
    # refVirus = align[refIdx].id
    # print(refVirus)
    
    positionData = []
    colorLabels = []
    # virusDict = {}
    strainNamesX = []
    
    for seqId, snp in virusDict.items():
        #if i != refIdx:
            #seqVirus = align[i].seq
            virusName = seqId#align[i].id
            virusNameT = virusName.split('|')
            virusName = virusNameT[0]+'|'+virusNameT[1]+'|'+virusNameT[3]
            # [dist,SNPs] = getMismatch(seqRef,seqVirus)
            # colorLabels.append(align[i].id)
            colorLabels.append(virusName)
            # print('positions',list(SNPs.keys()))
        
            positions = snp#list(SNPs.keys())
            positionData.append(positions)
            strainNamesX.append(virusName)
            
            # if len(positions)>0:
                # virusDict[virusName]=list(SNPs.keys())
    
    # strainNames = list(viruses.keys())
    colorsStrains,labelStrains  = getStrainColors(strainNamesX)
    refVirus = ''
    gg.eventPlot2(refVirus,colorLabels,positionData, colorsStrains,lineWidth=4)
    
    #return virusDict

msaName = './clustalW/GenomesGISAID_SARS-CoV-2_%s2020_MSA_%d.txt'%(SubDate, count)
# viruses = plotSNPs(msaName)

#====================================================
def getVirusSNPs(msaName): 

    align = AlignIO.read(msaName, "clustal")
    [N,refIdx] = gh.getReferenceMSA(msaName)
    
    seqRef = align[refIdx].seq
    refVirus = align[refIdx].id

    print('RefVirus:', refVirus)
    positionData = []
    colorLabels= []
    virusDict = {}
    strainNamesX = []
    
    for i in range(0,N):
        if i != refIdx:
            seqVirus = align[i].seq
            virusName = align[i].id
            # virusNameT = virusName.split('|')
            # virusName = virusNameT[0]+'|'+virusNameT[1]+'|'+virusNameT[3]
            [dist,SNPs] = getMismatch(seqRef,seqVirus)
            # colorLabels.append(align[i].id)
            colorLabels.append(virusName)
            print('positions',list(SNPs.keys()))
        
            positions = list(SNPs.keys())
            positionData.append(positions)
            strainNamesX.append(virusName)
            
            if len(positions)>0:
                virusDict[virusName] = list(SNPs.keys())
    
    # strainNames = list(viruses.keys())
    colorsStrains,labelStrains  = getStrainColors(strainNamesX)
    refVirus = ''

    return virusDict


import pickle

#-----------------------------------------------------------------------------
# Note: virus is a dictionary, the key is the strain name, the value is the position list of mutations (SNPs)
# print('Virus dictionary2',viruses)

# MDS Plotting
import numpy as np
from matplotlib import pyplot as plt
from sklearn import manifold

def plotMDS(viruses):
    strainNames = list(viruses.keys())
    n = len(strainNames)
    distM = np.zeros([n, n]) 
    distV = []
    
    for i in range(0,len(strainNames)):
        nameA = strainNames[i]
        posA = viruses[nameA]
        for j in range (0,n):
            nameB = strainNames[j]
            posB = viruses[nameB]
            dist =  distJaccard(posA ,posB)
            distV.append(dist)
            distM[i,j]= dist
            
    #-----------------------------------------------------------------------------
    import random
    random.seed(9000)
    
    mds = manifold.MDS(n_components = 2, max_iter = 500, dissimilarity = "precomputed", n_jobs = 1)
    pos = mds.fit(distM).embedding_
    
    # model = manifold.TSNE(metric='precomputed')
    # pos = model.fit_transform(distM) 
    # not good as MDS
    
    # import umap
    # U = umap.UMAP(metric='precomputed')
    # pos= U.fit_transform(distM)
    # NOT as good as MDS
    
    # shortProteinNames=[]
    # Plot the points
    plt.figure(figsize = (10,10))
    # m=len(strainNames)
    # colors = ['C{}'.format(i%10) for i in range(m)]
    
    # markerSizes=m*[2] #Illustration for realative lengths of amino acid sequences
    # v=[0.02*j for j in markerSizes]
    
    colorsStrains,labelStrains  = getStrainColors(strainNames)
    
    for i in range(len(pos)):
        plt.plot(pos[i, 0], pos[i, 1],marker = 'o',c=colorsStrains[i],markersize = 20)
        plt.text(pos[i, 0]+0.02, pos[i, 1] + 0.02, labelStrains[i], fontsize = 8,fontweight = 'bold') 
    
    plt.xticks(fontsize = 14, fontweight = 'bold')
    plt.yticks(fontsize = 14, fontweight = 'bold')
    plt.xlabel('Similitude longitude', labelpad = 5,fontsize = 12,fontweight = 'bold')
    plt.ylabel('Similitude latitude',labelpad = 5,fontsize = 12,fontweight = 'bold')
    plt.show()


# Get these virus trains that have at least 2 transmisstion
'''
snpLists=[]
transmits = 5
for key,value in virusDict.items():
    nexts = gh.getAllChildren(value,virusDict)
    if len(nexts.items())>=transmits:
        print(key,value)
        print('Number of nexts:',len(nexts.items()))
        snpLists.append(value)
    #print(firstChild)
    print('Total',len(virusDict))
    print('high SNP List',snpList)
'''
# for snp in snpLists:
#     nexts = gh.getAllChildren(snp,virusDict)
#     print(nexts.items())#,len(children.items()))
#     print('Numbers',len(nexts.items()))
#     print('\n')

#Get all nexts and plot them according the snp that is high transmitting mutations
def plotHighSNPs(snp,virusDict,genome_mutations):
    # snp=[8781, 28143]
    nexts = gh.getAllChildren(snp,virusDict)
    n = str(len(nexts.items()))
    
    strainNames = []
    strainIds = []
    positionData = []
    
    for key, positions in nexts.items():
        virusNameT=key.split('|')
        #strainIds =  key #virusNameT[0]+'|'+virusNameT[1]+'|'+virusNameT[3] #Do not use ids, too large
        strainIds.append(key)
        strainName = virusNameT[0]+'|'+virusNameT[1]+'|'+virusNameT[3] #Do not use ids, too large
        strainNames.append(strainName)
        positionData.append(positions)
    
    colorsStrains,labelStrains  = getStrainColors(strainNames)
    # print(strainNames)
    # print(colorsStrains)
    # mutationInfo = 'of ' + n +' high mutations ' +str(snp)
    mutationInfo =  n + 'high SNPs'+str(genome_mutations)
    gg.eventPlotSNPs(mutationInfo,strainNames,positionData,colorsStrains,lineWidth=4)
    
    return strainIds

def plotHighSNPsX(nexts,virusDict,genome_mutations):
    # snp=[8781, 28143]
    # nexts = gh.getAllChildren(snp,virusDict)
    # n = str(len(nexts.items()))
    
    strainNames = []
    strainIds =[]
    positionData = []
    
    for key, positions in nexts.items():
        virusNameT=key.split('|')
        # strainIds =  key #virusNameT[0]+'|'+virusNameT[1]+'|'+virusNameT[3] #Do not use ids, too large
        strainIds.append(key)
        strainName = virusNameT[0]+'|'+virusNameT[1]+'|'+virusNameT[3] #Do not use ids, too large
        strainNames.append(strainName)
        positionData.append(positions)
    
    colorsStrains,labelStrains  = getStrainColors(strainNames)
    # print('colors before plotting',colorsStrains)
    # print('strain names before plotting',strainNames)
    # print('All children:',nexts)
    
    mutationInfo = str(genome_mutations)
    mutationInfo =  'of SNPs '+str(genome_mutations)
    # print(mutationInfo)
    gg.eventPlotSNPs(mutationInfo,strainNames,positionData,colorsStrains,lineWidth=4)
    
    return strainIds

def getStrainIDSNPs(snp,virusDict):
    # snp=[8781, 28143]
    nexts = gh.getAllChildren(snp,virusDict)

    strainNames = []
    strainIds =[]
    positionData = []
    for key, positions in nexts.items():
        virusNameT=key.split('|')
        # strainIds =  key #virusNameT[0]+'|'+virusNameT[1]+'|'+virusNameT[3] #Do not use ids, too large
        strainIds.append(key)
        strainName = virusNameT[0]+'|'+virusNameT[1]+'|'+virusNameT[3] #Do not use ids, too large
        strainNames.append(strainName)
        positionData.append(positions)
    
    colorsStrains,labelStrains  = getStrainColors(strainNames)
    
    return strainIds
    
# plot all these who have children
# Get these virus trains that have at least 2 transmisstion
# We need to know the exact nucleotide on the reference genome and mutant nucleotide on the virus isolated
# according to the MSA position identified from the sequence Id and SNP location list

# Get the positions and NTs by SNP positions in the MSA file.
#
def getGenomePositions(msaName,seqId,snps):
    # msaName = 'GenomesGISAID_2019-nCoV_03022020_EA_MSA.txt'
    # refSeqId = '2019-nCoV|WH01|NC_045512|2020-01-05'
    # snps = [8789,28151] # This is the mismatched positions,i.e., posMSA, from aligned file
    
    genomeRef_PosNTs = {}
    genomeVirus_PosNTs = {}
    genome_mutations = {}
    snpsRef = []
    
    # print('SNP checked:',snps)
    for posMSA in snps:
        
        [posGenomeRef,ntRef] = gh.mapMSA2RefGenome(msaName,posMSA) #only find ref in the MSA file one #YES check it
        snpsRef.append(posGenomeRef)
        
        # print(seqId,posMSA,'Reference genome position:'+str(posGenomeRef),'NT:'+ntRef)
        genomeRef_PosNTs[posGenomeRef] = ntRef    
    
        # seqId = strainIds[0] #Assume all the mutation in the same position are the same mutation
        # print('posMSA?',posMSA)
        [posGenome,nt] = gh.mapMSA2VirusGenome(msaName,seqId,posMSA)
        # print(seqId,posMSA,'Virus genome position:'+str(posGenome),'NT:'+nt)
        genomeVirus_PosNTs[posGenome]= nt
    
        ntMutation = ntRef+'->'+nt
        genome_mutations[posGenomeRef] = ntMutation
    # print(genomeRef_PosNTs,genomeVirus_PosNTs)
    return genomeRef_PosNTs,genomeVirus_PosNTs,genome_mutations,snpsRef

#TEST: the virus record with deletion of three nucleotides actual genome position [1665-1667] deletion
'''
#2019-nCoV\Data\Alignment\Align0323\GenomesGISAID_SARS-CoV-2_EU_MSA.txt
msaName = 'GenomesGISAID_SARS-CoV-2_EU_MSA.txt'
seqId = 'UK|200690756|EPI_ISL_414044|2020-02-08'
snps =[6666, 13535]
genomeRef_PosNTs,genomeVirus_PosNTs,genome_mutations,snpsRef = getGenomePositions(msaName,seqId,snps)
print(genomeRef_PosNTs,genomeVirus_PosNTs,genome_mutations,snpsRef)
'''

#==============================================================================
refSeqId = '2019-nCoV|WH01|NC_045512|2020-01-05'
snpLists = []
transmits = 12

import csv
def msaSNP2Genome(msaNames):
    # msaName = 'GenomesGISAID_SARS-CoV-2_03062020_US_MSA.txt'
    snpRecords = []
    print('TESTX')
    for msaName in msaNames:
        virusDict = getVirusSNPs(msaName)
        print('Total genomes:',len(virusDict))
        
        for seqId,snps in virusDict.items():
            genomeRef_PosNTs, genomeVirus_PosNTs, genome_mutations, snpsRef = getGenomePositions(msaName,seqId,snps)
            print('msaName',msaName)
            print('SNPs'+str(snps)+' on reference genome:',refSeqId, genomeRef_PosNTs)
            print('SNPs'+str(snps)+' on virus genome:',seqId, genomeVirus_PosNTs)   
            record = {}
            record['seqId'] = seqId
            record['snps'] = snps # This is the positions in align entries in the clustal file
            record['snpsRef'] = snpsRef
            record['refNTs'] = genomeRef_PosNTs # This is the reference locations from the aligned positions
            record['mutatedNTs'] = genomeVirus_PosNTs
            mutationMap = gh.mapMutations(genomeRef_PosNTs,genomeVirus_PosNTs)
            record['mapNTs'] = mutationMap
            snpRecords.append(record) # A list of record dictionary
        print('Total genomes:',len(virusDict))
        
    snpRecords = gh.getAllNextsFromRecords(snpRecords) #ad numNexts (number of desenstants)
    
    pickle.dump(snpRecords, open('./snpRecords/snpRecords_%s2020_%d.pkl'%(SubDate,count), 'wb'))
    # records = pickle.load(open( "SARS-CoV-2_SNPs.pkl", "rb" ))

    csv_columns = ['seqId','snps', 'snpsRef','refNTs','mutatedNTs','mapNTs','numNexts']
    csv_file = './snpRecords/snpRecords_%s2020_%d.csv'%(SubDate,count)
    try:
        with open(csv_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writeheader()
            for record in snpRecords:
                writer.writerow(record)
    except IOError:
        print("I/O error")  
    

def nextsSNPs(msaName,virusDict,transmits):
    snpLists = []
    for seqId,snps in virusDict.items():
        # Find all the sequence snp records in the virusDict that are children of current sequence 
        nexts = gh.getAllChildren(snps,virusDict)
        
        nexts[seqId] = snps # include the transmission strain itself
        
        if len(nexts.items())>=transmits:
            print('This is a high transmission:\n',seqId,snps)
            print('Number of children:',len(nexts.items()))
            snpLists.append(snps)
            
            # msaName = gh.msaNameByStrainId(seqId)
            genomeRef_PosNTs,genomeVirus_PosNTs,genome_mutations = getGenomePositions(msaName,seqId,snps)
            print('SNPs'+str(snps)+' on reference genome:',refSeqId, genomeRef_PosNTs)
            print('SNPs'+str(snps)+' on virus genome:',seqId, genomeVirus_PosNTs)
            
            plotHighSNPsX(nexts,virusDict,genome_mutations)  
            
        print('Total genomes:',len(virusDict))
        print('high SNP List',snpLists)
        return snpLists


