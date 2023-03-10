# -*- coding: utf-8 -*-
'''
@Author: Rui Wang
@Date: 2020-05-03 08:27:49
@LastModifiedBy: Rui Wang
LastEditTime: 2021-06-13 11:37:14
@Email: wangru25@msu.edu
FilePath: /38_Influ/analysis/Animation/proteinMain.py
@Description: 
'''
import numpy as np
import pandas as pd 
from collections import Counter
import sys
import os
import re
import os
import sys
import json
import csv
from csv import DictReader
csv.field_size_limit(100000000)

def addQuotesDigits(strX):
    strY = re.sub(r"(\d+)", r'"\1"', strX)
    return strY

def snpRecordsFromCSV(csvSNPName):
    with open(csvSNPName, 'r') as recordsObj:
        dict_reader = DictReader(recordsObj)
        records = list(dict_reader)
        return records  

# Get a list of SNPs from the SNP csv data
# Each item is the SNP for the virus isolate
def parseSNPsRecords(records):
    SNPs = []
    # print(records)
    for record in records:
        seqId = record['seqId']
        refNTs = record['mapNTs']
        # print('SNP:',refNTs)   #Note: This is a dictionary like string, but it is a string
        refNTs = refNTs.replace('\'','\"')
        refNTs = addQuotesDigits(refNTs)
        if len(refNTs) > 1: # get rid of empty rows
            SNP = {}
            SNP['seqId'] = seqId
            SNP['date']= seqId.split('|')[3]
            refNTsDict = json.loads(refNTs) #Convert dictionary like string to dictionary
            # print(list(refNTsDict.keys())[-1])
            # print(refNTsDict)
            
            for pos, mutation in refNTsDict.items():                    
                SNP[pos] = mutation
                # print(SNP)
                # print('\n')
                if len(refNTsDict) > 1:
                    if int(pos) == int(list(refNTsDict.keys())[-1]):
                        SNPs.append(SNP)
                else:
                    SNPs.append(SNP)
    return SNPs

def uniqueSNPs(SNPs):
    snpNames = ['A->T','A->C','A->G','T->A','T->C','T->G','C->T','C->A','C->G','G->T','G->C','G->A']
    uniqueSNPsX = []
    for SNP in SNPs:
        for pos, snpType in SNP.items():
            if snpType in snpNames:
                uniqueSNPsX.append((pos,snpType))           
    uniqueSNPsX0 = list(set(uniqueSNPsX)) 
    return uniqueSNPsX0

def non_uniqueSNPs(SNPs):
    snpNames = ['A->T','A->C','A->G','T->A','T->C','T->G','C->T','C->A','C->G','G->T','G->C','G->A']
    uniqueSNPsX = []
    for SNP in SNPs:
        
        for pos, snpType in SNP.items():
            if snpType in snpNames:
                uniqueSNPsX.append((pos,snpType))           
    return uniqueSNPsX
    
def string2list(s):
    '''
    A string to list. For example:
    str([12,13,14,16]) ---> [12,13,14,16]
    '''
    S = []
    s = s.split(', ')
    for i in range(len(s)):
        if len(s) == 1:
            s[i] = s[i].replace('[','')
            s[i] = s[i].replace(']','')
            snpref = int(s[i])
            S.append(snpref)
        else:
            if i == 0:
                s[i] = s[i].replace('[','')
                snpref = int(s[i])
                S.append(snpref)
            elif i == len(s) - 1:
                s[i] = s[i].replace(']','')
                snpref = int(s[i])
                S.append(snpref)
            else:
                snpref = int(s[i])
                S.append(snpref)
    return S

def replaceString(instring):
    # instring = instring.replace(", 'm",'')
    instring = instring.replace(',','')
    instring = instring.replace('   ','')
    instring = instring.replace('  ','')
    instring = instring.replace(' ','')
    instring = instring.replace('\'','')
    instring = instring.replace(", 'm",'')
    # instring = instring.replace("'m",'')
    return instring

def getIndex(inList, inNum):
    return [i for i, x in enumerate(inList) if x == inNum]

def customLog(inList, inNum):
    return np.log(np.asarray(inList)+2)/np.log(inNum)

# a = [1,2,3,1000, 2000, 3000]
# print(customLog(a,1.3))

def getDNASequencesFasta(fastaFileName,idd):
    '''
    fastaFileName='./Protein/SARS-CoV-2_proteome.fasta'
    idd = 'nsp3'
    proteinSeq = getDNASequencesFasta(fastaFileName,idd)
    print(proteinSeq[0])
    '''
    from Bio.Seq import Seq
    from Bio import SeqIO
    seqR =''
    for record in SeqIO.parse(fastaFileName, "fasta"):
        #print(record.id) #header
        idx = record.id
        seq = record.seq
        if idx == idd:
            seqR = seq
    return seqR

def createProteinSite(inString, inNum):
    '''
    Input:
    '''
    return inString.replace(':',str(inNum))

# print(createProteinSite('L:V',727))

def removeString(inList):
    newList = []
    for i in range(len(inList)):
        new = int(inList[i][1:-1])
        newList.append(new)
    return newList

def hIndex(citations):
    """
    :type citations: List[int]
    :rtype: int
    """
    if not citations:
        return 0
    citations.sort()
    h_index = 0
    size = len(citations)
    for i in range(size):
        try_index = size - i
        if citations[i] >= try_index and (i < 1 or citations[i - 1] <= try_index):
            h_index = max(h_index , try_index)
    return h_index
#==============================================================================================

def prepareData(inFile,Date):
    idList = ['spike_surface_glycoprotein','3C-like_proteinase','envelope_protein','nucleocapsid_phosphoprotein','RNA-dependent-polymerase',
    'nsp3','nsp9', 'endoRNAse', "exonuclease", 'leader_protein_nsp1', 'nsp2', 'nsp4', 'nsp6', 'nsp7', 'nsp8', 'nsp9', 'nsp10', 'helicase',
    'methyltransferase', 'nsp11', 'ORF3a_protein', 'ORF6_protein', 'ORF7a_protein', 'ORF7b_protein', 'ORF8_protein', 'ORF10_protein','membrane_glycoprotein']

    for nameId in idList:
        proteinFile = open(inFile)
        if os.path.exists('./animation_%s'%Date) is False:
            os.mkdir('./animation_%s'%Date)

        if os.path.exists('./animation_%s/SNP'%Date) is False:
            os.mkdir('./animation_%s/SNP'%Date)

        extractFile = open("./animation_%s/SNP/%s_SNPs_%s.txt"%(Date,nameId,Date), 'w')
        fastaFile = 'SARS-CoV-2_proteome.fasta'
        # idd = '%s'%nameId
        idd = "%s"%nameId
        seqR = getDNASequencesFasta(fastaFile,idd)
        # print(seqR)
        # print(nameId)
        for line in proteinFile:
            if "%s"%nameId in line:
                tmp = int(re.findall(".*'position': (.*), 'mutation'.*",line)[0]) - 1
                if tmp < 0 or tmp > len(seqR)-1:
                    pass
                else:
                    print(seqR[tmp]+line)
                    extractFile.write('%s,'%seqR[tmp]+line)

        proteinFile.close()
        extractFile.close()




#================================================================================================
def prepareFreq(inFile, Date):
    records = snpRecordsFromCSV(inFile)
    SNPs = parseSNPsRecords(records)
    unique_SNPs = non_uniqueSNPs(SNPs)
    
    records = snpRecordsFromCSV(inFile)
    SNPs = parseSNPsRecords(records)
    nonunique_SNPs = non_uniqueSNPs(SNPs)
    frequencyDict = Counter(nonunique_SNPs)

    if os.path.exists('./animation_%s/MutationSummary'%Date) is False:
        os.mkdir('./animation_%s/MutationSummary'%Date)

    idList = ['spike_surface_glycoprotein','3C-like_proteinase','envelope_protein','nucleocapsid_phosphoprotein','RNA-dependent-polymerase',
    'nsp3','nsp9', 'endoRNAse', "exonuclease", 'leader_protein_nsp1', 'nsp2', 'nsp4', 'nsp6', 'nsp7', 'nsp8', 'nsp10', 'helicase',
    'methyltransferase', 'nsp11', 'ORF3a_protein', 'ORF6_protein', 'ORF7a_protein', 'ORF7b_protein', 'ORF8_protein', 'ORF10_protein','membrane_glycoprotein']



    dictList = {'spike_surface_glycoprotein': 'Spike','3C-like_proteinase':'3CL','envelope_protein': 'Envelope','nucleocapsid_phosphoprotein':'Nucleocapsid',
                'RNA-dependent-polymerase':'RdRp','nsp3': 'NSP3','nsp9':'NSP9', 'endoRNAse':'endoRNAse', "exonuclease":'Exonuclease', 'leader_protein_nsp1':'NSP1',
                'nsp2':'NSP2', 'nsp4':'NSP4', 'nsp6':'NSP6', 'nsp7':'NSP7', 'nsp8':'NSP8', 'nsp10':'NSP10', 'helicase':'Helicase',
                'methyltransferase':'2’-O-ribose MTases', 'nsp11':'NSP11', 'ORF3a_protein':'ORF3a', 'ORF6_protein':'ORF6', 'ORF7a_protein':'ORF7a', 
                'ORF7b_protein':'ORF7b', 'ORF8_protein':'ORF8', 'ORF10_protein':'ORF10','membrane_glycoprotein':'Membrane'}

    for nameId in idList:
        extractFile =  open('./animation_%s/SNP/%s_SNPs_%s.txt'%(Date,nameId,Date))
        dictnameID = dictList['%s'%nameId]

        frequencyFile = open('./animation_%s/MutationSummary/%s_frequency_%s.csv'%(Date,dictnameID,Date), 'w')
        frequencyFile.write('mutation_site,protein_site(actual),Total_frequency\n')
        for line in extractFile:
            if line[0] == line[-6]:
                # print(line[8+2:c+2])
                gene_pos = int(re.findall(".*'pos': (.*?), 'mutation'.*",line)[0])
                pro_pos = int(re.findall(".*'position': (.*), 'mutation'.*",line)[0])
                gene_mutation = re.findall(".*'mutation'(.*)},{'protein'.*",line)[0]
                if int(frequencyDict[('%s'%gene_pos, eval(gene_mutation.replace(": ",'')))]) != 0:
                    frequencyFile.write('%d%s,%s,%d\n'%(gene_pos,gene_mutation,
                    createProteinSite(line[-6:-3], pro_pos),
                    frequencyDict[('%s'%gene_pos, eval(gene_mutation.replace(": ",'')))]
                    ))
        frequencyFile.close()

def removeSymbol(inFile, outFile):
    readFile = open(inFile)
    writeFile = open(outFile,'w')
    lines = readFile.readlines()

    for item, line in enumerate(lines):
        if line[-4] != '#' and line[-4] != ',' and line[-4] != '_':
            writeFile.write(line)


if __name__ == "__main__":
    Date = sys.argv[1]  #09112020
    inFile_1 =  '/Users/rui/Dropbox/Linux_Backup/MSU/1_Training/38_Influ/analysis/proteinPlot/Protein/Protein_SNPs_%s.txt'%Date
    outFile = '/Users/rui/Dropbox/Linux_Backup/MSU/1_Training/38_Influ/analysis/proteinPlot/Protein/Protein_SNPs_%s_real.txt'%Date
    removeSymbol(inFile_1, outFile)
    
    inFile_2 = 'snpRecords_%s.csv'%Date
    prepareData(outFile,Date)
    prepareFreq(inFile_2, Date)
