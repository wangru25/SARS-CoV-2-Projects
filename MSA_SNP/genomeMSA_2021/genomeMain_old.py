# -*- coding: utf-8 -*- 
""" 
* @Author: Changchuan Yin
* @Date: 2020-04-06 09:27:09 
* @Last Modified by:   Rui Wang 
* @Last Modified time: 2020-04-06 09:27:09  
* @Contact: cyin1@uic.edu
* @Contact: wangru25@msu.edu 
"""
import os
import sys
import pickle
import genomeSNPs as gs
import genomeUtility as gh


SubDate = sys.argv[1]   # 0405: Represent for submission date: April 05
count = int(sys.argv[2])
numN = int(sys.argv[3])  # 10: Extract sequences with at least 10 descendants. A nonnegative integer.

msaName = './clustalW/GenomesGISAID_SARS-CoV-2_%s2020_MSA_%d.txt'%(SubDate, count)
msaNames = [msaName] 
#--------------------------------------------------------------------------------------
# snpLists = nextsSNPs(virusDict,transmits)
# virusSNPs={8782: 'C->T', 17747: 'C->T', 17858: 'A->G', 18060: 'C->T', 28144: 'T->C'}
# virusSNPs = list({v['pos']:v for v in snpss}.values()) #remove duplicated dictionaries
# print(virusSNPs)

def uniqueSNPfromRecords(records,numN):
    records = pickle.load(open("./snpRecords/snpRecords_%s2020_%d.pkl"%(SubDate, count), "rb" ))
    snpss = []
    for record in records:
        if record['numNexts'] >= numN:       
            print('High',record)
            # print(record['mapNTs'],record['numNexts'])
            mapNTs = record['mapNTs']
            snps = createSNP(mapNTs)
            snpss = snpss+snps
    virusSNPs = list({v['pos']:v for v in snpss}.values()) #remove duplicated dictionaries      
    return  virusSNPs    
    # print(snpss)

def snpRecordBySeqId(seqId,records):
    # seqId = 'DE|BavPat1|EPI_ISL_406862|2020-01-28'
    result = []
    for record in records:
        seqIdT = record['seqId']
        if seqIdT == seqId:
            result = record
            break
    return result

# Count the SNPS excludes those in the two ends (sequencing erros may occur in the two ends)
def confidentSNPs(snpsRef):
    # snpsRef =[10,20,300,1000,29600]
    confSNPs = []  
    for snp in snpsRef:
        if snp > 200 and snp < 29500:
            confSNPs.append(snp)
                
    numSNPs = len(confSNPs)
    # print(numSNPs)
    return numSNPs,confSNPs

# Function returns N largest elements 
def Nmaxelements(list1, N): 
    final_list = [] 
    for i in range(0, N):  
        max1 = 0
        for j in range(len(list1)):      
            if list1[j] > max1: 
                max1 = list1[j]; 
        list1.remove(max1); 
        final_list.append(max1) 
    return final_list 
    
def maxSNPMutation(records):
    # seqId = 'DE|BavPat1|EPI_ISL_406862|2020-01-28'
    snpLens = []
    numRecords = 0
    # records.remove()
    for record in records:
        seqId = record['seqId'] 
        snpsRef = record['snpsRef'] #this is a list of SNPs in ref genome
        # umSNPs = len(snpsRef)
        numSNPs, confSNPs = confidentSNPs(snpsRef) #only count the SNPs in the middle of the genomes
        if 'EPI_ISL_406592' not in seqId: #do not consider Shenzhen as it was probably mistake in sequencing?
            snpLens.append(numSNPs) 
            numRecords = numRecords+1
    maxSNPs =  max(snpLens) 
    idx = snpLens.index(maxSNPs)
    print('Maxium SNP numbers is at',idx)
    print('Maxium SNP numbers is at',records[idx])
    return numRecords, maxSNPs   

def createSNP(mapNTs):
    snps = []
    for pos,mutation in mapNTs.items():
        snp = {}
        snp['pos'] = pos
        snp['SNP'] = mutation
        snps.append(snp)
    return snps

# Get the records after a certain dateX
def recordsAfter(records,dateX):
    # records = pickle.load(open("./snpRecords/snpRecords_%s2020.pkl"%(SubDate), "rb" ))
    # print('TEST',records[0]['seqId'])
    recordsX=[]
    for record in records:
        seqId = record['seqId'] # 'FR|HF1805|EPI_ISL_414628|2020-03-02'
        seqIds = seqId.split('|')
        date = seqIds[3]
        if date > dateX: # '2020-03-05':
            # print('YES')
            recordsX.append(record)
    return recordsX

# Get the records before a certain dateX
def recordsBefore(records,dateX):
    # records = pickle.load(open("./snpRecords/snpRecords_%s2020.pkl"%(SubDate), "rb" ))
    # print('TEST',records[0]['seqId'])
    recordsX = []
    for record in records:
        seqId = record['seqId'] # 'FR|HF1805|EPI_ISL_414628|2020-03-02'
        seqIds = seqId.split('|')
        date =seqIds[3]
        if date < dateX: # '2020-03-05':
            # print('YES')
            recordsX.append(record)
    return recordsX

# seqId = 'NL|NA8|EPI_ISL_415497|2020-03-10'
# group = 'EU2'
# Get the records after a certain dateX
def recordsGroup(records,group):
    # records = pickle.load(open("./snpRecords/snpRecords_%s2020.pkl"%(SubDate), "rb" ))
    # print('TEST',records[0]['seqId'])
    recordsX = []
    for record in records:
        seqId = record['seqId']#'FR|HF1805|EPI_ISL_414628|2020-03-02'
        if gh.isInGroup(seqId,group):
            recordsX.append(record)
    return recordsX

# records = pickle.load(open("./snpRecords/snpRecords_%s2020.pkl"%(SubDate), "rb" ))
# recordsX = recordsAfter(records,'2020-02-28')
# print('Current',recordsX)
# print('dateTime',seqIds[3])

#-------------------------------------------------------------------
# Write all the SNPs into a csv file, save the SNP records as records.pkl 
# Generate all virus SNP records as snpRecors.pkl and csv file 'SARS-CoV-2_SNPs.csv'
# Get and plot all SNPs on the reference genome from all regions,
# and print the frequenciues

# Save all virus SNP records into a dictionary
mainData = True #True
if mainData:
    gs.msaSNP2Genome(msaNames) # generate snpRecors.pkl and SARS 'SARS-CoV-2_SNPs.csv'
    print('Completed in processing the MSA files')
    
#-------------------------------------------------------------------------
# Frequencies of SNPS
snpFrequency = False #True
if snpFrequency:
    records = pickle.load(open('./snpRecords/snpRecords_%s2020_%d.pkl'%(SubDate, count), "rb" ))
    virusDict = gh.getVirusDictFromRecords(records)
    gs.plotSNPs_noID(virusDict)
    
    freqDict = gh.plotSNPFrequencies(virusDict)
    print('Frequencies',freqDict)
    
    [numRecords, maxSNPs ] = maxSNPMutation(records)
    print('Total genome records and maximum mutations:',numRecords, maxSNPs)

# Frequencies of SNPS
snpFrequency_March = False
if snpFrequency_March:
    records = pickle.load(open("./snpRecords/snpRecords_%s2020_%d.pkl"%(SubDate, count), "rb" ))
    
    records = recordsAfter(records,'2020-03-01') # after 2020-02-20
    records = recordsGroup(records,'EU')

    virusDict = gh.getVirusDictFromRecords(records)
    gs.plotSNPs_noID(virusDict)
    gs.plotSNPs_ID(virusDict)
    
    freqDict = gh.plotSNPFrequencies(virusDict)
    print('Frequencies',freqDict)
    
    [numRecords, maxSNPs ] = maxSNPMutation(records)
    print('Total genome records and maximum mutations:',numRecords, maxSNPs)  
    
#-------------------------------------------------------------------------
# Protein mutation mapping
mapProteins = False
if mapProteins:
    # virusSNPs={8782: 'C->T', 17747: 'C->T', 17858: 'A->G', 18060: 'C->T', 28144: 'T->C'}
    # virusSNPs = list({v['pos']:v for v in snpss}.values()) #remove duplicated dictionaries 
    records = pickle.load(open("./snpRecords/snpRecords_%s2020_%d.pkl"%(SubDate, count), "rb" ))
    
    # records = recordsAfter(records,'2020-03-01') #after 2020-02-20
    # records = recordsGroup(records,'US')
    virusDict = gh.getVirusDictFromRecords(records)
    print(virusDict)
    virusSNPs = uniqueSNPfromRecords(records,numN)
    # print(virusSNPs)
    gh.mapSNPs2Proteins(virusSNPs, numN) #produce 'Protein_SNPs_04042020_numNexts_10.txt'
    
    print('Completed protein mapping')  
    
#--------------------------------------------------------------------------------------
# Plot the nexts of a given strain id
commonSNPs = False
if commonSNPs:
    records = pickle.load(open("./snpRecords/snpRecords_%s2020_%d.pkl"%(SubDate, count), "rb" )) 
    
    # seqId = 'ES|Valencia32|EPI_ISL_420117|2020-03-09|79' # No1
    # seqId = 'ES|Valencia38|EPI_ISL_420123|2020-03-09|74' # No2
    # seqId = 'ES|Valencia42|EPI_ISL_420126|2020-03-12|53' # No3
    # seqId = 'TW|125|EPI_ISL_420082|2020-03-19|27'        # No4   
    # seqId = 'ES|Valencia49|EPI_ISL_420130|2020-03-18|28' # No5
    
    recordQ = snpRecordBySeqId(seqId,records)
    
    # only search the unique records
    virusDictNexts = gh.nextsFromRecords(recordQ,records)
    print('Next SNPs:', virusDictNexts)

    # remove duplicated dictionary with the same snp list
    variants=[]
    uniqueVirusDictNexts={}
    
    for seqId,snp in virusDictNexts.items():
        isSNP = gh.isInList(snp,variants)
        if isSNP == False:
            variants.append(snp)
            uniqueVirusDictNexts[seqId]=snp
    print('Number of unique nexts',uniqueVirusDictNexts.keys())
    
    countries = []
    for seqId,snp in uniqueVirusDictNexts.items():
        ids = seqId.split('|')
        countryCode = ids[0]
        countries.append(countryCode)
        
    # gs.plotSNPs_noID(uniqueVirusDictNexts)
    gs.plotSNPs_ID(virusDictNexts)
    
    countryCodes = gh.getCountryCodes(countries) 
    gs.plotSNPs_Y(uniqueVirusDictNexts,countryCodes)
    # get MDS position all nodes
    mdsCoords = gh.snpMDS(uniqueVirusDictNexts) 
    
    # *******************************************************************
    # Plot connected network, non-duplicated nodes, and some nodes with 
    # closet descandants 
    gh.snpNetwork(uniqueVirusDictNexts,pos = mdsCoords)
    #*******************************************************************
    # Plot connected network (these nodes and their all descandants)
    gh.snpNetworkNexts(uniqueVirusDictNexts,pos = mdsCoords)
    
    print('Number of unique nexts',len(uniqueVirusDictNexts.keys()))
    
    
#--------------------------------------------------------------------------------------
# Get the SNPs for each region(done)
snpRegions = True
if snpRegions:
    
    virusDict_US = gs.getVirusSNPs(msaName)
    
    # viruses = gs.plotSNPs_noID(virusDict_EU)
    # viruses = gs.plotSNPs_noID(virusDict_US) 
    
    # virusDict = gs.plotSNPs(msaName0)
    '''
    virusDict_CN = gs.getVirusSNPs(msaName)
    gs.plotSNPs_noID(virusDict_CN)
    
    virusDict_US = gs.getVirusSNPs(msaName)
    gs.plotSNPs_noID(virusDict_US)
    
    stateName = 'CruiseA'
    virusDict = gs.plotSNPs_subset(msaName,stateName) #commented code in plotSNPs_subset
    '''
#--------------------------------------------------------------------------------------
test5 = False
if test5:
    virusDict =  gs.getVirusSNPs(msaName) 
    transmits = 10
    snpLists = gs.nextsSNPs(msaName,virusDict,transmits)
