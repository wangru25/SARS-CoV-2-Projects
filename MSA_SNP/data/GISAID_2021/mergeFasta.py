# -*- coding: utf-8 -*-
'''
@Author: Rui Wang
@Date: 2020-05-01 11:23:57
@LastModifiedBy: Rui Wang
@LastEditTime: 2020-05-15 01:19:43
@Email: wangru25@msu.edu
@FilePath: /38_Influ/data/GISAID/mergeFasta.py
@Description: Based on Dr. Changchuan Yin's code. Contact: yin1@uic.com
'''
import os
import sys
from Bio import SeqIO
import pycountry

'''
    merge fast files
    split data into set of 100 fasta files
    Output also has country.txt
    Make sure to check the country name and alpha-2-code matches
    If the name don't match, change def countryCode(self, countryT)
    If the country has multiple words, such as New Zealand, South Africa. Often misidentified
'''

class mergeFasta:
    def __init__(self, inFile, indexFile, numRecords, SubDate):
        self.inFile = inFile         # 'gisaid_cov2020_sequences_0416.fasta'
        self.indexFile = indexFile
        self.numRecords = numRecords
        self.SubDate = SubDate
        self.working_dir = './%s/'%self.SubDate
        self.refFasta = 'NC_045512.fasta'
        self.countryReplace = []
        
    def countryCode(self, countryT):
        countries = {}
        for country in pycountry.countries:
            countries[country.name] = country.alpha_2
        codeT = ''
        
        #if citys are given instead, make a list of citys for the country
        china = ['Beijing', 'Chongqing', 'Fujian', 'Fuzhou', 'Guangzhou', 'Guangdong', 'Harbin', 'Hangzhou', 'Jiujiang', 'Kashgar', 'Liaoning', "Lu`an", "Lu'an", 'Lishui', 'Shaoxing', 'Shanghai', 
                 'Shulan', 'Sichuan', 'Wuhan', 'Yichun', 'Yingtan', 'Yunnan', 'Zhejiang']
        southafrica = ['South Africa', 'SouthAfrica']
        hongkong = ['HongKong']
        saudia = ['SaudiArabia']
        SaintBarthelemy = ['SaintBarthelemy']
        newzealand = ['New', 'New Zealand', 'NewZealand']
        romania = ['Bucuresti', 'Teleorman']
        france = ['Felis', 'French_Guiana']
        germany = ['Berlin']
        dominican = ['DominicanRepublic']
        uk = ['Wales', 'England','Northern Ireland' , 'Scotland', 'Northern', 'uk', 'NorthernIreland']
        southkorea = ['SouthKorea']
        unitedarab = ['UnitedArabEmirates']
        CI = ["CotedIvoire", "Coted'Ivoire"]
        reunion = ['Reunion']
        sri = ['SriLanka']
        costa = ['CostaRica']
        faroe = ['FaroeIslands']
        czech = ["CzechRepublic"]
        saintmartin= ['SaintMartin']
        SierraLeone = ['SierraLeone']
        GuineaBissau = ['GuineaBissau']
        BQ = ["StEustatius"]
        CW = ['Curacao']
        GQ = ["EquatorialGuinea"]
        BF = ["BurkinaFaso"]
        SV = ["ElSalvador"]
        MK = ["NorthMacedonia"]
        PG = ["PapuaNewGuinea", "PapuaNewGuinea"]
        BA = ["BosniaandHerzegovina"]
        FR = ["FrenchGuiana"]
        XK = ["Kosovo"]
        netherlands = ["SintMaarten", 'Curacao', 'Aruba']
        for country,code in countries.items():
            if countryT.upper() in country.upper():  #The key in the countries are full country name
                codeT = code
                if countryT in southafrica:
                    codeT = 'ZA'
                elif countryT in uk:
                    codeT = 'UK'
                elif countryT in newzealand:
                    codeT = 'NZ'
                elif countryT in PG:
                    codeT = "PG"
                elif countryT in netherlands:
                    codeT = "NL"
                break
            elif countryT in hongkong:
                codeT = 'HK'
            elif countryT in southkorea:
                codeT = 'KR'
            elif countryT == 'USA':
                codeT = 'US'
            elif countryT == 'Vietnam':
                codeT = 'VN'
            elif countryT == 'Czech Republic':
                codeT == 'CZ'
            elif countryT == 'Lithuania':
                codeT = "LT"
            elif countryT == 'Estonia':
                codeT = "EE"
            elif countryT in china:
                codeT = 'CN'
            elif countryT == 'Bosnia-and-Herzegovina':
                codeT = 'BA'
            elif countryT in southafrica:   #south aftrica
                codeT = 'ZA'
            elif countryT in newzealand:
                codeT = 'NZ'
            elif countryT in romania: #city in Romaina
                codeT = 'RO'
            elif countryT in uk:
                codeT = 'UK'
            elif countryT == 'DRC':
                codeT = 'CD'
            elif countryT == 'Bahrein':
                codeT = 'BH'
            elif countryT in france:
                codeT = 'FR'
            elif countryT in saudia:
                codeT= 'SA'
            elif countryT in germany:
                codeT = 'DE'
            elif countryT in unitedarab:
                codeT = 'AE'
            elif countryT in dominican:
                codeT = 'DO'
            elif countryT in reunion:
                codeT = 'RE'
            elif countryT in sri:
                codeT = 'LK'
            elif countryT in costa:
                codeT = 'CR'
            elif countryT in faroe:
                codeT = 'FO'
            elif countryT in czech:
                codeT = "CZ"
            elif countryT in SaintBarthelemy:
                codeT = "BL"
            elif countryT in saintmartin:
                codeT = "MF"
            elif countryT in CI:
                codeT = "CI"
            elif countryT in BQ:
                codeT = "BQ"
            elif countryT in GQ:
                codeT = "GQ"
            elif countryT in BF:
                codeT = "BF"
            elif countryT in SV:
                codeT = "SV"
            elif countryT in MK:
                codeT = "MK"
            elif countryT in PG:
                codeT = "PG"
            elif countryT in BA:
                codeT = "BA"
            elif countryT in FR:
                codeT = "FR"
            elif countryT in XK:
                codeT = "XK"
            elif countryT in netherlands:
                codeT = "NL"
            elif countryT in SierraLeone:
                codeT = "SL"
            elif countryT in GuineaBissau:
                codeT = 'GW'
            else:
                codeT = countryT
        self.countryReplace.append(countryT + '->' + codeT + '\n')
        return codeT
                
    def writeFastaFile(self, headers, seqs, fileName):
        country =  list(set(self.countryReplace))
        with open('country.txt', 'w') as handle:
            for k in range((len(country))):
                handle.write(country[k])
            
        with open(fileName, 'w') as handle:
            for j in range(0, len(seqs)):
                header = headers[j]
                handle.write('>'+header+'\n')
                seq = seqs[j]
                width = 70
                for i in range(0, len(seq), width):
                    seqT = str(seq[i:i+width])
                    seqT = seqT.upper()
                    handle.write(seqT + '\n')
    
    def writeIndexFile(self, ids, prefix):
        newIds = []  
        for i in range(len(ids)):
            new = '%s_%s_%s'%(i, self.SubDate, prefix)   
            newIds.append(new)
        with open(self.working_dir + str(prefix) + '_' + self.indexFile,'w') as f:
            f.write('Full_Name,Map_Name\n')
            for i in range(len(ids)):
                f.write('%s,%s\n'%(ids[i], newIds[i]))
        return newIds

    def partition(self, ls, size):
        """
        Returns a new list with elements
        of which is a list of certain size.

            >>> partition([1, 2, 3, 4], 3)
            [[1, 2, 3], [4]]
        """
        return [ls[i:i+size] for i in range(0, len(ls), size)]

    def GISIADdata(self):
        seqFastaT = self.working_dir + 'temp_' + self.inFile
        openFile = self.working_dir + self.inFile
        with open(openFile, 'r') as file:
            filedata = ''
            filedata = file.read()
            filedata = filedata.replace('hCoV-19/', '')
            filedata = filedata.replace('mink/', '')
            filedata = filedata.replace('hCov-19|', '')
            filedata = filedata.replace('hCoV|', '')
            filedata = filedata.replace('/2020', '')
            filedata = filedata.replace('/', '|')
            filedata = filedata.replace('Hong ', 'Hong')
            filedata = filedata.replace('South ', 'South')
            filedata = filedata.replace('Dominican ', 'Dominican')
            filedata = filedata.replace('United ', 'United')
            filedata = filedata.replace('Arab ', 'Arab')
            filedata = filedata.replace('Sri ', 'Sri')
            filedata = filedata.replace('Faroe ', 'Faroe')
            filedata = filedata.replace('Costa ', 'Costa')
            filedata = filedata.replace('Czech ', 'Czech')
            filedata = filedata.replace('Northern ', 'Northern')
            filedata = filedata.replace('New ', 'New')
            filedata = filedata.replace('Saint Barthelemy', 'SaintBarthelemy')
            filedata = filedata.replace(' ', '')
        with open(seqFastaT, 'w') as file2:
            file2.write(filedata)

        genomeIds = []
        genomeSeqs = []
        numSamples = len(list(SeqIO.parse(seqFastaT, "fasta")))

        idx0 = ''
        seq0 = ''
        for record in SeqIO.parse(self.refFasta, "fasta"):
            idx0 = record.id
            seq0 = record.seq
            # genomeIds.append(idx0)
            # genomeSeqs.append(seq0)
        #print(idx0)
        temp_idA = []
        temp_seq = []
        for record in SeqIO.parse(seqFastaT, "fasta"):
            idx = record.id
                
            seq = record.seq
            country = idx.split('|')[0]
            idA = idx.replace(country, self.countryCode(country))
            idAT = idA.split('|')
            # dateId = idAT[3].replace('-','')
            # print(dateId)
            seqT = seq[500:29000]
            if 'NNNNN' in seqT:
                #print('Too many N in',idx)
                continue
            elif len(idAT) < 4:
                print('Problem with header format in ', idx)
                #continue
            #elif len(idAT[3].replace('-','')) != 8:
                #print('Problem with date format in ', idx)
            #    #continue
            #elif idAT[3].replace('-','').isdigit() is False:
                #print('Problem with date format in ', idx)
                #continue
            else: #only store the records having full 4 headers sections and having no 'NNNNN'
                genomeIds.append(idA) 
                genomeSeqs.append(seq)
                if 'KR' in idA:
                    temp_idA.append(idA)
                    temp_seq.append(seq)
                elif 'ZA' in idA:
                    temp_idA.append(idA)
                    temp_seq.append(seq)
                elif 'HK' in idA:
                    temp_idA.append(idA)
                    temp_seq.append(seq)
        
        if numSamples > self.numRecords:
            splitGenomeIds = self.partition(genomeIds, self.numRecords)
            splitGenomeSeqs = self.partition(genomeSeqs, self.numRecords)
        else: 
            splitGenomeIds = [genomeIds]
            splitGenomeSeqs = [genomeSeqs]
        for i in range(len(splitGenomeIds)):
            prefix = i + 1
            genomeIds = splitGenomeIds[i]
            genomeIds.append(idx0)
            genomeSeqs = splitGenomeSeqs[i]
            genomeSeqs.append(seq0)
            idxList = []
            genomeSeqs_nonDup = []
            genomeIds_nonDup = sorted(set(genomeIds), key=genomeIds.index)
            for i in range(len(genomeIds_nonDup)):
                if genomeIds_nonDup[i] in genomeIds:
                    idxList.append(i)
            for i in range(len(idxList)):
                seq_T = genomeSeqs[i]
                genomeSeqs_nonDup.append(seq_T)
            genomeIds_nonDup = self.writeIndexFile(genomeIds_nonDup, prefix)
            fileName = self.working_dir + str(prefix) + '_' + self.inFile
            self.writeFastaFile(genomeIds_nonDup, genomeSeqs_nonDup, fileName)
        
        
        splitGenomeIds = [temp_idA]
        splitGenomeSeqs = [temp_seq]
        for i in range(len(splitGenomeIds)):
            prefix = 'temp'
            genomeIds = splitGenomeIds[i]
            genomeIds.append(idx0)
            genomeSeqs = splitGenomeSeqs[i]
            genomeSeqs.append(seq0)
            idxList = []
            genomeSeqs_nonDup = []
            genomeIds_nonDup = sorted(set(genomeIds), key=genomeIds.index)
            for i in range(len(genomeIds_nonDup)):
                if genomeIds_nonDup[i] in genomeIds:
                    idxList.append(i)
            for i in range(len(idxList)):
                seq_T = genomeSeqs[i]
                genomeSeqs_nonDup.append(seq_T)
            genomeIds_nonDup = self.writeIndexFile(genomeIds_nonDup, prefix)
            fileName = self.working_dir + str(prefix) + '_' + self.inFile
            self.writeFastaFile(genomeIds_nonDup, genomeSeqs_nonDup, fileName)
        #os.remove(seqFastaT)

def main(inFile, indexFile, numRecords, SubDate):
    MergeFasta = mergeFasta(inFile, indexFile, numRecords, SubDate)
    MergeFasta.GISIADdata()
    
if __name__ == "__main__":
    SubDates = sys.argv[1]
    #SubDates = '0516,0517A'
    countryList = []
    SubDate_vec = SubDates.split(',')
    for SubDate in SubDate_vec:
        #SubDate = '0531'
        print("Analyzing " + SubDate)
        numRecords = 100
        inFile = 'gisaid_hcov-19_2021_%s.fasta'%SubDate
        indexFile = 'mapIndex_%s.csv'%SubDate
        absolute_dir = '/home/rui/Dropbox/Linux_Backup/MSU/1_Training/38_Influ/data/GISAID/%s'%SubDate
        main(inFile, indexFile, numRecords, SubDate)
        file = open('country.txt', 'r')
        country = file.readlines()
        file.close()
        countryList += country
        file2 = open('country2.txt', 'w')
        for line in countryList:
            file2.writelines(line)
        file2.close()
        # os.chdir('/home/rui/Dropbox/Linux_Backup/MSU/1_Training/38_Influ/data/GISAID/%s'%SubDate)
        # os.system('mafft --auto --clustalout --reorder "1_gisaid_cov2020_sequences_%s.fasta" > "../../../genomeMSA/clustalW/GenomesGISAID_SARS-CoV-2_%s2020_MSA_1.txt"'%(SubDate, SubDate))



