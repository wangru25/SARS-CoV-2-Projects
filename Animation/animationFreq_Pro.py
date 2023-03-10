# -*- coding: utf-8 -*-
'''
@Author: Rui Wang
@Date: 2020-05-11 13:20:40
@LastModifiedBy: Rui Wang
LastEditTime: 2022-08-22 16:25:17
@Email: wangru25@msu.edu
FilePath: /38_Influ/analysis/Animation/animationFreq_Pro.py
@Description: 
'''
import re
import os
import sys
import numpy as np
import pandas as pd
from collections import Counter
import plotly
import plotly.tools as tls 
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots
from datetime import datetime


class dataReader:
    def __init__(self, inFile, outFile):
        self.inFile = pd.read_csv(inFile)
        self.outFile = outFile

    def mapDict(self):
        keys = ['[266:805]', '[806:2719]', '[2720:8554]', '[8555:10054]', '[10055:10972]', '[10973:11842]', '[11843:12091]', '[12092:12685]', '[12686:13024]', '[13025:13441]', '[13442:16236]',
        '[16237:18039]', '[18040:19620]', '[19621:20658]', '[20659:21552]', '[13442:13480]', '[21563:25384]', '[25393:26220]', '[26245:26472]', '[26523:27191]', '[27202:27387]', '[27394:27759]', '[27756:27887]',
        '[27894:28259]', '[28274:29533]', '29558:29674']
        values = ['NSP1', 'NSP2', 'NSP3', 'NSP4', '3CL', 'NSP6', 'NSP7', 'NSP8', 'NSP9', 'NSP10', 'RdRp', 'Helicase', "Exonuclease", 'endoRNAse', 
        "2’-O-ribose MTases", 'NSP11','Spike', 'ORF3a', 'Envelope', 'Membrane', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'Nucleocapsid', 'ORF10']
        nameMap = dict(zip(keys,values))
        return nameMap

    def lableDict(self, inName):
        keys = ['NSP1', 'NSP2', 'NSP3', 'NSP4', '3CL', 'NSP6', 'NSP7', 'NSP8', 'NSP9', 'NSP10', 'RdRp', 'Helicase', "Exonuclease", 'endoRNAse', 
        "2’-O-ribose MTases", 'NSP11','Spike', 'ORF3a', 'Envelope', 'Membrane', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'Nucleocapsid', 'ORF10']
        values = range(26)
        nameMap = dict(zip(keys,values))
        return nameMap[inName]

    def getDate(self, seqId):
        '''
        Input:  seqId: TW|170|EPI_ISL_420084|2020-03-21|66
        Output: dateId: 20200321
        '''
        seqId = seqId.split('|')
        dateId = int(seqId[3].replace('-',''))
        return dateId

    def string2list(self, inString):
        '''
        A string to list. For example:
        str([12,13,14,16]) ---> [12,13,14,16]
        '''
        outString = []
        inString = inString.split(', ')
        for i in range(len(inString)):
            if len(inString) == 1:
                inString[i] = inString[i].replace('[','')
                inString[i] = inString[i].replace(']','')
                snpref = int(inString[i])
                outString.append(snpref)
            else:
                if i == 0:
                    inString[i] = inString[i].replace('[','')
                    snpref = int(inString[i])
                    outString.append(snpref)
                elif i == len(inString) - 1:
                    inString[i] = inString[i].replace(']','')
                    snpref = int(inString[i])
                    outString.append(snpref)
                else:
                    snpref = int(inString[i])
                    outString.append(snpref)
        return outString

    def getIndex(self,inList, item):
        return([index for (index,value) in enumerate(inList) if value == item])

    def getDateList(self):
        seqId = self.inFile['seqId']
        dateList = []
        for i in range(self.inFile.shape[0]):
            dateList.append(self.getDate(seqId[i]))
        dateList = np.sort(list(set(dateList)))
        return dateList

    def newFile(self, date):
        seqId = self.inFile['seqId']
        numSample = self.inFile.shape[0]
        selectList = []
        for i in range(numSample):
            if self.getDate(seqId[i]) <= date:
                selectList.append(i)
        newFile = self.inFile.iloc[selectList]
        newFile.index = range(len(newFile))
        return newFile
    
    def getSNPList(self, newFile):
        snpRef = newFile['snpsRef']
        numSample = newFile.shape[0]
        snpList = self.string2list(snpRef[0])
        for i in range(numSample-1):
            snpList_init = self.string2list(snpRef[i+1])
            snpList += snpList_init
        # snpList = np.sort(list(set(snpList))) 
        return snpList
    
    def getFullSNP_Nondup(self):
        snpRef = self.inFile['snpsRef']
        numSample = self.inFile.shape[0]
        snpFullList = self.string2list(snpRef[0])
        for i in range(numSample-1):
            snpList_init = self.string2list(snpRef[i+1])
            snpFullList += snpList_init
        snpFullList = np.sort(list(set(snpFullList))) 
        return len(snpFullList)

    def freqByDate(self, date):
        newFile = self.newFile(date)
        snpList = self.getSNPList(newFile)
        freqDict = Counter(snpList)
        return freqDict

    def dict2list(self, inDict):
        keyRangeList = []
        for i in range(len(inDict)):
            keyRange = list(inDict.keys())[i]
            keyRange = keyRange.replace('[','').replace(']','').split(':')
            keyRange = list(map(int,keyRange))
            keyRangeList.append(keyRange)
        return keyRangeList

    def mapProteinName(self, position):
        keyRangeList = self.dict2list(self.mapDict())
        proteinName = 'Others'
        for i in range(len(keyRangeList)):
            if position >= keyRangeList[i][0] and position <= keyRangeList[i][1]:
                proteinName =  list(self.mapDict().values())[i]
            else:
                continue
        return proteinName

    def dateList(self, beginDate, endDate):
        date_l=[int(datetime.strftime(x,'%Y-%m-%d').replace("-","")) for x in list(pd.date_range(start=beginDate, end=endDate, freq='60D'))]
        # string_l = [datetime.strftime(x,'%b %d %Y') for x in list(pd.date_range(start=beginDate, end=endDate, freq='10D'))]
        return date_l


    def freqWrite(self, position, freqFile):
        readFile = pd.read_csv(freqFile)
        positionList = []
        for i in range(readFile.shape[0]):
            positionList.append(int(readFile['mutation_site'][i].split(": ")[0]))
        # freq = Counter(positionList)[position]
        indexList = self.getIndex(positionList,position)
        if len(indexList) == 0:
            NTs_1 = 'None'; PTs_1 = 'None'; freq_1 = 0
            NTs_2 = 'None'; PTs_2 = 'None'; freq_2 = 0
            NTs_3 = 'None'; PTs_3 = 'None'; freq_3 = 0
            total_freq = freq_1+ freq_2 + freq_3
        elif len(indexList) == 1:
            NTs_1 = readFile['mutation_site'][indexList[0]].replace(": ",""); PTs_1 = readFile['protein_site(actual)'][indexList[0]]; freq_1 = readFile['Total_frequency'][indexList[0]]
            NTs_2 = 'None'; PTs_2 = 'None'; freq_2 = 0
            NTs_3 = 'None'; PTs_3 = 'None'; freq_3 = 0
            total_freq = freq_1+ freq_2 + freq_3
        elif len(indexList) == 2:
            NTs_1 = readFile['mutation_site'][indexList[0]].replace(": ",""); PTs_1 = readFile['protein_site(actual)'][indexList[0]]; freq_1 = readFile['Total_frequency'][indexList[0]]
            NTs_2 = readFile['mutation_site'][indexList[1]].replace(": ",""); PTs_2 = readFile['protein_site(actual)'][indexList[1]]; freq_2 = readFile['Total_frequency'][indexList[1]]
            NTs_3 = 'None'; PTs_3 = 'None'; freq_3 = 0
            total_freq = freq_1+ freq_2 + freq_3
        elif len(indexList) == 3:
            NTs_1 = readFile['mutation_site'][indexList[0]].replace(": ",""); PTs_1 = readFile['protein_site(actual)'][indexList[0]]; freq_1 = readFile['Total_frequency'][indexList[0]]
            NTs_2 = readFile['mutation_site'][indexList[1]].replace(": ",""); PTs_2 = readFile['protein_site(actual)'][indexList[1]]; freq_2 = readFile['Total_frequency'][indexList[1]]
            NTs_3 = readFile['mutation_site'][indexList[2]].replace(": ",""); PTs_3 = readFile['protein_site(actual)'][indexList[2]]; freq_3 = readFile['Total_frequency'][indexList[2]]
            total_freq = freq_1+ freq_2 + freq_3
        else: 
            print("Error!")
        return NTs_1,NTs_2,NTs_3,PTs_1,PTs_2,PTs_3,freq_1,freq_2,freq_3, total_freq
            


    def writeData(self,Date):
        # dateList = self.getDateList()
        dateList = self.dateList('2020-01-01','%s-%s-%s'%(Date[4:8], Date[0:2], Date[2:4]))
        print(dateList)
        writeFile = open(self.outFile,'w')
        # writeFile.write('Protein,Date,Nucleotide Position,ln(Frequency),Frequency,Protein Label,Unique Single Mutations,Num of Samples\n')
        writeFile.write('Protein,Date,Nucleotide Position,ln(Frequency),Frequency,Mutation Type 1/Freq,Mutation Type 2/Freq,Mutation Type 3/Freq,Unique Single Mutations,Num of Samples,Protein Label\n')
        uniqueSNPs = self.getFullSNP_Nondup()
        numSamples = self.inFile.shape[0]
        values = ['NSP1', 'NSP2', 'NSP3', 'NSP4', '3CL', 'NSP6', 'NSP7', 'NSP8', 'NSP9', 'NSP10', 'RdRp', 'Helicase', "Exonuclease", 'endoRNAse', 
        "2’-O-ribose MTases", 'NSP11','Spike', 'ORF3a', 'Envelope', 'Membrane', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'Nucleocapsid', 'ORF10']
        for i in range(len(dateList)):
            date = dateList[i]
            freqDict = self.freqByDate(date)
            positionList = list(freqDict.keys())
            for j in range(len(positionList)):
                position = positionList[j]

                if self.mapProteinName(position) in values:
                    freqFile = './animation_%s/MutationSummary/%s_frequency_%s.csv'%(Date,self.mapProteinName(position),Date)
                    readFile = pd.read_csv(freqFile)
                    tmpList = []
                    for k in range(readFile.shape[0]):
                        tmpList.append(int(readFile['mutation_site'][k].split(": ")[0]))
                    if position in tmpList:

                        NTs_1,NTs_2,NTs_3,PTs_1,PTs_2,PTs_3,freq_1,freq_2,freq_3,total_freq = self.freqWrite(position,freqFile)

                        if int(total_freq) != 0:
                    
                            writeFile.write('%s,%s,%d,%.4f,%d,%s-(%s): %s,%s-(%s): %s,%s-(%s): %s,%d,%d,%d\n'%(self.mapProteinName(position), date, position, np.log(total_freq),total_freq,
                            NTs_1,PTs_1,freq_1,
                            NTs_2,PTs_2,freq_2,
                            NTs_3,PTs_3,freq_3,
                            uniqueSNPs, numSamples, self.lableDict(self.mapProteinName(position))
                            ))
        writeFile.close()


def datelist(beginDate, endDate):
    date_l=[int(datetime.strftime(x,'%Y-%m-%d').replace("-","")) for x in list(pd.date_range(start=beginDate, end=endDate, freq='60D'))]
    # string_l = [datetime.strftime(x,'%b %d %Y') for x in list(pd.date_range(start=beginDate, end=endDate, freq='10D'))]
    return date_l


def getPeriodCSV(inCSV, Date):
    df = pd.read_csv(inCSV)
    dateList = datelist('2020-01-01','%s-%s-%s'%(Date[4:8], Date[0:2], Date[2:4]))
    maxDate = max(df['Date'].tolist())
    
    saveIndex = []
    for i in range(df.shape[0]):
        if df['Date'][i] in dateList:
            saveIndex.append(i)
        else:
            if df['Date'][i] == maxDate:
                saveIndex.append(i)

    newRecord = df.iloc[saveIndex]
    newRecord.index = range(len(newRecord))
    return newRecord
    # newRecord.to_csv('./data/Patients/snpRecords_07022020_age_valid.csv', index=False)

def plotHTML(inCSV, Date):
    import base64
    # image_filename = 'gisaid_logo.png'
    image_filename = 'logo_gisaid.png'
    plotly_logo = base64.b64encode(open(image_filename, 'rb').read())


    # df = pd.read_csv(inCSV)
    df = getPeriodCSV(inCSV,Date)    
    
    fig = px.scatter(df, x = 'Nucleotide Position', y = 'ln(Frequency)', color = 'Protein Label', animation_frame = "Date", animation_group="Nucleotide Position",
    hover_data=["Nucleotide Position","Frequency","Date","Protein","Mutation Type 1/Freq","Mutation Type 2/Freq","Mutation Type 3/Freq"], range_y=[0,20],
    color_continuous_scale=[
    "rgb(160,200,255)","yellow","rgba(13,143,129,0.7)", "rgba(119,74,1750,0.8)", "rgba(235,74,64,0.8)","rgb(62,109,178)", "rgb(160,200,255)", "rgb(253,174,107)", "rgba(119,74,1750,0.8)", 
    "rgb(160,200,255)","black","rgba(13,143,129,0.7)", "rgba(119,74,1750,0.8)", "rgba(235,74,64,0.8)","rgb(62,109,178)", "rgb(160,200,255)", "black", "rgb(260,200,255)",
    "rgb(160,200,255)","yellow","rgba(13,143,129,0.7)", "rgba(119,74,1750,0.8)", "rgba(235,74,64,0.8)","rgb(62,109,178)", "rgb(160,200,255)", "rgb(253,174,107)"])

    
    # tickvals = [10513.5, 14839, 23473.5, 26358.5, 28903.5]
    # ticktext = ['3CL', 'RdRp', 'S', 'E', 'N']
    tickvals = [535.5, 5637, 9304.5, 10513.5, 11407.5, 14839, 23473.5, 25806.5, 26358.5, 28903.5,17138]
    ticktext = ['NSP1', 'NSP3', 'NSP4', '3CL', 'NSP6', 'RdRp', 'S', 'ORF3a', 'E', 'N','Helicase']

    fig.add_layout_image(
    dict(
        source='data:image/png;base64,{}'.format(plotly_logo.decode()),
        # source='<a href="www.google.com/"><img src="data:image/png;base64,{}" /> </a>'.format(plotly_logo.decode()),
        xref="paper", yref="paper",
        x=0.83, y=0.90,
        sizex=0.14, sizey=0.14,
        xanchor="left", yanchor="bottom"
        )
    )

    fig.update_layout(template = 'simple_white', 
    # title = "%d Single Mutations in %d hCoV-19 Genomes"%(df['Unique Single Mutations'][0], df['Num of Samples'][0])
    title = ''
    + "<span style='font-size: 20px;'> %d Single Mutations in %d hCoV-19 Genomes</span>"%(df['Unique Single Mutations'][0], df['Num of Samples'][0])
    + '<br>'
    + "<span style='font-size: 17px;'> Relevant link: <a href='https://weilab.math.msu.edu/MutationAnalyzer/'> Analysis of S protein RBD mutations</a> </span>"
    + "</br>",
    # + "<span style='font-size: 11px;'> enabled by data from</span>"
    # + "<span style='font-size: 11px;'> <a href='https://www.gisaid.org'> GISAID</a> </span>",  
    title_x=0.55, title_y = 0.9,
    # + '<br>'
    # + "<span style='font-size: 17px;'> enabled by</span>"
    # + "<span style='font-size: 17px;'> <a href='https://www.gisaid.org'> GISAID's</a> </span>"
    # + "<span style='font-size: 17px;'> data which is subjected to</span>"
    # + '<span style="font-size: 17px;"> <a href="https://www.gisaid.org/registration/terms-of-use/"> Terms and Conditions</a> </span>',
    # title_x=0.5, title_y = 0.93,
    xaxis=dict(ticktext=ticktext,tickvals= tickvals),xaxis_title="GISAID data provided on this website is subject to GISAID’s <a href='https://www.gisaid.org/registration/terms-of-use/'> Terms and Conditions</a> <a href='https://users.math.msu.edu/users/weig/MutationSummary.zip'> <[Download Summary]> </a> ", yaxis_title="ln(Frequency)",
    font=dict(family='Times New Roman', size=14),
    # legend=dict(font_size = 14),
    # autosize = True
    # title = {
    # 'text': ''
    # + "<span style='font-size: 16px;'> %d Single Mutations in %d hCoV-19 Genomes</span>"%(df['Unique Single Mutations'][0], df['Num of Samples'][0])
    # + "<span style='font-size: 11px;'> enabled by data from</span>"
    # + "<span style='font-size: 11px;'> <a href='https://www.gisaid.org'> GISAID</a> </span>",
    # 'xref': 'paper',
    # 'yref': 'paper',
    # 'x' : 0.30, 
    # 'y' : 0.93,
    # 'xanchor': 'left',
    # 'yanchor': 'bottom'
    # }
    # annotations=[
    #     dict(
    #         x = 0.74,
    #         y = 0.90,
    #         xref="paper",
    #         yref="paper",
    #         text="<span style='font-size: 11px;'> enabled by data from</span>",
    #         # text="GISAID data provided on this website is subject to GISAID’s" + '<a href="https://www.gisaid.org/registration/terms-of-use/"> Terms and Conditions</a>',
    #         showarrow=False,
    #         font_size= 14,
    #         xanchor="left", yanchor="bottom",
    #         ax = 0,
    #         ay = -40
    #     )
    # ]
    )
    # fig.update_xaxes(tickfont=dict(family='Rockwell', size=14))
    # fig.update_yaxes(tickfont=dict(family='Rockwell', size=14))
    fig.update(layout_coloraxis_showscale=False)

    plotly.offline.plot(fig, filename='./animation_%s/SARS-CoV-2_Mutation_Tracker_%s_%s.html'%(Date,Date[0:2],Date[2:4]))

def main(inFile, outFile, Date):
    DataReader = dataReader(inFile, outFile)
    DataReader.writeData(Date)

if __name__ == "__main__":
    Date = sys.argv[1]  # 0511
    inFile = 'snpRecords_%s.csv'%Date
    outFile = './animation_%s/animation_%s.csv'%(Date,Date)
    # main(inFile,outFile,Date)
    plotHTML(outFile,Date)


