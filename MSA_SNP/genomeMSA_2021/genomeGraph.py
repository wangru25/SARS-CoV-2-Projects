# -*- coding: utf-8 -*- 
""" 
* @Author: Changchuan Yin
* @Date: 2020-04-06 17:01:44 
* @Last Modified by:   Rui Wang 
* @Last Modified time: 2020-04-06 17:01:44  
* @Contact: cyin1@uic.edu
* @Contact: wangru25@msu.edu 
"""


import matplotlib.pyplot as plt
import numpy as np

def eventPlot(virus,colorLabels,positionData,lineWidth=2):
    m = len(positionData)
    
    # virus = '2019-nCoV'
    # colorLabels = ['SLCoV/RaTG13','SLCoV/ZXC21','SLCoV/ZC45','SCoV/Tor2','MERS-CoV','SLCoV/RsSHC014','SLCoV/Rs3367','SLCoV/Shaanxi2011']
    positionData = np.array(positionData)
    
    colors1 = ['C{}'.format(i) for i in range(m)]
    
    lineSize = m * [0.55]#, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]                                 
    linewidths1 =m * [lineWidth]#, 0.9, 0.9, 0.9, 1, 1,1]
    
    fig = plt.figure(figsize = (10,8))
    # [markerline, stemlines, baseline] = plt.stem(x, nduAll, '-.')
    # plt.setp(markerline, 'markerfacecolor', 'b')
    # plt.setp(baseline, 'color', 'r', 'linewidth', 2)
    
    plt.eventplot(positionData, color=colors1, linelengths = lineSize,linewidths = linewidths1)    
    # Provide the title for the spike raster plot
    
    # plot.title('Spike raster plot')
    plt.xticks(fontsize = 12, fontweight = 'bold')
    plt.yticks(fontsize = 12, fontweight = 'bold')
    plt.yticks(np.arange(m), tuple(colorLabels),rotation = 45)
    # plt.ylabel('CoVs',fontsize = 16,fontweight = 'bold')
    # plt.legend(colorLabels, fontsize = 12,bbox_to_anchor=(0., 1.0, 1., .10), loc = 3,ncol = 3, mode="expand", borderaxespad=0.5)
    plt.xlabel('Nucleotide position '+virus,fontsize = 14,fontweight = 'bold')
    plt.show()


def eventPlot2(virus,colorLabels,positionData, colors1,lineWidth=2):
    m =len(positionData)
    
    positionData = np.array(positionData)
    
    # lineSize = [0.4, 0.3, 0.2, 0.8, 0.5, 0.6, 0.7, 0.9]    
    lineSize = m * [0.8]#, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]                                 
    linewidths1 = m * [lineWidth]#, 0.9, 0.9, 0.9, 1, 1,1]
    
    fig = plt.figure(figsize = (10,10))
    # [markerline, stemlines, baseline] = plt.stem(x, nduAll, '-.')
    # plt.setp(markerline, 'markerfacecolor', 'b')
    # plt.setp(baseline, 'color', 'r', 'linewidth', 2)
    
    plt.eventplot(positionData, color=colors1, linelengths = lineSize,linewidths = linewidths1)    
    # Provide the title for the spike raster plot
    
    #plot.title('Spike raster plot')
    plt.xticks(fontsize=10, fontweight = 'bold')
    plt.yticks(fontsize=7, fontweight = 'bold')
    plt.yticks(np.arange(m), tuple(colorLabels),rotation=45)
    # plt.ylabel('CoVs',fontsize = 16,fontweight = 'bold')
    # plt.legend(colorLabels, fontsize = 12,bbox_to_anchor=(0., 1.0, 1., .10), loc = 3,ncol = 3, mode="expand", borderaxespad=0.5)
    plt.xlabel('Nucleotide position '+virus,fontsize = 12,fontweight = 'bold')
    
    plt.show()


def eventPlotT(virus,colorLabels,positionData, colors1,lineWidth=2):
    m = len(positionData)
    
    # virus = '2019-nCoV'
    # colorLabels = ['SLCoV/RaTG13','SLCoV/ZXC21','SLCoV/ZC45','SCoV/Tor2','MERS-CoV','SLCoV/RsSHC014','SLCoV/Rs3367','SLCoV/Shaanxi2011']
    positionData = np.array(positionData)
    
    # lineSize = [0.4, 0.3, 0.2, 0.8, 0.5, 0.6, 0.7, 0.9]    
    lineSize = m * [0.8]#, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]                                 
    linewidths1 =m * [lineWidth]#, 0.9, 0.9, 0.9, 1, 1,1]
    
    fig = plt.figure(figsize = (10,10))
    # [markerline, stemlines, baseline] = plt.stem(x, nduAll, '-.')
    # plt.setp(markerline, 'markerfacecolor', 'b')
    # plt.setp(baseline, 'color', 'r', 'linewidth', 2)
    
    plt.eventplot(positionData, color=colors1, linelengths = lineSize,linewidths = linewidths1)    
    # Provide the title for the spike raster plot
    
    # plot.title('Spike raster plot')
    plt.xticks(fontsize = 14, fontweight = 'bold')
    plt.yticks(fontsize = 14, fontweight = 'bold')
    # plt.yticks(np.arange(m), tuple(colorLabels),rotation=45)
    plt.ylabel('SARS-CoV-2 mutations',fontsize = 14,fontweight = 'bold')
    # plt.legend(colorLabels, fontsize = 12,bbox_to_anchor=(0., 1.0, 1., .10), loc = 3,ncol = 3, mode="expand", borderaxespad=0.5)
    plt.xlabel('Nucleotide position '+virus,fontsize = 14,fontweight = 'bold')
    plt.show()    


def eventPlotY(yLabel,virus,colorLabels,positionData, colors1,lineWidth=2):
    m = len(positionData)
    
    # virus = '2019-nCoV'
    # colorLabels = ['SLCoV/RaTG13','SLCoV/ZXC21','SLCoV/ZC45','SCoV/Tor2','MERS-CoV','SLCoV/RsSHC014','SLCoV/Rs3367','SLCoV/Shaanxi2011']
    positionData = np.array(positionData)
    
    # lineSize = [0.4, 0.3, 0.2, 0.8, 0.5, 0.6, 0.7, 0.9]    
    lineSize = m * [0.8]#, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]                                 
    linewidths1 =m * [lineWidth]#, 0.9, 0.9, 0.9, 1, 1,1]
    
    fig = plt.figure(figsize = (10,10))
    # [markerline, stemlines, baseline] = plt.stem(x, nduAll, '-.')
    # plt.setp(markerline, 'markerfacecolor', 'b')
    # plt.setp(baseline, 'color', 'r', 'linewidth', 2)
    
    plt.eventplot(positionData, color=colors1, linelengths = lineSize,linewidths = linewidths1)    
    # Provide the title for the spike raster plot
    
    # plot.title('Spike raster plot')
    plt.xticks(fontsize = 14, fontweight = 'bold')
    plt.yticks(fontsize = 14, fontweight = 'bold')
    # plt.yticks(np.arange(m), tuple(colorLabels),rotation=45)
    plt.ylabel('SARS-CoV-2 in '+yLabel,fontsize = 8,fontweight = 'bold')
    # plt.legend(colorLabels, fontsize = 12,bbox_to_anchor=(0., 1.0, 1., .10), loc = 3,ncol = 3, mode="expand", borderaxespad=0.5)
    plt.xlabel('Nucleotide position '+virus,fontsize = 14,fontweight = 'bold')
    # Display the spike raster plot
    plt.show()    

def eventPlotSNPs(virus,seqIds,positionData, colors1, lineWidth=2):
    m = len(positionData)
    # print('StrainName in plot',seqIds)
    # print('Data coming:',positionData)
    positionData = np.array(positionData)
    
    lineSize = m * [0.8]                               
    linewidths1 = m * [lineWidth]
    
    fig = plt.figure(figsize = (10,10))
    # [markerline, stemlines, baseline] = plt.stem(x, nduAll, '-.')
    # plt.setp(markerline, 'markerfacecolor', 'b')
    # plt.setp(baseline, 'color', 'r', 'linewidth', 2)
    
    plt.eventplot(positionData, color = colors1, linelengths = lineSize, linewidths = linewidths1)    
    # Provide the title for the spike raster plot
    
    # plot.title('Spike raster plot')
    plt.xticks(fontsize = 14, fontweight = 'bold')
    plt.yticks(fontsize = 14, fontweight = 'bold')
    plt.yticks(np.arange(m), tuple(seqIds),rotation=45)
    # plt.ylabel('CoVs',fontsize = 16,fontweight = 'bold')
    # plt.legend(colorLabels, fontsize = 12,bbox_to_anchor=(0., 1.0, 1., .10), loc = 3,ncol = 3, mode="expand", borderaxespad=0.5)
    plt.xlabel('Nucleotide position '+virus, fontsize = 14,fontweight = 'bold')
    # Display the spike raster plot
    plt.show()
