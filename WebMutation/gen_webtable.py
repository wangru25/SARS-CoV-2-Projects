# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2021-06-01 13:07:07
LastModifiedBy: Rui Wang
LastEditTime: 2021-09-17 14:12:38
Email: wangru25@msu.edu
FilePath: /38_Influ/analysis/WebMutation/gen_webtable.py
Description: 
'''
import os
from datetime import timezone
from os import write
import sys
import math
import heapq
import numpy as np
import pandas as pd

def get_index(lst=None, item=''):
    return [index for (index,value) in enumerate(lst) if value == item]

def extractRBDGene(Date):
    working_dir = '../Animation/animation_%s/MutationSummary/'%Date
    binding_dir = './data/'
    spikeFile = pd.read_csv(working_dir + 'Spike_frequency_%s.csv'%Date) 
    bindingFile = pd.read_csv(binding_dir + '6M0J_Spike_C4_0.results', header = None)
    writeFile_gene = open(working_dir + 'Spike_RBD_Gene_%s.csv'%Date,'w')
    writeFile_gene.write('mutation_site,protein_site(actual),BFE changes,Total_frequency\n')

    for i in range(spikeFile.shape[0]):
        a = spikeFile['protein_site(actual)'][i][1:-1]
        if int(a) >= 333 and int(a) <= 526:
            if spikeFile['protein_site(actual)'][i][-1] == spikeFile['protein_site(actual)'][i][0]:
                BFE = 0.0
            else:
                index_1 = get_index(bindingFile[5].tolist(), int(a))
                index_2 = get_index(bindingFile[6].tolist(), str(spikeFile['protein_site(actual)'][i][-1]))
                unionSet = [i for i in index_1 if i in index_2]
                BFE = bindingFile[8][unionSet[0]]
            writeFile_gene.write('%s,%s,%.4f,%s\n'%(spikeFile['mutation_site'][i],spikeFile['protein_site(actual)'][i],BFE,spikeFile['Total_frequency'][i]))
    
def extractRBDResidue(Date):
    working_dir = '../Animation/animation_%s/MutationSummary/'%Date
    binding_dir = './data/'
    RBDFile = pd.read_csv(working_dir + 'Spike_RBD_Gene_%s.csv'%Date) 
    bindingFile = pd.read_csv(binding_dir + '6M0J_Spike_C4_0.results', header = None)
    writeFile_residue = open(working_dir + 'Spike_RBD_NonDegenate_%s.csv'%Date,'w')
    writeFile_residue.write('Wild type,Residue position,Mutant type,Total_frequency,BFE changes,protein_site(actual)\n')

    mutantList = RBDFile['protein_site(actual)'].tolist()
    mutantList_Nondup = list(dict.fromkeys(mutantList))
    for item in mutantList_Nondup:
        if item[0] != item[-1]:
            index_1 = get_index(mutantList, item)
            BFE = RBDFile['BFE changes'][index_1[0]]
            freq = 0
            for i in index_1:
                freq +=  RBDFile['Total_frequency'][i]
            writeFile_residue.write('%s,%s,%s,%d,%.4f,%s\n'%(item[0],item[1:-1],item[-1],freq,BFE,item))


def extractSResidue(Date):
    working_dir = '../Animation/animation_%s/MutationSummary/'%Date
    binding_dir = './data/'
    RBDFile = pd.read_csv(working_dir + 'Spike_frequency_%s.csv'%Date) 
    writeFile_residue = open(working_dir + 'Spike_frequency_NonDegenate_%s.csv'%Date,'w')
    writeFile_residue.write('Wild type,Residue position,Mutant type,Total_frequency,protein_site(actual)\n')

    mutantList = RBDFile['protein_site(actual)'].tolist()
    mutantList_Nondup = list(dict.fromkeys(mutantList))
    for item in mutantList_Nondup:
        if item[0] != item[-1]:
            index_1 = get_index(mutantList, item)
            freq = 0
            for i in index_1:
                freq +=  RBDFile['Total_frequency'][i]
            writeFile_residue.write('%s,%s,%s,%d,%s\n'%(item[0],item[1:-1],item[-1],freq,item))

def create_disruptive_table(Date):
    antibody_namefile = pd.read_csv('./data/antibodies.csv',header=None)
    antibody_list = []
    for i in range(antibody_namefile.shape[0]):
        antibody_list.append(antibody_namefile[0][i][0:4])

    RBD_residue_file = pd.read_csv('../Animation/animation_%s/MutationSummary/Spike_RBD_NonDegenate_%s.csv'%(Date,Date))
    write_antibody_disrupt = open('./data/distructed_antibodies.csv','w')
    write_antibody_disrupt.write('Mutation,BFE changes,counts,rates\n')

    for i in range(len(RBD_residue_file['protein_site(actual)'].tolist())):
        mutation = RBD_residue_file['protein_site(actual)'][i]
        print(mutation)
        bfe = RBD_residue_file['BFE changes'][i]
        antibodies_disrupt_list = []
        for antibody in antibody_list:
            antibodies_corr =  [txt for txt in os.listdir('../CoMutation/antibodies_BFE/results') if txt.startswith('%s_RBD_C3'%antibody)]
            antibody_file = pd.read_csv('../CoMutation/antibodies_BFE/results/%s'%antibodies_corr[0],header=None)
            for i in range(antibody_file.shape[0]):
                if antibody_file[4][i] == mutation[0] and antibody_file[6][i] == mutation[-1] and int(antibody_file[5][i]) == int(mutation[1:-1]):
                    if antibody_file[8][i] < -0.3:
                        antibodies_disrupt_list.append(antibody)
        write_antibody_disrupt.write('%s,%s,%d,%.2f\n'%(mutation,bfe,len(antibodies_disrupt_list), 100*len(antibodies_disrupt_list)/antibody_file.shape[0]))


def refine_disruptive_table(date):
    antibodies_file = pd.read_csv('./data/destructed_antibodies.csv')
    stats_file = pd.read_csv('../Animation/animation_%s/MutationSummary/Spike_RBD_NonDegenate_%s.csv'%(date,date))
    write_file = open('./data/destructed_antibodies_%s.csv'%date,'w')
    write_file.write('Mutation,BFE changes,counts,rates\n')

    mutation_date = stats_file['protein_site(actual)'].tolist()
    mutation_disruptive = antibodies_file['Mutation'].tolist()
    print(antibodies_file)
    BFE_changes = stats_file['BFE changes'].tolist()
    for i in range(len(mutation_date)):
        if mutation_date[i] not in mutation_disruptive:
            write_file.write("%s,%.4f,%d,%.2f\n"%(mutation_date[i],BFE_changes[i],int(0),0.00))
        else:
            idx = get_index(mutation_disruptive,mutation_date[i])[0]
            write_file.write("%s,%.4f,%d,%.2f\n"%(mutation_date[i],BFE_changes[i],antibodies_file['counts'].tolist()[idx],antibodies_file['rates'].tolist()[idx]))
    write_file.close()


def creat_webtable(top_num, date):
    stats_file = pd.read_csv('../Animation/animation_%s/MutationSummary/Spike_RBD_NonDegenate_%s.csv'%(date,date))
    antibodies_file = pd.read_csv('./data/destructed_antibodies_%s.csv'%date)
    write_file = open('./data/web_table.txt','w')
    write_file.write('  <table class="sortable" width="100%" border="0" cellspacing="1" style="width:77%; margin-left:auto;margin-right:auto;">\n')
    write_file.write('  <tr><th>Rank/Mutation</th><th>Worldwdie counts</th><th>BFE change</th></tr>\n')

    mutation_list = stats_file['protein_site(actual)'].tolist()
    BFE_changes = stats_file['BFE changes'].tolist()
    antibodies_mutation = antibodies_file['Mutation'].tolist()
    antibodies_counts = antibodies_file['counts'].tolist()
    antibodies_rates = antibodies_file['rates'].tolist()
    freq_list = stats_file['Total_frequency'].tolist()
    freq_sortidx = np.argsort(np.array(freq_list)).tolist()[::-1]

    BFE_changes_top = []
    antibodies_top = []
    for i in range(top_num):
        idx = freq_sortidx[i]
        mutation = mutation_list[idx]
        tmp_idx = get_index(antibodies_mutation,mutation)[0]
        BFE_changes_top.append(BFE_changes[idx])
        if isinstance(antibodies_counts[tmp_idx],int):
            antibodies_top.append(antibodies_counts[tmp_idx])  
        else:
            antibodies_top.append(-math.inf)  
    BFE_sortidx = np.argsort(np.array(BFE_changes_top)).tolist()[::-1]
    antibodies_sortidx = np.argsort(np.array(antibodies_top)).tolist()[::-1]

    BFE_rank = np.zeros((1,len(BFE_sortidx)))
    for i in range(len(BFE_sortidx)):
        BFE_rank[0,BFE_sortidx[i]] = i+1

    antibodies_rank = np.zeros((1,len(antibodies_sortidx)))
    for i in range(len(antibodies_sortidx)):
        antibodies_rank[0,antibodies_sortidx[i]] = i+1

    for i in range(top_num):
        idx = freq_sortidx[i]
        mutation = mutation_list[idx]
        tmp_idx = get_index(antibodies_mutation,mutation)[0]

        BFE_idx = BFE_rank[0,i]
        antibodies_idx = antibodies_rank[0,i]
        write_file.write('    <tr><td>%s</td><td>%d</td><td>%d</td><td>%.4f</td><td>%d</td><td>%s</td><td>%s</td><td>%d</td></tr>\n'%(mutation,freq_list[idx],i+1,BFE_changes[idx],BFE_idx,antibodies_counts[tmp_idx],antibodies_rates[tmp_idx],antibodies_idx))
    write_file.write('  </table>')
    write_file.close()



date = sys.argv[1]

extractRBDGene('%s'%date) 
extractRBDResidue('%s'%date) 
# create_disruptive_table(date)
refine_disruptive_table(date)
stats_file = pd.read_csv('../Animation/animation_%s/MutationSummary/Spike_RBD_NonDegenate_%s.csv'%(date,date))
total_num = stats_file.shape[0]
creat_webtable(total_num,'%s'%date)



# extractSResidue('%s'%date) 
