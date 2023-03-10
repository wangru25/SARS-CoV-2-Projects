# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 18:49:18 2021

@author: yutah
"""

import os as os
import sys as sys
from os import listdir
from os.path import isfile, join
#date = sys.argv[1]
#end = int(sys.argv[2])


def run(date_all):
    for date in date_all:
        inpath = "%s/"%date
        files = [f for f in listdir(inpath) if isfile(join(inpath, f))]
        end = int(len(files)/2)
        clustalWdownload_path = "../../genomeMSA_2021/clustalW/%s/"%(date)
        
        for idx in range(1, end+1):
            clustal_name = "GenomesGISAID_SARS-CoV-2_%s2021_MSA_%s.txt"%(date, idx)
            fasta = inpath + "%s_gisaid_hcov-19_2021_%s.fasta"%(str(idx), date)
            if not os.path.exists(clustalWdownload_path + clustal_name):
                job_name = "%s_gisaid_hcov-19_2021_%s"%(str(idx), date)
                fasta = inpath + "%s_gisaid_hcov-19_2021_%s.fasta"%(str(idx), date)
                print("sending ", fasta)
                
                #os.system("python clustalo.py --email hozumiyu@msu.edu --stype dna  --asyncjob --title %s --sequence %s "%(job_name, fasta) )
                os.system("python clustalo.py --email hfeng5@crimson.ua.edu --stype dna  --asyncjob --title %s --sequence %s "%(job_name, fasta) )
            else:
                print(fasta, "already finished")
    
if __name__ == "__main__":
    folders = [x[0] for x in os.walk('./')]
    SubDates = []
    bad_name = ["./", './clust', './jobId', './Unporcessed_Result']
    for folder in folders:
        if folder not in bad_name:
            SubDates.append(folder[2:])
    run(SubDates)
    
    
    '''
    date = sys.argv[1]
    date_all = date.split(',')
    run(date_all)
    '''