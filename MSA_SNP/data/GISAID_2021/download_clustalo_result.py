# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 12:19:16 2021

@author: yutah
"""
import os
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
        


def obtain_all_jobId():
    inpath = "jobId/"
    jobId_ls_pre = [f for f in listdir(inpath) if isfile(join(inpath, f))]
    jobId_ls = []
    for jobId in jobId_ls_pre:
        jobId = jobId.replace(".txt", '')
        jobId_ls.append(jobId)
    
    return jobId_ls

def download_finished_clustatlo(jobId_ls):
    for jobId in jobId_ls:
        os.system("python clustalo.py --polljob --jobid %s"%jobId)
        
def deletejobId(jobId):
    jobpath = "jobId/"
    os.remove(jobpath + jobId + ".txt")
    print("--------------- finished removing jobId --------------")        
    
def deleteResult(jobId):
    inpath = "Unporcessed_Result/"
    downloaded_files = [f for f in listdir(inpath) if isfile(join(inpath, f))]
    
    for files in downloaded_files:
        if jobId in files:
            os.remove(inpath + files)
    print("--------------- finished removing excess files --------------")    
    
def obtain_all_finished(jobId_ls):
    inpath = "Unporcessed_Result/"
    jobpath = "jobId/"
    downloaded_files = [f for f in listdir(inpath) if isfile(join(inpath, f))]
    print(downloaded_files)
    outpath = "../../genomeMSA_2021/clustalW/"
    for jobId in jobId_ls:
        for files in downloaded_files:
            
            if jobId in files:
                try:
                    print("----------- job found: %s ----------------"%jobId)
                    fin = open(inpath + "%s.aln-clustal_num.clustal_num"%jobId, 'r')
                    lines = fin.readlines()
                    line = lines[3]
                    name = line.split(' ')[0].split('_')
                    idx = name[-1]
                    date = name[1]
                    
                    path = outpath + "%s/"%date
                    createFolder(path)
                    filename = "GenomesGISAID_SARS-CoV-2_%s2021_MSA_%s.txt"%(date, idx)
                    fout = open( path + filename, "w")
                    fout.writelines(lines)
                    
                    fin.close()
                    fout.close()
                    
                    print("---------- finished transfering alignment to", path + filename, "------------------")
                    
                    print("---------- Deleting jobId and excess materia ----------l")
                    deletejobId(jobId) ; deleteResult(jobId)
                    break
                except:
                    print("Error", "%s.aln-clustal_num.clustal_num"%jobId)
            

def runMain():
    jobId_ls = obtain_all_jobId()
    download_finished_clustatlo(jobId_ls)
    obtain_all_finished(jobId_ls)
    
    
if __name__ == "__main__":
    runMain()