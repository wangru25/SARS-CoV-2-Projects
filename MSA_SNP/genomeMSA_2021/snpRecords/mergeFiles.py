# -*- coding: utf-8 -*-
'''
@Author: Rui Wang
@Date: 2020-05-02 21:41:38
@LastModifiedBy: Rui Wang
@LastEditTime: 2020-05-03 02:10:28
@Email: wangru25@msu.edu
@FilePath: /38_Influ/genomeMSA/snpRecords/mergeFiles.py
@Description: 
'''
import pandas as pd 
import numpy as np
import os
import sys
import removedate as rd
from os import listdir
from os.path import isfile, join

#change this if needed
mutation_threshold = 200

month2021 = ['1', '01', '2', '02', '3', '03', '4', '04', '5', '05', '6', '06']
        

def checkDate(seqId, month2021 = month2021):
    '''
        Checks date to make sure the date is valid
        month2021: shows the month we are 
        output: True/False, seqId
            True if the date is good
            False if the date is bad
            seqId: empty if bad, adjusted if good

    '''
    thirty = ['1', '01', '2', '02', '3', '03', '4', '04', '5', '05', '6', '06', '7', '07', '8', '08', '9', '09', 
          '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', 
          '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30']
    thirtyone = ['1', '01', '2', '02', '3', '03', '4', '04', '5', '05', '6', '06', '7', '07', '8', '08', '9', '09',
                 '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', 
                 '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31']
    twenty = ['1', '01', '2', '02', '3', '03', '4', '04', '5', '05', '6', '06', '7', '07', '8', '08', '9', '09', 
              '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', 
              '20', '21', '22', '23', '24', '25', '26', '27', '28', '29']
    thirtyone_month = ['1','3', '5', '7', '01', '03', '05', '07', '8', '08', '10', '12']
    thirty_month = ['4', '6', '04', '06', '9', '09', '11']
    single = ['1', '2', '3', '4' ,'5', '6', '7', '8', '9']
    
    valid_year = ['2020', '2021']
    
    #extract the date
    splitName = seqId.split('|')
    fulldate = splitName[-1]
    try:
        year, month, date = fulldate.split('-')
    except:
        return False, seqId
    
    #check year
    if year not in valid_year:
        return False, seqId
    
    #check for month-date combo
    if month in thirty_month:
        if date not in thirty:
            return False, seqId
    elif month in thirtyone_month:
        if date not in thirtyone:
            return False, seqId
    elif month in ['2', '02']:
        if date not in twenty:
            return False, seqId
    else: 
        return False, seqId
    
    #if the year is 2021, make sure that the data is not in the future
    if year == "2021":
        if month not in month2021:
            return False, seqId
    
    #fix the date
    fulldate_fixed = year
    if month in single:
        fulldate_fixed += '-0' + month
    else:
        fulldate_fixed += '-' + month
    
    if date in single:
        fulldate_fixed += '-0' + date
    else:
        fulldate_fixed += '-' + date
    
    seqIdFixed = ''
    for i in range(len(splitName) - 1):
        seqIdFixed += splitName[i]
        seqIdFixed += '|'
    seqIdFixed += fulldate_fixed
    
    return True, seqIdFixed
    

def getAllFolder(subMonth, subMinDate, subMaxDate):
    '''
        Obtain all the folders associated with the date
        ex: [0410A, 0410B, 0411]
    '''
    allDate = list(range(int(subMinDate) , int(subMaxDate) + 1)) 
    
    allFolders = [x[0] for x in os.walk("./")]
    
    folder_ls = []
    for date in allDate:
        if date < 10:
            date_temp = subMonth + str(0) + str(date)  #full date
        else:
            date_temp = subMonth + str(date)  #full date
        #get all the folder associated with the date
        for folder in allFolders:
            if date_temp in folder:
                folder_ls.append(folder) 
    return folder_ls


def checkOneDate(folder, mutation_threshold = mutation_threshold):
    '''
        Check each file in te folder to make sure that the mutation makes sense.
        If there are too many mutation (more than the threshold), add to badSeqId and badFile.
        badSeqId: The SeqId of the bad file
        badFile: the file nam. May need to rerun clustal for this one

    '''
    inpath = folder[2:]  + '/'
    allfiles = [f for f in listdir(inpath) if isfile(join(inpath, f))]
    
    badSeqId = []
    badFile = []
    for file in allfiles:
        if '.txt' not in file:
            if "merged" not in file:
                if "bad" not in file:
                    snpRecords = pd.read_csv(inpath + file)
                    seqId_ls = snpRecords["seqId"]
                    snpsRef_ls = snpRecords["snpsRef"]
                    total = len(snpsRef_ls)
                    for idx in range(total):
                        snpsRef = snpsRef_ls[idx]
                        if len(eval(snpsRef)) > mutation_threshold:
                            badSeqId.append(seqId_ls[idx])
                            badFile.append(file)
                
    return badSeqId, badFile

def mergeOneData(folder, badSeqId, badFile):
    '''
        Merge the files in one folder.
        Removes any record with too many mutation, which may need to be realigned
    '''
    
    inpath = folder[2:]  + '/'
    allfiles = [f for f in listdir(inpath) if isfile(join(inpath, f))]
    
    total = 0   #count the number of relavant files
    for f in allfiles:
        if "merge" not in f:
            if "new" not in f:
                if "bad" not in f:
                    total +=1
    
    fileNumbers = list(range(1, total + 1))   #change the second number (always the max + 1)
    frames = []
    for i in fileNumbers:
        frames.append(pd.read_csv('./%s/snpRecords_%s2021_%s.csv'%(folder[2:], folder[2:],str(i))))
        
    concactenated_record = {"seqId": [], "snps": [], "snpsRef": [], "refNTs":[], "mutatedNTs": [], "mapNTs": [] , "numNexts": []}
    bad_concactenated_record = {"seqId": [], "snps": [], "snpsRef": [], "refNTs":[], "mutatedNTs": [], "mapNTs": [] , "numNexts": []}
    bad_date_record = {"seqId": [], "snps": [], "snpsRef": [], "refNTs":[], "mutatedNTs": [], "mapNTs": [] , "numNexts": []}
    concactenate_temp = pd.concat(frames, ignore_index = True)
    
    seqId_temp = concactenate_temp["seqId"]
    print(badSeqId)
    totalId = len(seqId_temp)
    for idx in range(totalId):
        
        if seqId_temp[idx] in badSeqId:
            print("Too many mutation in ", seqId_temp[idx])
            print(concactenate_temp["snpsRef"][idx])
            bad_concactenated_record["seqId"].append(concactenate_temp["seqId"][idx])
            bad_concactenated_record["snps"].append(concactenate_temp["snps"][idx])
            bad_concactenated_record["snpsRef"].append(concactenate_temp["snpsRef"][idx])
            bad_concactenated_record["refNTs"].append(concactenate_temp["refNTs"][idx])
            bad_concactenated_record["mutatedNTs"].append(concactenate_temp["mutatedNTs"][idx])
            bad_concactenated_record["mapNTs"].append(concactenate_temp["mapNTs"][idx])
            bad_concactenated_record["numNexts"].append(concactenate_temp["numNexts"][idx])
        else:
            valid, seqId = checkDate(seqId_temp[idx])
            if valid:
                concactenated_record["seqId"].append(seqId)
                concactenated_record["snps"].append(concactenate_temp["snps"][idx])
                concactenated_record["snpsRef"].append(concactenate_temp["snpsRef"][idx])
                concactenated_record["refNTs"].append(concactenate_temp["refNTs"][idx])
                concactenated_record["mutatedNTs"].append(concactenate_temp["mutatedNTs"][idx])
                concactenated_record["mapNTs"].append(concactenate_temp["mapNTs"][idx])
                concactenated_record["numNexts"].append(concactenate_temp["numNexts"][idx])
            else:
                bad_date_record["seqId"].append(seqId)
                bad_date_record["snps"].append(concactenate_temp["snps"][idx])
                bad_date_record["snpsRef"].append(concactenate_temp["snpsRef"][idx])
                bad_date_record["refNTs"].append(concactenate_temp["refNTs"][idx])
                bad_date_record["mutatedNTs"].append(concactenate_temp["mutatedNTs"][idx])
                bad_date_record["mapNTs"].append(concactenate_temp["mapNTs"][idx])
                bad_date_record["numNexts"].append(concactenate_temp["numNexts"][idx])
            
    concactenated_record = pd.DataFrame(concactenated_record, columns=concactenated_record.keys())
    concactenated_record.to_csv('./%s/snpRecords_%s2021_merged.csv' %(folder[2:],folder[2:]), index=False)
    
    bad_date_record = pd.DataFrame(bad_date_record, columns=bad_date_record.keys())
    bad_date_record.to_csv('./%s/snpRecords_%s2021_merged_badDates.csv' %(folder[2:],folder[2:]), index=False)
    
    bad_concactenated_record = pd.DataFrame(bad_concactenated_record, columns=bad_concactenated_record.keys())
    bad_concactenated_record.to_csv('./%s/snpRecords_%s2021_badMutation.csv' %(folder[2:],folder[2:]), index=False)
    
    
def mergeOneData_v2(folder, badFile):
    '''
        Merge the files in one folder. 
        Removes any file with too many mutation, which may need to be realigned
    '''
    badFile = list(set(badFile))
    
    inpath = folder[2:]  + '/'
    allfiles = [f for f in listdir(inpath) if isfile(join(inpath, f))]
    
    total = 0   #count the number of relavant files
    for f in allfiles:
        if "merge" not in f:
            if "new" not in f:
                if "bad" not in f:
                    total +=1
    
    fileNumbers = list(range(1, total + 1))   #change the second number (always the max + 1)
    frames = []
    badframes = []
    for i in fileNumbers:
        currentFile = 'snpRecords_%s2021_%s.csv'
        if currentFile not in badFile:
            frames.append(pd.read_csv('./%s/snpRecords_%s2021_%s.csv'%(folder[2:], folder[2:],str(i))))
        else:
            badframes.append(currentFile)
        
    concactenated_record = {"seqId": [], "snps": [], "snpsRef": [], "refNTs":[], "mutatedNTs": [], "mapNTs": [] , "numNexts": []}
    bad_date_record = {"seqId": [], "snps": [], "snpsRef": [], "refNTs":[], "mutatedNTs": [], "mapNTs": [] , "numNexts": []}
    concactenate_temp = pd.concat(frames, ignore_index = True)
    
    seqId_temp = concactenate_temp["seqId"]
    totalId = len(seqId_temp)
    for idx in range(totalId):
        valid, seqId = checkDate(seqId_temp[idx])
        if valid:
            concactenated_record["seqId"].append(seqId)
            concactenated_record["snps"].append(concactenate_temp["snps"][idx])
            concactenated_record["snpsRef"].append(concactenate_temp["snpsRef"][idx])
            concactenated_record["refNTs"].append(concactenate_temp["refNTs"][idx])
            concactenated_record["mutatedNTs"].append(concactenate_temp["mutatedNTs"][idx])
            concactenated_record["mapNTs"].append(concactenate_temp["mapNTs"][idx])
            concactenated_record["numNexts"].append(concactenate_temp["numNexts"][idx])
        else:
            bad_date_record["seqId"].append(seqId)
            bad_date_record["snps"].append(concactenate_temp["snps"][idx])
            bad_date_record["snpsRef"].append(concactenate_temp["snpsRef"][idx])
            bad_date_record["refNTs"].append(concactenate_temp["refNTs"][idx])
            bad_date_record["mutatedNTs"].append(concactenate_temp["mutatedNTs"][idx])
            bad_date_record["mapNTs"].append(concactenate_temp["mapNTs"][idx])
            bad_date_record["numNexts"].append(concactenate_temp["numNexts"][idx])
            
    concactenated_record = pd.DataFrame(concactenated_record, columns=concactenated_record.keys())
    concactenated_record.to_csv('./%s/snpRecords_%s2021_merged_v2.csv' %(folder[2:],folder[2:]), index=False)
    
    bad_date_record = pd.DataFrame(bad_date_record, columns=bad_date_record.keys())
    bad_date_record.to_csv('./%s/snpRecords_%s2021_badDates_v2.csv' %(folder[2:],folder[2:]), index=False)
    file = open("./%s/badFiles_%s2021.txt"%(folder[2:],folder[2:]), "w")
    for line in badframes:
        file.writelines(line + "\n")
    print("bad mutation file number:", len(badframes))
    file.close()
    
    
def mergeAllGood(folder_ls, subMonth, subMinDate, subMaxDate):
    
    frames = []
    for date in folder_ls:
        frames.append(pd.read_csv('./%s/snpRecords_%s2021_merged.csv' %(date[2:],date[2:])))
    concactenatedFrames = pd.concat(frames, ignore_index = True)
    concactenatedFrames = pd.DataFrame(concactenatedFrames, columns=concactenatedFrames.keys())
    concactenatedFrames.to_csv('./MergedFiles/snpRecords_%s-%s2021_new.csv' %(subMonth+subMinDate, subMonth+subMaxDate), index=False)
    
    print("Number of good records:", len(concactenatedFrames) )
    return

def mergeBadDate(folder_ls, subMonth, subMinDate, subMaxDate):
    frames = []
    for date in folder_ls:
        frames.append(pd.read_csv('./%s/snpRecords_%s2021_merged_badDates.csv' %(date[2:],date[2:])))
    concactenatedFrames = pd.concat(frames, ignore_index = True)
    concactenatedFrames = pd.DataFrame(concactenatedFrames, columns=concactenatedFrames.keys())
    concactenatedFrames.to_csv('./MergedFiles_badDate/snpRecords_%s-%s2021_badDates.csv' %(subMonth+subMinDate, subMonth+subMaxDate), index=False)
    
    print("Number of badDates:", len(concactenatedFrames))
    
    
def mergeBadMutation(folder_ls, subMonth, subMinDate, subMaxDate):
    #MergedFiles_badMutation
    frames = []
    for date in folder_ls:
        frames.append(pd.read_csv('./%s/snpRecords_%s2021_badMutation.csv' %(date[2:],date[2:])))
    concactenatedFrames = pd.concat(frames, ignore_index = True)
    concactenatedFrames = pd.DataFrame(concactenatedFrames, columns=concactenatedFrames.keys())
    concactenatedFrames.to_csv('./MergedFiles_badMutation/snpRecords_%s-%s2021_badMutation.csv' %(subMonth+subMinDate, subMonth+subMaxDate), index=False)
    print("Number of bad records:", len(concactenatedFrames))
    
    
def mergeAllGood_v2(folder_ls, subMonth, subMinDate, subMaxDate):
    
    frames = []
    for date in folder_ls:
        frames.append(pd.read_csv('./%s/snpRecords_%s2021_merged_v2.csv' %(date[2:],date[2:])))
    concactenatedFrames = pd.concat(frames, ignore_index = True)
    concactenatedFrames = pd.DataFrame(concactenatedFrames, columns=concactenatedFrames.keys())
    concactenatedFrames.to_csv('./MergedFiles_v2/snpRecords_%s-%s2021_new_v2.csv' %(subMonth+subMinDate, subMonth+subMaxDate), index=False)
    
    print("Number of good records:", len(concactenatedFrames) )
    return

def mergeBadDate_v2(folder_ls, subMonth, subMinDate, subMaxDate):
    frames = []
    for date in folder_ls:
        frames.append(pd.read_csv('./%s/snpRecords_%s2021_badDates_v2.csv' %(date[2:],date[2:])))
    concactenatedFrames = pd.concat(frames, ignore_index = True)
    concactenatedFrames = pd.DataFrame(concactenatedFrames, columns=concactenatedFrames.keys())
    concactenatedFrames.to_csv('./MergedFiles_badDate_v2/snpRecords_%s-%s2021_badDates_v2.csv' %(subMonth+subMinDate, subMonth+subMaxDate), index=False)
    
    print("Number of badDates:", len(concactenatedFrames))
    
    
def mergeBadMutation_v2(folder_ls, subMonth, subMinDate, subMaxDate):
    #MergedFiles_badMutation
    frames = []
    for date in folder_ls:
        file = open(("./%s/badFiles_%s2021.txt"%(date[2:],date[2:])))
        lines = file.readlines()
        frames.append(lines)
        file.close()
    
    file = open("./MergedFiles_badMutation_v2/badMutation_%s-%s_2021.txt"%(subMonth+subMinDate, subMonth+subMaxDate), "w")
    for line in lines:
        file.writelines(lines + "\n")
    file.close()


if __name__ == "__main__":
    subMonth = "06"
    subMinDate = '01'
    subMaxDate = '07'
    folder_ls = getAllFolder(subMonth, subMinDate, subMaxDate)
    for folder in folder_ls:
        print("Merging:", folder[2:])
        badSeqId, badFile = checkOneDate(folder)
        mergeOneData(folder, badSeqId, badFile)
    print("Merging all Folders")
    mergeAllGood(folder_ls, subMonth, subMinDate, subMaxDate)
    mergeBadDate(folder_ls, subMonth, subMinDate, subMaxDate)
    mergeBadMutation(folder_ls, subMonth, subMinDate, subMaxDate)
    
    for folder in folder_ls:
        print("Merging:", folder[2:])
        badSeqId, badFile = checkOneDate(folder)
        mergeOneData_v2(folder, badFile)
    mergeAllGood_v2(folder_ls, subMonth, subMinDate, subMaxDate)
    mergeBadDate_v2(folder_ls, subMonth, subMinDate, subMaxDate)
    mergeBadMutation_v2(folder_ls, subMonth, subMinDate, subMaxDate)
    
'''
if __name__ == "__main__":
    #subMonth = sys.argv[1]
    #subMaxMonth = int(sys.argv[2])
    subMonth = "05"
    subMinDate = 15
    subMaxDate = 31
    date_num = list(range(subMinDate, subMaxDate + 1))
    date_vec = []
    
            
    numFolder = 0
    numFile = 0
    #get all the folder names
    allFolders = [x[0] for x in os.walk("./")]
    
    
    #iterate over all the dates
    for date in date_vec:
        date_temp = subMonth + date  #full date
        folder_ls = []
        #get all the folder associated with the date
        for folder in allFolders:
            if date_temp in folder:
                folder_ls.append(folder)
        
        #iterate through the folder with the same date ie 0301US, 0301EU
        for folder in folder_ls:
            numFolder += 1
            newDate = folder[2:]
            inpath = folder + '/'
            
            #get all the files in the folder
            allfiles = [f for f in listdir(inpath) if isfile(join(inpath, f))]
            
            total = 0   #count the number of relavant files
            for f in allfiles:
                if "merge" not in f:
                    if "new" not in f:
                        total +=1
            number = list(range(1, total + 1))   #change the second number (always the max + 1)
            print(total)
            frames = []
            if number != []:
                for i in number:
                    numFile += 1
                    frames.append(pd.read_csv('./%s/snpRecords_%s2021_%s.csv'%(newDate, newDate,str(i))))
                    print('snpRecords_%s2021_%s.csv'%(newDate,str(i)))
                
                newsnp = pd.concat(frames)
                newsnp.to_csv('./%s/snpRecords_%s2021_merged_unchecked.csv' %(newDate,newDate), index=False)
                a = pd.read_csv('./%s/snpRecords_%s2021_merged_unchecked.csv'%(newDate,newDate))
                newdf = a.drop_duplicates(subset=None, keep='first')
                newdf.to_csv('./%s/snpRecords_%s2021_merged_unchecked.csv'%(newDate,newDate), index=False)
                print("finished merging %s"%newDate)
                print("removing bad dates")
                infile = './%s/snpRecords_%s2021_merged_unchecked'%(newDate,newDate)
                outfile = './%s/snpRecords_%s2021_merged'%(newDate,newDate)
                rd.remove_bad_dates(infile, outfile)
            
            
            
            
            if oldDate == None:
                records = pd.read_csv('./%s/snpRecords_%s2021_merged.csv'%(newDate,newDate))
                records.to_csv('./%s/snpRecords_%s2021_new.csv'%(newDate,newDate), index=False)
                oldDate = newDate
                records_p = pd.read_csv('./%s/snpRecords_%s2021_merged_problem.csv'%(newDate,newDate))
                records_p.to_csv('./%s/snpRecords_%s2021_new_problem.csv'%(newDate,newDate), index=False)
                
                
            else:
                frames = [] 
                frames.append(pd.read_csv('./%s/snpRecords_%s2021_new.csv'%(oldDate,oldDate)))
                frames.append(pd.read_csv('./%s/snpRecords_%s2021_merged.csv'%(newDate,newDate)))
                
                
                newsnp = pd.concat(frames)
                newsnp.to_csv('./%s/snpRecords_%s2021_new.csv'%(newDate,newDate), index=False)
                a = pd.read_csv('./%s/snpRecords_%s2021_new.csv'%(newDate,newDate))
                newdf = a.drop_duplicates(subset=None, keep='first')
                newdf.to_csv('./%s/snpRecords_%s2021_new.csv'%(newDate,newDate), index=False)
                
                frames = []
                frames.append(pd.read_csv('./%s/snpRecords_%s2021_new_problem.csv'%(oldDate,oldDate)))
                frames.append(pd.read_csv('./%s/snpRecords_%s2021_merged_problem.csv'%(newDate,newDate)))
                
                
                newsnp = pd.concat(frames)
                newsnp.to_csv('./%s/snpRecords_%s2021_new_problem.csv'%(newDate,newDate), index=False)
                a = pd.read_csv('./%s/snpRecords_%s2021_new_problem.csv'%(newDate,newDate))
                newdf = a.drop_duplicates(subset=None, keep='first')
                newdf.to_csv('./%s/snpRecords_%s2021_new_problem.csv'%(newDate,newDate), index=False)
                
                print("finished processing")
                
                oldDate = newDate
    print("total number of dates merged:", numFolder)
    print("total number of file:", numFile)
'''