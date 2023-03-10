# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 09:40:48 2020

@author: Yuta
@email: hozumiyu@msu.edu
"""

import pandas as pd
import csv

'''
Remove bad dates, and store it into a csv
Just change the date for the current one to run this code
Input example:
    snpRecords_08072020_merged_unchecked.csv  
output is stored into 2 files
Example: 
    The valid dates: snpRecords_08072020_merged.csv  
    Invalid Dates: snpRecords_08072020_merged_problem.csv  

'''

def remove_bad_dates(infile, outfile):
    records = pd.read_csv(infile+'.csv')
    
    csv_columns = ["seqId", "snpsRef", "refNTs", "mutatedNTs", "mapNTs", "numNexts"]
    
    seqId = records["seqId"]
    
    thirty = []
    thirtyone = []
    twenty = []
    
    for i in range(1,32):
        if i < 31:
            if i < 10:
                thirty.append(str(i))
                thirty.append('0' + str(i))
            else:
                thirty.append(str(i))
        if i < 30:
            if i < 10:
                twenty.append(str(i))
                twenty.append('0' + str(i))
            else:
                twenty.append(str(i))
        
        if i < 10:
            thirtyone.append(str(i))
            thirtyone.append('0' + str(i))
        thirtyone.append(str(i))
        
    thirtyone_month = ['1','3', '5', '7', '01', '03', '05', '07', '8', '08', '10', '12']
    thirty_month = ['4', '6', '04', '06', '9', '09', '11']
    single = ['1', '2', '3', '4' ,'5', '6', '7', '8', '9']
    
    valid_year = ['2020', '2021']
    with open(outfile + '_problem.csv', 'w', newline='') as csvfile_problem:
        writerproblem = csv.DictWriter(csvfile_problem, fieldnames=csv_columns)
        writerproblem.writeheader()
        with open(outfile + '.csv', 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
                writer.writeheader()
                errorcount = 0
                for idx in range(len(seqId)):
                    
                    
                    validDate = True
                    full_name = seqId[idx]
                    splitname = full_name.split('|')
                    full_date = splitname[-1]
                    try:
                        year, month, date = full_date.split('-')
                    except:
                        validDate = False
                    
                    if validDate == True:
                        
                        #check that year is 2020 or 2021
                        if year not in valid_year:
                            validDate = False
                        
                        #check that the months and date is ok
                        if month in thirty_month:
                            if date not in thirty:
                                validDate = False
                        
                        elif month in thirtyone_month:
                            if date not in thirtyone:
                                validDate = False
                                
                        elif month in ['2', '02']:
                            if date not in twenty:
                                validDate = False
                        else: 
                            validDate = False
                        #check the 2021 
                        if year == "2021":
                            months2021 = ["01", "02", "03", "1", "2", "04", "4", "05", "5"]
                            if month not in months2021:
                                validDate = False
                        
                    if validDate == True:
                        data = {}
                        if month in single:
                            #adjust the dates if the date is single digit 4 -> 04
                            full_date = year + '-0' + month + '-' + date
                            name = ''
                            for i in range(len(splitname) - 1):
                                name = name + splitname[i]
                                name = name + '|'
                            name += full_date
                            print('Adjusting bad dates/month.')
                            print('Original:', records["seqId"][idx], 'Adjusted:', name)
                        else:
                            name = records["seqId"][idx]
                        #name = name.replace('|2021|', '|')
                        #name = name.replace('|2020|', '|')
                        data["seqId"] = name
                        data["snpsRef"] = records["snpsRef"][idx]
                        data["refNTs"] = records["refNTs"][idx]
                        data["mutatedNTs"] = records["mutatedNTs"][idx]
                        data["mapNTs"] = records["mapNTs"][idx]
                        data["numNexts"] = records["numNexts"][idx]
                        writer.writerow(data)
                    elif validDate == False:
                        errorcount += 1
                        print('idx=', idx, 'Problem with date:', full_name)
                        data = {}
                        data["seqId"] = records["seqId"][idx]
                        data["snpsRef"] = records["snpsRef"][idx]
                        data["refNTs"] = records["refNTs"][idx]
                        data["mutatedNTs"] = records["mutatedNTs"][idx]
                        data["mapNTs"] = records["mapNTs"][idx]
                        data["numNexts"] = records["numNexts"][idx]
                        writerproblem.writerow(data)
                        
                        
        print('Total of %s records removed' %str(errorcount))

if __name__ == "__main__":
    date = '0310'   #change the date
    infile = './%s/snpRecords_%s2021_merged_unchecked'%(date,date)
    outfile = './%s/snpRecords_%s2021_merged'%(date,date)
    remove_bad_dates(infile, outfile)
