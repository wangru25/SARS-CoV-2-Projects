# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 11:21:36 2020

@author: Yutaho
"""


import pandas as pd
import csv 

date = "0329"

snpRecords = pd.read_csv("snpRecords/%s/snpRecords_%s2021_new.csv"%(date, date))

a = [4,5]
seqId = snpRecords["seqId"]

csv_columns = ["seqId", "snpsRef", "refNTs", "mutatedNTs", "mapNTs", "numNexts"]
with open("snpRecords/%s/snpRecords_%s2021_new2.csv"%(date, date), "w", newline = '') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
    writer.writeheader()
    for idx, s in enumerate(seqId):
        data = {}
        s = s.replace('env|', '')
        s = s.replace('||', '|')
        s = s.split('|')
        name = ''
        for ss in s:
            if ss == "2021":
                continue
            elif ss == "2020":
                continue
            elif ss == "hCoV19":
                continue
            elif ss == "France":
                name = name + "FR" + "|"
            elif ss == "Hebei" or ss == "Henan":
                name = name + "CN" + "|"
            elif ss == "Brasil":
                name = name + "BR" + "|"
            elif ss == "GuyaneFrancaise":
                name = name + "FR" + "|"
            elif ss == "Austria":
                name = name + "AT" + "|"
            elif ss == "Curacao":
                name = name + "NL" + "|"
            elif ss == "USA":
                name = name + "US" + "|"
            else:
                name = name + ss + "|"
        name = name[:-1]
        
        data["seqId"] = name
        data["snpsRef"] = snpRecords["snpsRef"][idx]
        data["refNTs"] = snpRecords["refNTs"][idx]
        data["mutatedNTs"] = snpRecords["mutatedNTs"][idx]
        data["mapNTs"] = snpRecords["mapNTs"][idx]
        data["numNexts"] = snpRecords["numNexts"][idx]
        writer.writerow(data)
    

    
snpRecords = pd.read_csv("snpRecords/%s/snpRecords_%s2021_new2.csv"%(date, date))

a = [4,5]
seqId = snpRecords["seqId"]
country = []
for idx, s in enumerate(seqId):
    s = s.split('|')
    if len(s) != 4:
        print(idx+2, "incorrect length")
    d = s[-1]
    if len(d) != 10:
        print(idx+2, "incorrect date")
    if len(s[0]) != 2:
        print(idx+2, "incorrect length")
    country.append(s[0])
    
country = list(set(country))
