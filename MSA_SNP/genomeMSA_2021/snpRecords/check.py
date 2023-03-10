# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 11:21:36 2020

@author: Yutaho
"""


import pandas as pd
import csv 


def adjustCountry(country):
    if "Austria" in  country:
        country = "AT"
    elif "Brasil" in country:
        country = "BR"
    elif "Bealrus" in country:
        country = "BP"
    elif "Beijing" in country:
        country = "CN"
    elif "Hebei" in country:
        country = "CN"
    elif "Henan" in country:
        country = "CN"
    elif "CotedIvoire" in country:
        country = "CI"
    elif "Estonia" in country:
        country = "EE"
    elif "GuyaneFrancaise" in country:
        country = "FR"
    elif "France" in country:
        country = "FR"
    elif "Ireland" in country:
        country = "IE"
    elif "Libia" in country:
        country = "LY"
    elif "Liechtenstein" in country:
        country = "LI"
    elif "NorthernMarianaIslands" in country:
        country = "MP"
    elif "Curacao" in country:
        country = "NL"
    elif "Sweden" in country:
        country = "SE"
    elif "USA" in country:
        country = "US"
    elif "Germany" in country:
        country = "DE"
    elif "Switzerland" in country:
        country = "CH"
    elif "Liaoning" in country:
        country = "CN"
    elif "Italy" in country:
        country = "IT"
    elif "NorthAmerica" in country:
        country = "US"
    elif "Greece" in country:
        country = "GR"
    elif "Tahiti" in country:
        country = "FR"
    elif "CaymanIslands" in country:
        country = "UK"
        
    return country

def adjustLen5(splitId):
    country = splitId[0]
    date = splitId[-1]
    input_vec = splitId[1:-1]
    
    country = adjustCountry(country)
    
    for idx in range(len(input_vec)):
        if "EPI_" in input_vec[idx]:
            EPI_idx = idx
            break
    inputname = ''
    for idx in range(len(input_vec)):
        if idx != EPI_idx:
            inputname += input_vec[idx] + "-"
    inputname = inputname[:-1]
    seqId = country + "|" + inputname + "|" + input_vec[EPI_idx] + "|" + date
    return seqId
        


def adjustLen3(splitId):
    country, EPI, date = splitId
    country = adjustCountry(country)
    seqId = country + "|" + date[:4] + "|" + EPI + "|" + date
    return seqId

def adjustLen4(splitId):
    country ,input1, input2, date = splitId
    country = adjustCountry(country)
    if "EPI" in input1:
        seqId = country + "|" + input2 + "|" + input1 + "|" + date
    else:
        seqId = country + "|" + input1 + "|" + input2 + "|" + date
    return seqId


date = "0601-0607"

snpRecords = pd.read_csv("./MergedFiles/snpRecords_%s2021_new.csv"%( date))


seqId = snpRecords["seqId"]

csv_columns = ["seqId", "snpsRef", "refNTs", "mutatedNTs", "mapNTs", "numNexts"]
count = 0
with open("./MergedFiles/snpRecords_%s2021_new2.csv"%( date), "w", newline = '') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
    writer.writeheader()
    for idx, s in enumerate(seqId):
        data = {}
        s = s.replace('env|', '')
        s = s.replace('||', '|')
        s = s.replace('hCov-19|', '')
        s = s.replace('hCoV19|', '')
        
        
        
        if s[:3] == "cat":
            print(s, "removed")
        elif s[:3] == "dog":
            print(s, "removed")
        elif s[:5] == "mouse":
            print(s, "removed")
        elif s[:7] == "gorilla":
            print(s, "removed")
        elif s[:5] == "tiger":
            print(s, "removed")
        elif s[:7] == "leopard":
            print(s, "removed")
        elif s[:4] == "lion":
            print(s, "removed")
            
        else:
            splitId = s.split('|')
            if len(splitId) == 3:
                name = adjustLen3(splitId)
            elif len(splitId) == 4:
                name = adjustLen4(splitId)
            elif len(splitId) >= 5:
                name = adjustLen5(splitId)
            else:
                print(idx, "Length of name too long")
            
            data["seqId"] = name
            data["snpsRef"] = snpRecords["snpsRef"][idx]
            data["refNTs"] = snpRecords["refNTs"][idx]
            data["mutatedNTs"] = snpRecords["mutatedNTs"][idx]
            data["mapNTs"] = snpRecords["mapNTs"][idx]
            data["numNexts"] = snpRecords["numNexts"][idx]
            writer.writerow(data)
            count += 1
    

    
snpRecords2 = pd.read_csv("./MergedFiles/snpRecords_%s2021_new2.csv"%(date))


seqId2 = snpRecords2["seqId"]
country = []
for idx3, s in enumerate(seqId2):
    s = s.split('|')
    if len(s) != 4:
        print(idx3+2, "incorrect length")
        print(s)
    d = s[-1]
    if len(d) != 10:
        print(idx3+2, "incorrect date")
        print(s)
    if len(s[0]) != 2:
        print(idx3+2, "incorrect country length")
        print(s)
    country.append(s[0])

snpsRef = snpRecords2["snpsRef"]
for idx4, snps in enumerate(snpsRef):
    if len(eval(snps)) > 200:
        print("too many snp in idx", idx4+2)
    
country = list(set(country))
