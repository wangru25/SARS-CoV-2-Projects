import sys
import os


Date = sys.argv[1]
OldDate = sys.argv[2]

name_id = 'Protein_SNPs_%s.txt '%OldDate
for i in range(7):
    name_id += 'Protein_SNPs_%s_new_%d.txt '%(Date,i)
 
os.system('cat %s | sort |uniq > Protein_SNPs_%s.txt'%(name_id, Date))
#print('cat %s | sort |uniq > Protein_SNPs_%s.txt'%(name_id, Date))

#cat Protein_SNPs_04302021.txt Protein_SNPs_05242021_new_0.txt Protein_SNPs_05242021_new_1.txt Protein_SNPs_05242021_new_2.txt Protein_SNPs_05242021_new_3.txt Protein_SNPs_05242021_new_4.txt Protein_SNPs_05242021_new_5.txt Protein_SNPs_05242021_new_6.txt Protein_SNPs_05242021_new_7.txt Protein_SNPs_05242021_new_8.txt Protein_SNPs_05242021_new_9.txt Protein_SNPs_05242021_new_10.txt Protein_SNPs_05242021_new_11.txt Protein_SNPs_05242021_new_12.txt Protein_SNPs_05242021_new_13.txt Protein_SNPs_05242021_new_14.txt Protein_SNPs_05242021_new_15.txt Protein_SNPs_05242021_new_16.txt Protein_SNPs_05242021_new_17.txt Protein_SNPs_05242021_new_18.txt Protein_SNPs_05242021_new_19.txt Protein_SNPs_05242021_new_20.txt Protein_SNPs_05242021_new_21.txt | sort |uniq > Protein_SNPs_05242021.txt



