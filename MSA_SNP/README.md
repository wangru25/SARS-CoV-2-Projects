# How to Decode Mutations

Preprocessing genome information has two primary steps.

#### Step 1: Perform multiple sequence alignment

- Fasta file is stored in 2019-nCov/data/GISAID_2021/date/
  	date: submission date
  	The file should be name gisaid_hcov-19_2021_date.fasta  
  		with the appropriate date (ex: 0603A, 0511) *no space in between
  	You should also keep the .tsv file: gisaid_hcov-19_2021_date.TSV in the same folder (this is not used in our analysis but useful for future citation)

- In 2019-nCov/data/GISAID_2021/, run 

  ```python
  python mergeFasta.py date1,date2,..
  ```

- And 

  ```python
  python mergeFasta.py
  ```

  splits data for each date into a group of 100 sequence each. For example: mergeFasta.py 0401,0402 will split data for dates 0401, 0402 into 100 sequence each
  **IMPORTANT**: NC_045512.fasta must be in 2019-nCov/data/GISAID_2021 -- This is the reference genome

	This program will give you an output 2019-nCov/data/GISAID_2021/country2.txt
	country2.txt has all the mapping from country name to its ISO-2-letter code. Please make sure that the mapping is correct
	Sometimes, cities are listed instead of the country. You can add exception in def countryCode(self, countryT).
	Each date will have mapping informations, such as 1_mapIndex_0603B.csv. You will need this later, so do not delete

- Next, run

```
python runClustalo.py
```

 in 2019-nCov/data/GISAID_2021/. This sends data to clusal omega's server. Only the unaligned data will be sent. This takes quite some time. The jobid number will be saved in 2019-nCov/data/GISAID_2021/jobId

- Finally, run

  ```
  python download_clustalo_result.py 
  ```

  This downloads the result. When you send alot of data, this could take somewhere between 3 hours to 3 days. The downloaded data will be moved to 2019-nCoV/genomeMSA_2021/clustalW/date

#### Step 2: Find the mutation mapping of the aligned sequence.

- In 2019-nCoV/genomeMSA_2021/, run 

```
python genomeMain.py
```

​	This processes mutation mapping for all the downloaded files saved in 2019-nCoV/genomeMSA_2021/clustalW/date/
​	This will use the mapping data from 2019-nCov/data/GISAID_2021/date
​	The aligned data will be saved in 2019-nCoV/genomeMSA_2021/snpRecords/date/

- In 2019-nCoV/genomeMSA_2021/snpRecords/ run 

  ```
  python mergeFiles.py subMonth subMinDate subMaxDate
  ```

  ​	For example, if you want to merge dates 05/01-05/15. python mergeFiles.py 05 01 15
  ​	Merged records will be stored in 2019-nCoV/genomeMSA_2021/snpRecords/MergedFiles/snpRecords_0501-0515_new.csv There are other folders: badRecords: store all the bad records (such as improper dates)

  ​                *badMutation: if the data is corrupt, remove the data   *most likely not needed

  ​                *badMutation2: if the data is corrupt and it is causing problem with the entire alignment data *most likely not needed

- In 2019-nCoV/genomeMSA_2021/snpRecords/, run 

  ```
  python check.py 
  ```

  ​	*Open this in a compiler and change the subMonth, subMinDate, subMaxDate.
  ​	This program double checks to make sure seqId is correct. 
​	seqId format should be ISO-2-letter-countrycode|extra stuff|AccessionID (EPI_xxx)|date (year-month-date)
  ​	*Date and country code is the most important.
​	Sometimes, animal data will be mixed. Please remove these files.

​	



