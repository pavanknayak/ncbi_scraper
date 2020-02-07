import requests, bs4
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm



#This script converts all ensembl IDs into conventional gene names
#input is the name of the excel file containing the ensembl ID list

def ensembl_id_finder(idlist):

#read the excel file as a dataframe and convert to list
    ens_df = pd.read_excel(idlist)
    ens_dfli = ens_df.values.tolist()
#preallocate list
    ens_dfli_offgene = [0] * len(ens_dfli)

#loop through the list of ensembl IDs
    for k in tqdm(range(0, len(ens_dfli))):

#access ensembl API and look up each ensembl ID on the server
        server = "http://rest.ensembl.org"
        ext = "/lookup/id/" + ens_dfli[k][0]

#get the gene info for the ensembl ID
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
         
    #    if not r.ok:
    #      r.raise_for_status()
    #      sys.exit()
    #

#convert to json file
        decoded = r.json()
#if the ensembl ID is not found or is outdated, print error next to the respective ensembl ID
        if 'error' in decoded:
            ens_gene_name = 'ERROR:NO ENSEMBL ID FOUND'
            ens_dfli_offgene[k] = ens_gene_name
#else, put the correct gene name in the respective spot in the conventional gene list
        else:
            ens_gene_name = decoded['display_name']
            ens_dfli_offgene[k] = ens_gene_name
        

#append the conventional gene list to the original ensembl ID dataframe
    ens_df.insert(len(ens_df.columns), 'Gene Name', ens_dfli_offgene, True)

#output as an excel file
    return ens_df.to_excel('ensembl-id-output1.xlsx')


#This script takes a gene list plus a keyword of your choice and searches each gene plus your keyword into pubmed and returns the number
# of articles published for that particular search term

#The input gene list file has to have a header in the first row (Name) and a list of genes


#import packages and read excel file of genes, then convert to list

#The inputs  *inputfile* which is a string of the excel file with the gene list
#the input *key* is a string of keywords to search on pubmed for each gene
def pubmedscrape( inputfile, key):

    file = inputfile
    df = pd.read_excel(file)
    dfli = df.values.tolist()
    
    #if the list is an ensembl ID list, first find the conventional names of all the genes
    ens_id_str = dfli[0][0]
    ens_id_str1 = dfli[1][0]
    ens_id_str2 = dfli[2][0]
    if ens_id_str.startswith('ENS') & ens_id_str1.startswith('ENS') & ens_id_str2.startswith('ENS'):
        
        df = pd.read_excel('ensembl-id-output.xlsx')
        dfli = df.values.tolist()
        #preallocate a list of counts with zeros

        dfli_count = [0] * len(dfli)

        #loop through gene list
        for i in range(0, len(dfli)):

        #with each gene plus your keywords
            gene = dfli[i][2]
            
            if 'ERROR' in gene:
                continue
            keyword = key

        #remove all spaces before and after the gene and keywords
            gene = gene.strip()
            keyword = keyword.strip()
        #replace the spaces in the keyword with '+' so that they fit in the ncbi search URL format
            keyword = keyword.replace(' ', '+')

        #get the html page for your search term and convert to text
            res = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=' + gene + '+'+ keyword)

            ncbi = bs4.BeautifulSoup(res.text)

        #find all the HTML elements with <count></count> and take the first one (this is the one with the actual number of papers)
            count = ncbi.select('Count')
            if count == []:
                dfli_count[i] = 0
                continue
            once = str(count[0])
        #remove the <count></count> tags
            once = once.replace('<count>', '')
            once = once.replace('</count>', '')
        #add to the respective place in list of counts    
            dfli_count[i] = once

        #Change all str to numbers in the count list
        for j in range(0, len(dfli_count)):
            dfli_count[j] = int(dfli_count[j])


        #append the list of counts column to the end of the original gene list column
        df.insert(len(df.columns), 'Pubmed Count', dfli_count, True)
        df.drop(df.columns[[0]], axis=1, inplace=True)
        #output to an excel doc named gene-list-output + the keyword searched
        searchterm = keyword.replace(' ', '-')
        searchterm = keyword.strip()
        output = df.to_excel('gene-list-output' + '-' + searchterm + '.xlsx')
        
        return output
        
    else:
    
        #preallocate a list of counts with zeros

        dfli_count = [0] * len(dfli)

        #loop through gene list
        for i in range(0, len(dfli)):

        #with each gene plus your keywords
            gene = dfli[i][0]
            keyword = key

        #remove all spaces before and after the gene and keywords
            gene = gene.strip()
            keyword = keyword.strip()
        #replace the spaces in the keyword with '+' so that they fit in the ncbi search URL format
            keyword = keyword.replace(' ', '+')

        #get the html page for your search term and convert to text
            res = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=' + gene + '+'+ keyword)

            ncbi = bs4.BeautifulSoup(res.text)

        #find all the HTML elements with <count></count> and take the first one (this is the one with the actual number of papers)
            count = ncbi.select('Count')
            if count == []:
                dfli_count[i] = 0
                continue
            once = str(count[0])
        #remove the <count></count> tags
            once = once.replace('<count>', '')
            once = once.replace('</count>', '')
        #add to the respective place in list of counts    
            dfli_count[i] = once

        #Change all str to numbers in the count list
        for j in range(0, len(dfli_count)):
            dfli_count[j] = int(dfli_count[j])


        #append the list of counts column to the end of the original gene list column
        df.insert(len(df.columns), 'Pubmed Count', dfli_count, True)
        
        #output to an excel doc named gene-list-output + the keyword searched
        searchterm = keyword.replace(' ', '-')
        searchterm = keyword.strip()
        output = df.to_excel('gene-list-output' + '-' + searchterm + '.xlsx')
        
        return output
    


