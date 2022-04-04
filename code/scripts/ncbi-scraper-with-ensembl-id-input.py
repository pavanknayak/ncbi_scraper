import requests, bs4
import sys
import numpy as np
import os
import pandas as pd
from tqdm import tqdm
import argparse


#This script converts all ensembl IDs into conventional gene names
#input is the name of the excel file containing the ensembl ID list

def ensembl_id_finder(idlist):

#read the excel file as a dataframe and convert to list
    ens_df = pd.read_excel(idlist)
    ens_dfli = ens_df.values.tolist()
#preallocate list
    ens_dfli_offgene = [0] * len(ens_dfli)

#loop through the list of ensembl IDs
    print('CONVERTING ENSEMBL IDs TO CONVENTIONAL GENE NAMES...')
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
    return ens_df.to_excel('ensembl-id-temp.xlsx')


#This script takes a gene list plus a keyword of your choice and searches each gene plus your keyword into pubmed and returns the number
# of articles published for that particular search term

#The input gene list file has to have a header in the first row (Name) and a list of genes


#import packages and read excel file of genes, then convert to list

#The inputs  *inputfile* which is a string of the excel file with the gene list
#the input *key* is a string of keywords to search on pubmed for each gene
def pubmedscrape( inputfile, key, outname):

    file = inputfile
    df = pd.read_excel(file)
    if 'Unnamed: 0' in df.columns:
        df.drop('Unnamed: 0', axis=1)
        dfli = df.values.tolist()



        #if the list is an ensembl ID list, first find the conventional names of all the genes
        ens_id_str = dfli[0][1]
        ens_id_str1 = dfli[1][1]
        ens_id_str2 = dfli[2][1]
        if ens_id_str.startswith('ENS') & ens_id_str1.startswith('ENS') & ens_id_str2.startswith('ENS'):
            dfli = df.values.tolist()
            #preallocate a list of counts with zeros

            dfli_count = [0] * len(dfli)
            #loop through gene list
            print('SEARCHING GENES AGAINST PUBMED DATABASE...')
            for i in tqdm(range(0, len(dfli))):
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

                ncbi = bs4.BeautifulSoup(res.text,"html.parser")

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

                #sometimes if there are no articles, the API returns a count in the thousands. This filters those out.
                if int(once) > 5000:
                    once = 0
                dfli_count[i] = once

            #Change all str to numbers in the count list
            for j in range(0, len(dfli_count)):
                dfli_count[j] = int(dfli_count[j])


            #append the list of counts column to the end of the original gene list column
            df.insert(len(df.columns), 'Pubmed Count', dfli_count, True)
            df.drop(df.columns[[0]], axis=1, inplace=True)
            df.sort_values(by="Pubmed Count", inplace=True, ascending=False)
            df.reset_index(inplace=True)
            df.drop(df.columns[0], axis=1, inplace=True)
            #output to an excel doc named gene-list-output + the keyword searched
            searchterm = keyword.replace(' ', '-')
            searchterm = keyword.strip()
            output = df.to_excel( outname + '.xlsx')

            return output
        
    else:
        dfli = df.values.tolist()
        #preallocate a list of counts with zeros

        dfli_count = [0] * len(dfli)

        #loop through gene list
        print('SEARCHING GENES AGAINST PUBMED DATABASE...')
        for i in tqdm(range(0, len(dfli))):
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
            ncbi = bs4.BeautifulSoup(res.text,"html.parser")
        #find all the HTML elements with <count></count> and take the first one (this is the one with the actual number of papers)
            count = ncbi.select('Count')
            if count == []:
                dfli_count[i] = 0
                continue
            once = str(count[0])
        #remove the <count></count> tags
            once = once.replace('<count>', '')
            once = once.replace('</count>', '')
        # sometimes if there are no articles, the API returns a count in the thousands. This filters those out.
            if int(once) > 5000:
                once = 0
        #add to the respective place in list of counts    
            dfli_count[i] = once

        #Change all str to numbers in the count list
        for j in range(0, len(dfli_count)):
            dfli_count[j] = int(dfli_count[j])


        #append the list of counts column to the end of the original gene list column
        df.insert(len(df.columns), 'Pubmed Count', dfli_count, True)
        df.sort_values(by="Pubmed Count", inplace=True, ascending=False)
        df.reset_index(inplace=True)
        df.drop(df.columns[0], axis=1, inplace=True)
        #output to an excel doc named gene-list-output + the keyword searched
        searchterm = keyword.replace(' ', '-')
        searchterm = keyword.strip()
        output = df.to_excel(outname + '.xlsx')
        
        return output



def main():
    my_parser = argparse.ArgumentParser(description="This script converts all ensemble IDs into conventional"
                                                    "gene names, then takes a gene list plus a keyword of your choice"
                                                    "and searches each gene plus your keyword into pubmed and "
                                                    "returns the number of articles published for that particular search term")
    my_parser.add_argument('Path', metavar='path', type=str, help='the path to the gene list excel file')
    my_parser.add_argument('Search Term', metavar='search-term', type=str,
                           help='the search term to search pubmed alongside each gene in the list')
    my_parser.add_argument('Output File Name', metavar='output-file-name', type=str, help='the name of the excel file you would like to output')
    args = my_parser.parse_args()
    input_path = args.Path
    ens = pd.read_excel(sys.argv[1])
    if ens.iloc[0, 0].startswith('ENS'):
        ensembl_id_finder(sys.argv[1])
        pubmedscrape('ensembl-id-temp.xlsx', str(sys.argv[2]), str(sys.argv[3]))
        os.remove('ensembl-id-temp.xlsx')
    else:
        pubmedscrape(str(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3]))
main()

