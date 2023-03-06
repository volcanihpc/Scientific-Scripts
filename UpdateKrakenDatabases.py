#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
## Kraken2 genome and protein download script using Multiprocessing.
This script takes an optional command-line argument for the number
of CPUs for multiprocessing and for building the database.
By default, 8 CPUs will be allocated.

Standard - archaea, bacteria, viral, human, UniVec_Core
"""
#import modules
import os
import re
import sys
import subprocess
import multiprocessing
import pandas as pd
from Bio import SeqIO

def download_file(url):
    ## download ~31654 files in parallel by the amount of allocated cpus,and in silent mode
    subprocess.call(['wget', '-N', url])

def main(cpus=None):
    try:
        libraries = ["archaea" , "bacteria" , "viral" , "human" , "UniVec_Core"]
        for l in libraries:
            os.makedirs(f'library/{l}')
    except FileExistsError:
        pass
       
    ## get the current working directory
    if cpus == None:
        cpus = 8

    cwd = os.getcwd()

    ## function to process ftp url file that is created from assembly files
    def process_url_file(inputurlfile, file_suffix):
        url_file=open(inputurlfile,'r')
        pool = multiprocessing.Pool(processes=cpus)
        url_list = []
        for line in url_file:
            url=line.rstrip('\n').split(',')
            ftp_url= url[0]+'/'+url[1]+'_'+url[2]+'_'+file_suffix
            url_list.append(ftp_url)
        ## Download the files in the fna format
        print(len(url_list))
        pool.map(download_file, url_list)
        ## unzip the files
        subprocess.call("gunzip *.gz",shell=True)
        pool.close()
        pool.join()
        return
    
    ## Download Assembly File and return a Pandas DataFrame:
    def download_assembly_file(url, filename):
        ## Download the file using wget system call
        subprocess.call("wget "+url, shell=True)
        ## Reformat the file to pandas-friendly format
        subprocess.call(f"sed -i '1d' {filename}",shell=True)
        subprocess.call(f"sed -i 's/^# //' {filename}", shell=True)
        ## Read the file as a dataframe - using read_table
        ## Use read_table if the column separator is tab
        assembly_sum = pd.read_csv(filename, sep='\t', dtype='unicode')

        return assembly_sum


    ## function to download bacterial sequences
    def download_bacterial_genomes(outfile='outfile.txt'):
        assembly_file = 'assembly_summary.txt'
        ncbi_url = f'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/{assembly_file}'
        if os.path.exists(assembly_file):
           os.remove(assembly_file)
        ## Download the file using wget sysyem call
        assembly_df = download_assembly_file(ncbi_url, assembly_file)

        ## filter the dataframe and save the URLs of the complete genomes in a new file
        my_df=assembly_df[(assembly_df['version_status'] == 'latest') &
                    (assembly_df['assembly_level']=='Complete Genome') 
                    ]
        my_df=my_df[['ftp_path','assembly_accession','asm_name']]
        ## output_file.write
        my_df.to_csv(outfile,mode='w',index=False,header=None)
        print(my_df)
        process_url_file(outfile, 'genomic.fna.gz')
        process_url_file(outfile, 'protein.faa.gz')
        return
        
    ## function to download reference genomes 
    ## this function downloads latest version human reference genome by default 
    def download_refseq_genome(path, taxid=None,outfile='refseq_genome.txt'):
        os.chdir(path)
        assembly_file='assembly_summary.txt'
        ncbi_url=f"https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/{assembly_file}"

        if os.path.exists(assembly_file):
            os.remove(assembly_file)
        assembly_df = download_assembly_file(ncbi_url, assembly_file)

        my_df=assembly_df[(assembly_df['taxid'] == taxid) &
                        ((assembly_df['assembly_level'] == 'Complete Genome') 
                        )]
        my_df=my_df[['ftp_path','assembly_accession','asm_name']]
        ## Process the newly created file and download genomes from NCBI website
        my_df.to_csv(outfile,mode='w',index=False,header=None)
        print(my_df)
        process_url_file(outfile, 'genomic.fna.gz')
        process_url_file(outfile, 'protein.faa.gz')
        return

    ## Human
    print('Downloading human genome'+'\n')
    ## change argument in the following function if you want to download other reference genomes
    ## taxonomy ID 9606 (human) should be replaced with taxonomy ID of genome of interest
    download_refseq_genome('library/human', '9606','human_genome_url.txt')
    
    ## viral
    os.chdir('../../')
    os.chdir('library/viral')
    print('Downloading viral genomes and protein'+'\n')
    subprocess.call('wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz',shell=True)
    subprocess.call('wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz',shell=True)
    subprocess.call('gunzip *.gz',shell=True)

    ## UniVec_Core
    os.chdir('../../')
    os.chdir('library/UniVec_Core')
    print('Downloading UniVec_Core'+'\n')
    subprocess.call('wget https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core', shell=True)

    ## archaea
    os.chdir('../../')
    os.chdir('library/archaea')
    print('Downloading archaeal genomes'+'\n')
    for r in range(1, 4):
        for k in range(1, 3):
            subprocess.call(f'wget https://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/archaea.{r}.{k}.genomic.fna.gz',shell=True)
    for r in range(1, 6):
        subprocess.call(f"wget https://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/archaea.nonredundant_protein.{r}.protein.faa.gz", shell=True)
    subprocess.call('gunzip *.gz',shell=True)

    ## bacteria
    os.chdir('../../')
    os.chdir('library/bacteria')
    print('Downloading bacterial genomes'+'\n')
    download_bacterial_genomes('bacterial_complete_genome_url.txt')

    os.chdir('../../')

if __name__ == "__main__":
    if len(sys.argv) > 1:
        try:
            main(int(sys.argv[1]))
        except ValueError as e:
            print("Enter the number of cpus to use. default is 8.")
            exit()
    else:
        main()
    
