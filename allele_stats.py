# !/usr/bin/env python
# encoding: utf-8

"""
Goal: To calculate various allele count and frequency statistics described in 
Racimo et al. 2017 (Mol. Biol. Evol.). These stats include U20, 
U50, and Q95, which are all used to detect candidate sites of adaptive 
introgression.
________________________________________________________________________________
Input: VCF file (with the expected nine intro columns from 'CHROM' through 
'FORMAT'), population file, genomic windows file
________________________________________________________________________________
Setup:
Use grep -n "#CHROM" filename.vcf to identify the line number of the VCF file 
header line. Subtract one and use that value for the --skipRows flag.

Create bed file with desired genomic windows (use -s flag if you want sliding 
windows):
*create genome file with structure <chromName><TAB><chromSize>. I show how to
create this file from the reference genome fasta index file (fasta.fai).

awk -v OFS='\t' {'print $1,$2'} genome.fasta.fai > genome_file.txt
awk -v OFS='\t' {'print $1,$2'} xipho_ref.fasta.fai > genome_file.txt

From the fasta index file for the reference genome, 
*create windows bed file
*input windows to program
bedtools makewindows -g genome_file.txt -w 2000000 -s 1900000 > windows.bed
bedtools makewindows -g genome_file.txt -w 2000000 > windows_test.bed

________________________________________________________________________________
Output: A text file with VCF data for "qualifying sites" (those sites that meet 
the conditions necessary to calculate U20, U50, and Q95 stats) and new columns 
with allele count and frequency data for the target population at each site.
________________________________________________________________________________
Note: This script filters out sites with >= 50% missingness PER POPULATION to 
help ensure data quality.
________________________________________________________________________________
Formally, in Racimo et al. 2017:
***U STATS***
"U_A_B_C(w,x,y): Number of sites in which any allele is at a frequency lower 
than w in panel A, higher than x in panel B and equal to y in panel C." In this 
script, we use w=1% and y=100%.
***Q95 STAT***
"Q95_A_B_C(w,y): 95% quantile of the distribution of derived allele frequencies 
in panel B, for sites where the derived allele is at a frequency lower 
than w in panel A and equal to y in panel C." In this script, we use w=1% and 
y=100%.
________________________________________________________________________________
System requirements:
Currently, this scripts loads the whole VCF file into memory, so you need a lot 
of memory for large VCF files. I ran this program with a 25 GB VCF file on the 
Louisiana State University cluster in about 20 minutes using ~155 GB of memory. 
________________________________________________________________________________
Command to run on LSU cluster:

#!/bin/bash
#PBS -A hpc_genomes3
#PBS -l nodes=1:ppn=28
#PBS -l walltime=72:00:00
#PBS -q bigmem
#PBS -N allele_stats

cd $PBS_O_WORKDIR

python allele_stats.py --vcfFile test.vcf --skipRows 81 \
--windowFile windows.bed --popKey popKey.txt --popA belem \
--popB xingu --popC tapajos

________________________________________________________________________________
Copyright 2024 Andre E. Moncrieff. All rights reserved.

"""

import argparse
#from ast import AsyncFunctionDef
#from curses.ascii import alt
#from logging import setLogRecordFactory
#from turtle import up, update
import pandas
from collections import Counter
import numpy

pandas.options.mode.chained_assignment = None  # default='warn'


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcfFile", required=True,
                        help="Enter the file name (+ .vcf)",
                        type=str)
    parser.add_argument("--skipRows", required=True,
                        help="Enter number of rows to skip before VCF header",
                        type=int)
    parser.add_argument("--windowFile", required=True,
                        help="Enter the file name (+.bed)",
                        type=str)
    parser.add_argument("--popKey", required=True,
                        help="Enter the population file name (+ .txt)",
                        type=str)
    parser.add_argument("--popA", required=True,
                        help="Enter the population A basename",
                        type=str)
    parser.add_argument("--popB", required=True,
                        help="Enter the population B basename",
                        type=str)
    parser.add_argument("--popC", required=True,
                        help="Enter the population C basename",
                        type=str) 
    args = parser.parse_args()
    return args


def read_in_vcf(txt_file_dataframe, rows_to_skip):
    raw_dataframe = pandas.read_csv(txt_file_dataframe,
                                    sep='\t',
                                    encoding = "ISO-8859-1",
                                    dtype=str, skiprows=rows_to_skip)
    updated_dataframe = raw_dataframe.rename(columns={'#CHROM': 'CHROM'})
    return updated_dataframe


def read_in_popKey(txt_file_dataframe):
    raw_dataframe = pandas.read_csv(txt_file_dataframe,
                                    sep='\t',
                                    encoding = "ISO-8859-1",
                                    dtype=str)
    return raw_dataframe

def read_in_windows(windowFile, column_names):
    raw_dataframe = pandas.read_csv(windowFile,
                                    sep='\t',
                                    encoding = "ISO-8859-1",
                                    dtype=str,
                                    names=column_names)
    return raw_dataframe



def pop_dict_maker(pop_df):
    keys = pop_df['Individual'].tolist()
    values = pop_df['Population'].tolist()    
    dictionary = dict(zip(keys, values))
    return dictionary


def rows_to_list(vcf_df):
    list_of_rows = vcf_df.values.tolist()
    return list_of_rows


def trim_lists(list_of_rows):
    #takes out the first 9 columns of VCF: '#CHROM' through 'FORMAT'
    #The point is to get a dataframe of the individuals only
    trimmed_list_of_rows = []
    for item in list_of_rows:
        trimmed_row= item[9:]
        trimmed_list_of_rows.append(trimmed_row)
    return trimmed_list_of_rows


def genotype_list(list_of_rows):
    #get the genotypes only for each individual
    genotypes_by_row = []
    for row in list_of_rows:
        rowlist = []
        for item in row:
            genotype = item.split(':')[0]
            rowlist.append(genotype)
        genotypes_by_row.append(rowlist)
    return genotypes_by_row


def find_pop_slices(sorted_vcf_df3, popA_string, popB_string, popC_string):
    #print(sorted_vcf_df3)
    untrimmed_list_of_basenames = []
    for header in sorted_vcf_df3.columns:
        pop_basename = header.split('_')[0]
        untrimmed_list_of_basenames.append(pop_basename)
    trimmed_list_of_basenames = untrimmed_list_of_basenames[9:]
    popA_first_index = trimmed_list_of_basenames.index(popA_string)
    popB_first_index = trimmed_list_of_basenames.index(popB_string)
    popC_first_index = trimmed_list_of_basenames.index(popC_string)
    popA_count = trimmed_list_of_basenames.count(popA_string)
    popB_count = trimmed_list_of_basenames.count(popB_string)
    popC_count = trimmed_list_of_basenames.count(popC_string)
    popA_second_index = popA_first_index + popA_count
    popB_second_index = popB_first_index + popB_count
    popC_second_index = popC_first_index + popC_count
    return popA_first_index, popA_second_index, \
           popB_first_index, popB_second_index, \
           popC_first_index, popC_second_index


def polarize_popA_popC(list_of_genotypes, slices):   
    popA_first_index = slices[0]
    popA_second_index = slices[1]
    popB_first_index = slices[2]
    popB_second_index = slices[3] 
    popC_first_index = slices[4]
    popC_second_index = slices[5]
    row_indicator = []
    for item in list_of_genotypes:
        #each item is a list representing an entire row of genotypes from all 
        #individuals in the dataset
        popA_genotypes = item[popA_first_index:popA_second_index]
        popB_genotypes = item[popB_first_index:popB_second_index]
        popC_genotypes = item[popC_first_index:popC_second_index]
        number_of_popA_alleles = len(popA_genotypes)*2
        number_of_popB_alleles = len(popB_genotypes)*2
        number_of_popC_alleles = len(popC_genotypes)*2
        popA_set = set(popA_genotypes)
        popC_set = set(popC_genotypes)
        homosetref = set(['0/0', '0|0', "./."])
        homosetalt = set(['1/1', '1|1', "./."])
        popA_indiv_alleles = []
        popB_indiv_alleles = []
        popC_indiv_alleles = []        
        for unit in popA_genotypes: #unit = genotype
            popA_allele1 = unit.split('/')[0]
            popA_allele2 = unit.split('/')[1]
            popA_indiv_alleles.append(popA_allele1)
            popA_indiv_alleles.append(popA_allele2)
        for block in popB_genotypes: #block = genotype
            popB_allele1 = block.split('/')[0]
            popB_allele2 = block.split('/')[1]
            popB_indiv_alleles.append(popB_allele1)
            popB_indiv_alleles.append(popB_allele2)
        for element in popC_genotypes: #element = genotype
            popC_allele1 = element.split('/')[0]
            popC_allele2 = element.split('/')[1]
            popC_indiv_alleles.append(popC_allele1)
            popC_indiv_alleles.append(popC_allele2)
        if popA_set.issubset(homosetref) == True and \
            popC_set.issubset(homosetalt) == True:
            popA_indiv_allele_counts = (Counter(popA_indiv_alleles))
            popB_indiv_allele_counts = (Counter(popB_indiv_alleles))
            popC_indiv_allele_counts = (Counter(popC_indiv_alleles))
            popA_missing_allele_count = popA_indiv_allele_counts['.']
            popB_missing_allele_count = popB_indiv_allele_counts['.']
            popC_missing_allele_count = popC_indiv_allele_counts['.']
            if popA_missing_allele_count < (number_of_popA_alleles/2) and \
            popB_missing_allele_count < (number_of_popB_alleles/2) and \
            popC_missing_allele_count < (number_of_popC_alleles/2): 
                # requirement of less than 50% missingness per population
                row_indicator.append("popA_popC_polarized_site")
            else:
                row_indicator.append("no")
        else: 
            row_indicator.append("no")
    return row_indicator


def add_column_pandas(vcf_df, indicator_row):
    vcf_df['Indicator'] = indicator_row
    return vcf_df


def filter_dataframe(new_vcf_df):
    filtered_df = new_vcf_df[new_vcf_df['Indicator'].isin(
                                                  ["popA_popC_polarized_site"])]
    return filtered_df


def delete_column(filtered_df):
    del filtered_df['Indicator']
    return filtered_df


def popB_allele_stats(final_list_of_genotypes, slices):
    popB_alt_allele_counts = []
    popB_missing_allele_counts = []
    popB_alt_allele_freq = []
    popB_alt_allele_adj_freq = []
    popB_first_index = slices[2]
    popB_second_index = slices[3]
    for row in final_list_of_genotypes:
        popB_genotypes = row[popB_first_index:popB_second_index]
        popB_indiv_alleles = []
        for genotype in popB_genotypes:
            allele1 = genotype.split('/')[0]
            allele2 = genotype.split('/')[1]
            popB_indiv_alleles.append(allele1)
            popB_indiv_alleles.append(allele2)
        total_allele_count_popB = len(popB_indiv_alleles)
        counts = (Counter(popB_indiv_alleles))
        alternate_allele_count = counts['1']
        missing_allele_count = counts['.']
        alternate_allele_freq = round(alternate_allele_count/
                                      total_allele_count_popB, 3)
        alternate_allele_adj_freq = round(alternate_allele_count/
                            (total_allele_count_popB-missing_allele_count), 3)
        popB_alt_allele_counts.append(alternate_allele_count)
        popB_missing_allele_counts.append(missing_allele_count)
        popB_alt_allele_freq.append(alternate_allele_freq)
        popB_alt_allele_adj_freq.append(alternate_allele_adj_freq)
    return popB_alt_allele_counts, popB_missing_allele_counts, \
           popB_alt_allele_freq, popB_alt_allele_adj_freq


def create_full_window_lists(final_vcf_dataframe, window_df):
    list_of_window_rows = window_df.values.tolist()
    
    #print(list_of_window_rows)
    updated_list_window_rows = []
    CHROM = final_vcf_dataframe['CHROM'].tolist()  
    POS = final_vcf_dataframe['POS'].tolist()  
    allele_freq = final_vcf_dataframe['popB_alt_allele_freq'].tolist()
    #allele_adj_freq = final_vcf_dataframe['popB_alt_allele_adj_freq'].tolist()
    updated_list_window_rows = []
    for row in list_of_window_rows:
        row_list = []
        row_list.extend(row)
        for i in range(len(CHROM)):
            if CHROM[i] == row[0] and int(POS[i]) in range(int(row[1]),int(row[2])):
                row_list.append(str(allele_freq[i]))
            else:
                pass
        updated_list_window_rows.append(row_list)
    return updated_list_window_rows


def stats_from_window_lists(updated_window_rows):
    trimmed_window_rows = []
    for row in updated_window_rows:
        trimmed_row = row[3:]
        trimmed_window_rows.append(trimmed_row)
    rows_with_final_stats = []
    for row in trimmed_window_rows:
        row_values = []
        windows_with_no_informative_alleles = ['0','0','0', '0']
        if len(row) == 0:
            row_values.extend(windows_with_no_informative_alleles)
        elif len(row) > 0:
            float_row = [float(i) for i in row]
            sites_over_20p_freq = [i for i in float_row if i>0.2]
            U20 = len(sites_over_20p_freq)
            sites_over_50p_freq = [i for i in float_row if i>0.5]
            U50 = len(sites_over_50p_freq)
            np_array = numpy.array(float_row)
            Q95 = round(numpy.percentile(np_array, 95), 3)
            informative_sites = len(row)
            row_values.append(U20)
            row_values.append(U50)
            row_values.append(Q95)
            row_values.append(informative_sites)
        rows_with_final_stats.append(row_values)
    return rows_with_final_stats
            
        
def assemble_final_list_of_rows(updated_window_rows, rows_with_final_stats):
    final_list_of_rows = []
    trimmed_window_rows = []
    for row in updated_window_rows:
        trimmed_row = row[:3]
        trimmed_window_rows.append(trimmed_row)
    for i in range(len(updated_window_rows)):
        chrom_and_window_list = trimmed_window_rows[i]
        final_stats_list = rows_with_final_stats[i]
        complete_row = chrom_and_window_list + final_stats_list
        final_list_of_rows.append(complete_row)
    return final_list_of_rows


def list_to_df(final_list_of_rows):
    colnames = ['chromosome', 'start', 'end', 'U20', 'U50', 'Q95', 
                'informative_sites']
    dataframe = pandas.DataFrame(final_list_of_rows, columns = colnames)
    return dataframe


def main():
    #create args object
    args = parser()    
    rows_to_skip = args.skipRows
    print("\n" + "RUNNING: reading in VCF file")
    vcf_df = read_in_vcf(args.vcfFile, rows_to_skip)
    pop_df = read_in_popKey(args.popKey)
    column_names = ['chrom', 'start', 'end']
    window_df = read_in_windows(args.windowFile, column_names) 
    popA_string = args.popA
    popB_string = args.popB
    popC_string = args.popC
    pop_dict = pop_dict_maker(pop_df) # make a dictionary of old (keys) and new 
    # (values) population names based on the popKey file
    renamed_vcf_df2 = vcf_df.rename(columns=pop_dict) # rename dataframe headers
    sorted_vcf_df3 = renamed_vcf_df2.reindex(sorted(renamed_vcf_df2.columns), 
                                             axis=1) # sort the columns in 
                                                     # dataframe by population
    
    list_of_rows = rows_to_list(sorted_vcf_df3) # convert each dataframe row 
    # into a list
    trimmed_list_of_rows = trim_lists(list_of_rows) # remove first nine elements 
    # from each row list (corresponding to the CHROM through FORMAT elements of 
    # the VCF file)
    list_of_genotypes = genotype_list(trimmed_list_of_rows) # create a list of 
    # all genotypes for each site
    slices = find_pop_slices(sorted_vcf_df3, popA_string, popB_string, 
                             popC_string)
    print("\n" + "RUNNING: identifying sites with fixed reference alleles in")
    print("         population A and fixed alternate alleles in population C")
    row_indicator_list = polarize_popA_popC(list_of_genotypes, slices)
    indicator_column_vcf_df4 = add_column_pandas(sorted_vcf_df3, 
                                                 row_indicator_list)
    filtered_vcf_df5 = filter_dataframe(indicator_column_vcf_df4)
    cleaned_vcf_df6 = delete_column(filtered_vcf_df5)
    filtered_list_of_rows = rows_to_list(cleaned_vcf_df6)
    filt_trimmed_list_of_rows = trim_lists(filtered_list_of_rows)
    final_list_of_genotypes = genotype_list(filt_trimmed_list_of_rows)
    print("\n" + "RUNNING: calculating allele statistics for population B ")
    freq_and_counts_quad = popB_allele_stats(final_list_of_genotypes, slices)
    popB_alt_allele_count = freq_and_counts_quad[0]
    popB_missing_allele_count = freq_and_counts_quad[1]
    popB_alt_allele_freq = freq_and_counts_quad[2]
    popB_alt_allele_adj_freq = freq_and_counts_quad[3]
    cleaned_vcf_df6['popB_alt_allele_count'] = popB_alt_allele_count
    cleaned_vcf_df6['popB_missing_allele_count'] = popB_missing_allele_count
    cleaned_vcf_df6['popB_alt_allele_freq'] = popB_alt_allele_freq
    cleaned_vcf_df6['popB_alt_allele_adj_freq'] = popB_alt_allele_adj_freq
    print("\n"+"RUNNING: saving 1st output file -> allele_stats_by_site.csv")
    cleaned_vcf_df6.to_csv('allele_stats_by_site.csv', sep=',', index=False)
    print("Preview of allele_stats_by_site.csv:")
    print(cleaned_vcf_df6)
    updated_window_rows = create_full_window_lists(cleaned_vcf_df6, window_df)
    rows_with_final_stats = stats_from_window_lists(updated_window_rows)
    final_list_of_rows = assemble_final_list_of_rows(updated_window_rows, 
                                                     rows_with_final_stats)
    allele_stats_by_window_df = list_to_df(final_list_of_rows)
    print("\n"+"RUNNING: saving 2nd output file -> allele_stats_by_window.csv")
    
    allele_stats_by_window_df.to_csv('allele_stats_by_window.csv', sep=',', 
                                     index=False)
    print("Preview of allele_stats_by_window.csv:")
    print(allele_stats_by_window_df)
    print("\n"+"DONE"+"\n")

if __name__ == '__main__':
    main()