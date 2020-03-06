#!/usr/bin/env python

# Gene rate Submission.sh within working dir
__author__ = "Niu Du"
__email__ = "ndu@lji.org"

import argparse,csv
import json
import pandas as pd
from RNA_SEQ_Func import gen_Submission

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--conf', type=str)
    parser.add_argument('-i', '--fastq_table',type=str)
    parser.add_argument('-n', '--threads', type=int, default = 4)
    opts = parser.parse_args()

    conf_file = opts.conf
    fastqfile = opts.fastq_table
    n_cpu = opts.threads

    dict_conf = read_json(conf_file)
    fastq_table = pd.read_csv(fastqfile)
    fastq_table = fastq_table.set_index(fastq_table.columns[0])

    gen_Submission(fastq_table,dict_conf,thread = n_cpu)   
    
def read_json(file = 'conf_RNA_Seq.json'):
    with open(file) as json_file:
        conf = json.load(json_file)
    return conf


if __name__ == "__main__":
    main()