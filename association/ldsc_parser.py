import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(prog="LDSC-parser", description=("Parses "
    "ldsc.py --rg analysis.log files for Lambda GC, Intercept, Ratio and "
    "summary table."))

parser.add_argument(
        '-f', '--file_log',
        action="store",
        dest='f',
        required=True,
        type=str,
        help=('Path to LDSC .log file. Default: %(default)s'))

parser.add_argument(
        '-o', '--outdir',
        action="store",
        dest='o',
        required=True,
        type=str,
        help=('Directory where results will be saved. '
            'Default: %(default)s'))

parser.add_argument(
        '-n', '--name',
        action="store",
        dest='name',
        required=True,
        type=str,
        help=('Analysis name. Default: %(default)s'))

options = parser.parse_args()
fname = options.f
directory = options.o
name = options.name

res = []
summary = []
summary_indicator = False
header = False


#fname = '/homes/hannah/data/ukbb/ukb-hrt/gwas/ldsc_slices_fd.log'
with open(fname) as f:
    content = f.readlines()
    content = [x.strip() for x in content]
    for line in content:
        if line.startswith('Heritability of phenotype'):
            phenotype = line.split(" ")[-1].split("/")[0]
        if line.startswith('Total'):
            h2 = line.split(" ")[-2]
        if line.startswith('Lambda GC:'):
            lambdaGC = line.split(" ")[-1]
        if line.startswith('Mean Chi'):
            meanChi = line.split(" ")[-1]
        if line.startswith('Intercept:'):
            intercept = line.split(" ")[-2]
        if line.startswith('Ratio'):
            ratio = line.split(" ")[-2]
            res.append((phenotype, h2, lambdaGC, meanChi, intercept, ratio))
        if line.startswith('Analysis'):
            summary.pop()
            summary_indicator = False
        if summary_indicator:
            if header:
                summary_columns = line.split()
                header = False
            else:
                summary.append(line.split())
        if line.startswith('Summary of Genetic Correlation Results'):
            summary_indicator = True
            header = True
    ldsc = pd.DataFrame.from_records(data=res,
            columns=['Phenotype', 'observed h2', 'LambdaGC', 'MeanChi2',
                'Intercept', 'Ratio'])
    summary = pd.DataFrame(data=summary, columns=summary_columns)

summary.to_csv('{}/LDSC_genetic_correlationP1vsAll_{}.csv'.format(directory,
    name), header=True, sep=",", index=False)

ldsc.to_csv('{}/LDSC_GC_{}.csv'.format(directory, name), header=True, sep=",",
        index=False)
