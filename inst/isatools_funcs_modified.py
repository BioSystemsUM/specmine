# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:44:18 2017

@author: sara
"""

#pip install isatools

from __future__ import absolute_import
import ftplib
import glob
import logging
import os
import pandas as pd
import tempfile
import shutil
import re

from isatools import config
from isatools import isatab
from isatools.convert import isatab2json
from isatools.model import OntologyAnnotation
from isatools.net import mtbls

EBI_FTP_SERVER = 'ftp.ebi.ac.uk'
MTBLS_BASE_DIR = '/pub/databases/metabolights/studies/public'
INVESTIGATION_FILENAME = 'i_Investigation.txt'

logging.basicConfig(level=config.log_level)
log = logging.getLogger(__name__)

# REGEXES
_RX_FACTOR_VALUE = re.compile('Factor Value\[(.*?)\]')

def get_factor_names_m(directory):
    #Get the metadata of the study:
    factors = set()
    for table_file in glob.iglob(os.path.join(directory, '[a|s]_*')):
        with open(os.path.join(directory, table_file), encoding='utf-8') as fp:
            df = isatab.load_table(fp)
            factors_headers = [header for header in list(df.columns.values) if _RX_FACTOR_VALUE.match(header)]
            for header in factors_headers:
                factors.add(header[13:-1])
    return factors

def get_factors_summary_m(ISA):
    #Get full metadata in a list of dictionaries: each dictionary is the 
    #metadata information for a sample
    all_samples = []
    for study in ISA.studies:
        all_samples.extend(study.samples)
    samples_and_fvs = []
    for sample in all_samples:
        sample_and_fvs = {
                "sources": ';'.join([x.name for x in sample.derives_from]),
                "sample": sample.name,
            }
        for fv in sample.factor_values:
            if isinstance(fv.value, (str, int, float)):
                fv_value = fv.value
            elif isinstance(fv.value, OntologyAnnotation):
                fv_value = fv.value.term
            sample_and_fvs[fv.factor_name.name] = fv_value
        samples_and_fvs.append(sample_and_fvs)
    df = pd.DataFrame(samples_and_fvs)
    nunique = df.apply(pd.Series.nunique)
    cols_to_drop = nunique[nunique == 1].index
    df = df.drop(cols_to_drop, axis=1)
    return df.to_dict(orient='records')

def slice_data_files(dir, factor_selection=None):
    #Gets list of dictionaries with each one being the data file string(s) for
    #the sample
    results = list()
    # first collect matching samples
    for table_file in glob.iglob(os.path.join(dir, '[a|s]_*')):
        log.info("Loading {}".format(table_file))
        with open(table_file, encoding='utf-8') as fp:
            df = isatab.load_table(fp)
            if factor_selection is None:
                matches = df['Sample Name'].items()
                for indx, match in matches:
                    sample_name = match
                    if len([r for r in results if r['sample'] == sample_name]) == 1:
                        continue
                    else:
                        results.append(
                            {
                                "sample": sample_name,
                                "data_files": []
                            }
                        )
            else:
                for factor_name, factor_value in factor_selection.items():
                    if 'Factor Value[{}]'.format(factor_name) in list(df.columns.values):
                        matches = df.loc[df['Factor Value[{}]'.format(factor_name)] == factor_value]['Sample Name'].items()
                        for indx, match in matches:
                            sample_name = match
                            if len([r for r in results if r['sample'] == sample_name]) == 1:
                                continue
                            else:
                                results.append(
                                    {
                                        "sample": sample_name,
                                        "data_files": [],
                                        "query_used": factor_selection
                                    }
                                )
    # now collect the data files relating to the samples
    for result in results:
        sample_name = result['sample']
        for table_file in glob.iglob(os.path.join(dir, 'a_*')):
            with open(table_file, encoding='utf-8') as fp:
                df = isatab.load_table(fp)
                data_files = list()
                table_headers = list(df.columns.values)
                sample_rows = df.loc[df['Sample Name'] == sample_name]
                if 'Raw Spectral Data File' in table_headers:
                    data_files = sample_rows['Raw Spectral Data File']
                elif 'Free Induction Decay Data File' in table_headers:
                    data_files = sample_rows['Free Induction Decay Data File']
                result['data_files'] = [i for i in list(data_files) if str(i) != 'nan']
    return results