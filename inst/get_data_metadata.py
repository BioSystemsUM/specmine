# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 11:07:39 2017

@author: sara
"""

#####IMPORT METABOLIGHTS FILES AND METADATA#####

#List of studies in metabolights:
def metabolights_studies_list():
  from isatools.net import mtbls
  metabolights_studies_list=mtbls.get_mtbls_list()
  return(metabolights_studies_list)

#GET METADATA INFO INTO A FILE:
def get_metabolights_metadata(studyID, study_directory, investigation_file):
    from isatools import isatab
    from os.path import join
    
    with open(investigation_file, encoding='utf-8') as f:
        ISA_object = isatab.load(f)
    
    met_vars=get_factor_names_m(study_directory) #Metadata classes

    #Get full metadata values:
    metadata_full=get_factors_summary_m(ISA_object)
    s=","
    header="sample,"+s.join(met_vars)
    f=open(join(study_directory,'metadata.csv'),'w')
    f.write(header)
    f.write('\n')
    for metadata_sample in metadata_full:
        line="\""+metadata_sample["sample"]+"\""
        for met_var in met_vars:
            line=line+","+"\""+metadata_sample[met_var]+"\""
        f.write(line)
        f.write('\n')
    f.close()
    return (join(study_directory,'metadata.csv'))


#LOAD DATA FILES:
def load_metabolights_dataFiles(studyID, directory):
    from os.path import join
    from os import rename, listdir
    import ftplib
    
    #Get data files names
    data_files_str=slice_data_files(directory)
    data_files=list()
    for file in data_files_str:
        for f in file['data_files']:
            data_files.append(f)
            
    #Download data files
    ftp = ftplib.FTP('ftp.ebi.ac.uk')
    response = ftp.login()
    if '230' in response:  # 230 means Login successful
        ftp.cwd('pub')
        ftp.cwd('databases')
        ftp.cwd('metabolights')
        ftp.cwd('studies')
        ftp.cwd('public')
        ftp.cwd(studyID)
        for file_name in data_files:
            print(join(directory, file_name))
            with open(join(directory, file_name), 'wb') as out_file:
                ftp.retrbinary('RETR ' + file_name, out_file.write)
        ftp.close()
    
    #Write file that says what files belong to what samples:
    f=open(join(directory, 'samples_files.csv'),'w')
    for file_name_sample in data_files_str:
        file_name=file_name_sample["data_files"]
        file_sample=file_name_sample["sample"]
        line=file_sample
        for n in file_name:
            line=line+","+n
        f.write(line)
        f.write("\n")
    f.close()
    return directory


#LOAD COMPLETE DATA :
def load_data_metadata_metabolights(studyID, directory):
    from isatools.net import mtbls
    from os.path import join, exists
    from os import makedirs
    
    #Create directory to save files
    directory=join(directory, studyID)
    if not exists(directory):
        makedirs(directory)
    
    #Get information files on the study
    investigation_file=join(directory, "i_Investigation.txt")
    mtbls.get(studyID, target_dir=directory)
    
    #Download data files:
    print(load_metabolights_dataFiles(studyID, directory))
    
    #Create metadata file:
    print(get_metabolights_metadata(studyID, directory, investigation_file))
    
    return()

            
