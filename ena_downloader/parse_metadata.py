import xml.etree.ElementTree as xml
import pandas as pd
import urllib

def get_study_metadata_request(study_summary_file):

    study_ids = pd.read_table(study_summary_file, header = 0, index_col = None, sep = '\t')
    ERSs = list(study_ids['secondary_sample_accession'])
    ERSs = [int(x[3:]) for x in ERSs]
    ERSs.sort()

    batch = ''
    i = 0
    while i < len(ERSs):
        left = ERSs[i]
        right = ERSs[i]
        if ERSs[i] != ERSs[-1]:
            while ERSs[i+1] == ERSs[i]+1:
                right = ERSs[i + 1]
                i+=1
                if ERSs[i] == ERSs[-1]:
                    break
        if left == right:
            batch += 'ERS{},'.format(left)
        else:
            batch += 'ERS{}-ERS{},'.format(left,right)
        i += 1

    batch = batch.strip(',')
    request = 'http://www.ebi.ac.uk/ena/data/view/{}&display=xml'.format(batch)
    metadata = urllib.request.urlopen(request).read()
    file = open('.'.join(study_summary_file.split(sep='.')[:-1]) + '_metadata.txt', 'wb')
    file.write(metadata)
    file.close()
    return(request)

def determine_runs_sample_parameter(metadata_table, variable_name, ID = {}):
    root = xml.parse(metadata_table).getroot()
    samples = root.findall("./SAMPLE")
    vars = [variable_name]

    for i in samples:
        links = i.findall('SAMPLE_LINKS/')
        err = [x[0][1].text for x in links][2]
        ID[err] = {variable_name: None}
        if err in ID.keys():
            meta = i.findall('SAMPLE_ATTRIBUTES/')
            for j in meta:
                var = j[0].text
                if var in vars:
                    ID[err][var] = j[1].text
    return(ID)

def show_values_of_variables(IDs, variable_name):
    values = [IDs[k][variable_name] for k in IDs.keys()]
    values = {x:values.count(x) for x in set(values)}
    return(values)

def filter_runs_based_on_sample_parameter(IDs, variable_name, variable_value):
    IDs = [k for k in IDs.keys() if IDs[k][variable_name] == variable_value]
    return(IDs)

def prepare_api_requests(ERRs):
    ERRs = [x.split(sep=',') for x in ERRs]
    ERRs = [x for y in ERRs for x in y]
    reqs = ['ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'+x[:6]+'/00'+x[-1]+'/'+x+'/*' for x in ERRs]
    with open('API_REQUESTS.txt', 'a') as f:
        for el in reqs:
            f.write(el + '\n')
    return()

# 'PRJEB13679_metadata.txt'

# get_study_metadata_request('PRJEB13679.txt')
# temp = determine_runs_sample_parameter('PRJEB13679_metadata.txt', 'diagnosis')
# print(show_values_of_variables(temp, 'diagnosis'))
# cd =  filter_runs_based_on_sample_parameter(temp, 'diagnosis', 'CD')
# prepare_api_requests(cd)
# hc =  filter_runs_based_on_sample_parameter(temp, 'diagnosis', 'no')
# prepare_api_requests(hc)

