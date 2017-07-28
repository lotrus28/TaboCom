import xml.etree.ElementTree as xml

def get_api_requests():

    # Parse meta data
    root = xml.parse('../AGP/american_gut_meta').getroot()
    samples =  root.findall("./SAMPLE")

    # Classify all entries based on ibd-diagnosis
    diagnosis = {}
    for i in samples:
        meta = i.findall('SAMPLE_ATTRIBUTES/')
        for j in meta:
            if j[0].text == 'ibd':
                if j[1].text in diagnosis:
                    diagnosis[j[1].text].append(i[0][0].text)
                else:
                    diagnosis[j[1].text] = [i[0][0].text]

    # Classify all entries based on sampling site
    sites = {}
    for i in samples:
        meta = i.findall('SAMPLE_ATTRIBUTES/')
        for j in meta:
            if j[0].text == 'body_site':
                if j[1].text in sites:
                    sites[j[1].text].append(i[0][0].text)
                else:
                    sites[j[1].text] = [i[0][0].text]

    # Filter out only samples from feces and healthy/ill ppl
    sample_ibd_ID = set(sites['UBERON:feces']) & \
             set(diagnosis['Diagnosed by a medical professional (doctor, physician assistant)'])
    sample_healthy_ID = set(sites['UBERON:feces']) & \
             set(diagnosis['I do not have this condition'])

    # Find links to data for all samples
    sample_errs = {}
    for i in samples:
        meta = i.findall('SAMPLE_LINKS/SAMPLE_LINK/XREF_LINK')
        for j in meta:
            if j[0].text == 'ENA-RUN':
                sample_errs[i[0][0].text] = j[1].text

    # Filter out the links we need
    ibd_errs = []
    for el in sample_ibd_ID:
        ibd_errs.append(sample_errs[el])
    healthy_errs = []
    for el in sample_healthy_ID:
        healthy_errs.append(sample_errs[el])

    # make api requests
    ibd_errs = [x.split(sep=',') for x in ibd_errs]
    ibd_errs = sum(ibd_errs, [])
    ibd_errs = ['ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'+x[:6]+'/00'+x[-1]+'/'+x+'/*' for x in ibd_errs]
    healthy_errs = [x.split(sep=',') for x in healthy_errs]
    healthy_errs = sum(healthy_errs, [])
    healthy_errs = ['ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'+x[:6]+'/00'+x[-1]+'/'+x+'/*' for x in healthy_errs]

    # Write all requests to files
    with open('api_ibd.txt','a') as f:
        for el in ibd_errs:
            f.write(el+'\n')
    with open('api_healthy.txt','a') as f:
        for el in healthy_errs:
            f.write(el+'\n')

    return()
	# while read a; do wget "$a"; done < api_ibd.txt
