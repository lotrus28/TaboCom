= ena_download_agp.py
A script to create ENA API download requests for sequenator runs on IBD or healthy samples.
AGP metadata required
To actually download fastq-files use:
while read a; do wget "$a"; done < api_ibd.txt

= parse_metadata.py
More generic script that

= api_ibd.txt
= api_heal.txt
These file contain a list of ENA download requests

= american_gut_meta.7z
Arrchived metadata on AGP runs published in ENA as of 24.11.2016
Full size: 736MB
