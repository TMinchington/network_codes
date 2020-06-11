"""
compare crms
"""

crm1 = open('/Users/mqbpktm3/Dropbox (The University of Manchester)/0-PhD/1-ChIP Network/17-CRM/input_files/ReMap2_allPeaks-REFORMATED-crm_filtered.bed')
crm2 = open('/Users/mqbpktm3/Dropbox (The University of Manchester)/0-PhD/1-ChIP Network/100-check_network_codes/input_files/remap2018_all_macs2_hg38_v1_2-REFORMATED-crm_filteredNew.bed')
x=0
for line1, line2 in zip(crm1, crm2):
    x+=1
    # print(line1.split()[0:3], line2.split()[0:3])
    if line1.split('\t')[0:3] != line2.split('\t')[0:3]:
        print(x, line1, line2)
        exit()
