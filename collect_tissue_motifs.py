"""
collect_tissues
"""

import os

tissue_path = "/Users/mqbpktm3/Dropbox (The University of Manchester)/0-PhD/1-ChIP Network/check_numbers/tissues"

tissue_list = [x for x in os.listdir(tissue_path) if '.' not in x and x not in ['motif_master']]

print(tissue_list)

motif_master_path = os.path.join(tissue_path, 'motif_master')

if not os.path.isdir(motif_master_path):
    os.makedirs(motif_master_path)

first_pass = True

for tissue in tissue_list:
    tiss_mo_path = os.path.join(tissue_path, f'{tissue}/fast_motifs')
    print(tiss_mo_path)

    mo_files = [x for x in os.listdir(tiss_mo_path) if '.txt']

    if first_pass:
        first_pass = False
        mo_dic = {}
        for mo in mo_files:
            mo_dic[mo] = open(os.path.join(motif_master_path, mo.replace('.t', '-master.t')), 'w')

    for mo in mo_files:

        with open(os.path.join(tiss_mo_path, mo)) as open_mo:
            for line in open_mo:
                mo_dic[mo].write(f'{tissue}\t{line}')