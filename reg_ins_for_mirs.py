"""
regulatory input counts
"""

import os

map_dir = "/Users/mqbpktm3/Dropbox (The University of Manchester)/0-PhD/11-Writing/2019_network_paper/test_host_vs_own/map/"

host_file = "/Users/mqbpktm3/Dropbox (The University of Manchester)/0-PhD/1-ChIP Network/17-CRM/input_files/mirHOST/hsa_hosts.txt"

host_dic = {}

with open(host_file) as open_host:
    for line in open_host:
        split_line = line.strip().split('\t')
        # exit(split_line[6])
        if line.startswith('chr\t') or split_line[6] == 'miRNA':
            continue
        try:
            host_dic[split_line[1]] +=1
        except KeyError:
            host_dic[split_line[1]] = 1

count_dic = {}
total = 0
# print(host_dic)
# exit()
for file in os.listdir(map_dir):
    print(file)

    with open(os.path.join(map_dir, file)) as open_file:
        first_line = True
        for line in open_file:
            if first_line:
                first_line =False
                continue
            split_line = line.strip().split('\t')
            if split_line[-1] == 'na':
                continue
            total += 1
            # print(split_line[8])
            # exit()
            try:

                count_dic[split_line[8]][split_line[-1]] += 1

            except KeyError:

                try:
                    count_dic[split_line[8]][split_line[-1]] = 1

                except KeyError:

                    count_dic[split_line[8]] = {split_line[-1]: 1}

# print(total)
outfile = open('/Users/mqbpktm3/Dropbox (The University of Manchester)/0-PhD/11-Writing/2019_network_paper/test_host_vs_own/mir_coutns.txt', 'w')
outfile.write('mir\thost_type\tcount\tnumber_of_hosts\tintra_inter\n')
for key in count_dic:
    for key2 in count_dic[key]:
        # print(key, key2)
        try:

            print(key, key2, count_dic[key][key2], host_dic[key], 'intragenic')
            outfile.write('\t'.join([str(x) for x in [key, key2, count_dic[key][key2], host_dic[key], 'intragenic']])+'\n')
        except KeyError:
            try: 
                print(key, key2, count_dic[key][key2], 0, 'intragenic', 'intergenic')
                outfile.write('\t'.join([str(x) for x in [key, key2, count_dic[key][key2], 0, 'intergenic']])+'\n')
            except KeyError:
                continue
            