# limits maps to just the TFs which regualte the network

from sys import argv

map_file = open(argv[1])

reg_dic = {}

tf_dic = {}

line_count = 0
head_line = 0

for line in map_file:

    if line_count == 0:

        line_count += 1
        head_line = line
        continue

    split_line = line.strip().split('\t')

    try:

        reg_dic[split_line[4]].append(line)

    except KeyError:

        reg_dic[split_line[4]] = [line]


    tf_dic[split_line[1]] = 0

outfile = open(argv[1].replace('.', '-just_interactions.'), 'w')
outfile.write(head_line)

for key in list(tf_dic):

    try:

        for x in reg_dic[key]:

            outfile.write(x)

    except KeyError:

        continue

map_file.close()
outfile.close()
exit()