import pandas as pd
from sys import argv
from sys import stdout

# this approach is too slow

'''

map = pd.DataFrame(pd.read_csv(argv[1], sep='\t', header=0, dtype={'score': str, 'distance_abs': str,
                                                                  'distance_act': str, 'num_sites': str}))

new_map = open(argv[1].replace('.', '-collapsed.'), 'w')

old_head = '\t'.join(list(map.columns))
old_head += '\tweight\n'
new_map.write(old_head)

map['join'] = map['regulator_name'] + '$' + map['target_name']

# print(map.head())

nodes_ls = list(set(list(map['join'])))

line_count = 0

for node in nodes_ls:

    check_sub = pd.DataFrame(map.loc[map['join'] == node])

    # print(check_sub.head())

    new_map.write('\t'.join(check_sub.values.tolist()[0][:-1]) + '\t' + str(len(check_sub)) + '\n')

    line_count += 1

    print('{}\r'.format(line_count))
'''


# use dictionary

map_file = open(argv[1])

map_dic = {}

line_count = 0

header = 0

print('\n\n')

for line in map_file:

    line = line.strip().split('\t')

    if line_count == 0:

        header = '\t'.join(line) + '\t' + 'count\n'

        line_count +=1

        continue

    if line[0].startswith('DONT_USE'):
        continue

    # print(line)

    try:

        map_dic[line[0] + '$' + line[4]][1] += 1
        if line[-2] == 'no_data':
            pass
        elif int(line[-3]) < int(map_dic[line[0] + '$' + line[4]][0][-3]):
            map_dic[line[0] + '$' + line[4]][0] = line

    except KeyError:

        map_dic[line[0] + '$' + line[4]] = [line, 1]


    stdout.write('\rLines completed: {}\r'.format(line_count))

    stdout.flush()

    line_count += 1


new_map = open(argv[1].replace('.', '-collapsed.'), 'w')

new_map.write(header)

for key in list(map_dic):

    if map_dic[key][0][2] == 'miRNA':
        print('\t'.join(map_dic[key][0]))

    new_map.write('\t'.join(map_dic[key][0]) + '\t' + str(map_dic[key][1]) + '\n')

stdout.write('\rLines completed: {}\r'.format(line_count))

print('\n\n')