from sys import argv
import os
import pandas as pd

map_dir = argv[1]

out_dir = (map_dir.split('/')[:-2])

print('/'.join(out_dir))
out_dir = '/'.join(out_dir)

file_name = '{}/map-0.tsv'.format(out_dir)

def do_not_overwrite(file_name):

    import os

    while os.path.isfile(file_name):
        print('Found match')

        file = os.path.split(file_name)[1]
        direct = os.path.split(file_name)[0]

        current = int(file.split('-')[1].split('.')[0])

        file = file.replace(str(current), str(current + 1))

        file_name = direct + '/' + file

    return file_name

file_name = do_not_overwrite(file_name)


map_file = open(file_name, 'w')
map_file.write('regulator\tregulator_name\treg_type\ttarget\ttarget_name\ttarget_type\tscore\tdistance_act\tdistance_abs\tnum_sites\n')
map_files = os.listdir(map_dir)
chip_files = [x for x in map_files if not x.startswith('.') and 'ccm' in x]
mir_files = [x for x in map_files if not x.startswith('.') and 'mir_temp' in x]

for file in chip_files:
    print(file)

    work_file = open('{}/{}'.format(map_dir, file))

    for line in work_file:



        line = line.strip()
        line = line.split()

        if line[3] == 'chip_gene_id' or line[3] == 'no_match':

            continue
        # print(line[0])
        map_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line[3], line[2], line[4], line[6], line[8], line[10], line[14], line[12], line[11], 'no_data'))

    work_file.close()




for file in mir_files:

    work_file = open('{}/{}'.format(map_dir, file))
    print(file)
    for line in work_file:

        line = line.strip()
        line = line.split()

        if line[0] == 'regulator' or line[3] == 'no_match':

            continue

        map_file.write("{}\t{}\tmiRNA\t{}\t{}\t{}\t{}\tno_data\tno_data\tno_data\n".format(line[2], line[0], line[7], line[6], line[9], line[4], line[9]))

    work_file.close()

map_file.close()



map = pd.read_csv(file_name, sep='\t', header=0)

print(map.head())

print(len(map))
map2 = map.drop_duplicates(inplace=False)
print(map.columns)
print(len(map2))

pd.DataFrame.to_csv(map2, file_name,  sep='\t', index=False)
