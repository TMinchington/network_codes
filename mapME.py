# Uses a file system structure to identify files required and run all of the codes


from sys import argv
import os
import subprocess

experiment_dir = argv[1]

# produce file tree

test_mode = False

subprocess.check_output(['ulimit', '-n', '1000'], shell=True)# remove open file limit on mac os

if os.path.split(experiment_dir)[1] == 'test_peak':

    print('\n\n--------------------------------------------------------------------------------------------'
          '\n\n\t\tRUNNING IN TEST MODE'
          '\n\n--------------------------------------------------------------------------------------------')

    test_mode = True

ex_dir_ls = os.listdir(experiment_dir)

ex_dir_folders = ['work_files', 'map_file', 'input_files', 'log_files']


for folder in ex_dir_folders:

    test_dir = '{}/{}'.format(experiment_dir, folder)

    if not os.path.isdir(test_dir):

        os.makedirs(test_dir)


work_file_folders = ['chip_all_match', 'chip_close_match', 'chip_split', 'genome_split']

work_dir = '{}/work_files'.format(experiment_dir)


for folder in work_file_folders:

    test_dir = '{}/{}'.format(work_dir, folder)

    if not os.path.isdir(test_dir):

        os.makedirs(test_dir)

input_files = ['chip', 'genome', 'host', 'vicious']

input_dir = '{}/input_files'.format(experiment_dir)

input_dir_ls = os.listdir(input_dir)

input_dic = {}


for file in input_files:


    for xfile in input_dir_ls:

        print(file, xfile)

        if file in xfile.lower():

            input_dic[file] = xfile

            break

        else:

            continue


print(input_dic)

for file in input_files:

    if file not in list(input_dic):

        print('{} is missing please move file to input_files directory at : {}\n\n'
              'Note that the file name should contain \"{}\" to be detected'.format(file, input_dir, file))

        exit()

    else:

        continue


genome_file = '{}/{}'.format(input_dir, input_dic['genome'])
chip_file = '{}/{}'.format(input_dir, input_dic['chip'])
pass_this = "{}/peakMULTI.py".format(os.getcwd())

# subprocess.run(["python", "peakMULTI.py", "{}".format(genome_file), "{}".format(chip_file), "{}".format(experiment_dir)])

print('1-------------------------------------------')

subprocess.run(["python", "mir_chip_multi.py", "{}/chip_close_match".format(work_dir), "{}/{}".format(input_dir, input_dic['host']), genome_file, "{}/TFBS_ids.txt".format(work_dir), chip_file])

vicious = '{}/{}'.format(input_dir, input_dic['vicious'])
hosts = '{}/{}'.format(input_dir, input_dic['host'])

print('2-------------------------------------------')

if not test_mode:

    subprocess.run(["python", "mir_mirs_multi.py", vicious, genome_file, "{}/microRNA".format(work_dir), hosts])

    '\nSKIPPING MIR MIRS DUE TO NON COMPATABILITY WITH TEST MODE\n\n'


print('3-------------------------------------------')


subprocess.run(["python", "join_map.py", "{}/map".format(work_dir)])

if test_mode:

    def get_map_file(work_dir):

        map_dir = os.path.split(work_dir)[0]

        map_list = [x for x in os.listdir(map_dir) if 'map-' in x]

        print(map_list)
        print(sorted(map_list))

        return sorted(map_list)[-1]


    map_file = get_map_file(work_dir)


    try:

        out_check = subprocess.check_output(['diff', '{}{}'.format(work_dir.replace('work_files', ''), map_file), '{}/map-compare.tsv'.format(work_dir.replace('work_files', 'control_files'))])

    except subprocess.CalledProcessError as diffex:

        # print(diffex.returncode)
        # print(str(diffex.output))
        print('\n\nTEST FAILED -- MAP DOES NOT MATCH CONTROL MAP\n\n')
        exit()

    print('\nTEST COMPLETED -- NO ERRORS')


