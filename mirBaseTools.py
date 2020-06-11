'''-------------------------------tools for looking though mirBase files-------------------------------------------'''

'''------------------------------------------Get microRNA pos data-----------------------------------------------------'''

def mir_pos(species, build='22'):

    from ftplib import FTP

    print('contacting miRBase')

    mirbase = FTP('mirbase.org')
    mirbase.login()
    mirbase.getwelcome()

    mirbase.cwd('/pub/mirbase/')

    mirbase.cwd('{}/genomes/'.format(build))
    # mirbase.retrlines('LIST')
    gen_file = []
    print('collecting data')
    mirbase.retrlines('RETR {}.gff3'.format(species), gen_file.append)
    print('data collected')
    mirbase.quit()

    return gen_file

'''-----------------------------------------------------------------------------------------------------------------'''

def find_all_reads(experiment_file):

    #   this takes the mature read count by experiment and sums the reads

    infile = open(experiment_file)
    outfile = open(experiment_file.replace('.txt', '_summed.txt'), 'w')
    microRNA_ls = []
    read_count_ls = []
    # experiment_ls = []

    for line in infile:

        line = line.strip()
        line = line.split('\t')

        microRNA_ls.append(line[1])
        read_count_ls.append(line[2])
        # experiment_ls.append(line[3])


    unique_mir_ls = list(set(microRNA_ls))

    print(len(unique_mir_ls), len(microRNA_ls))

    read_sum_ls = []

    for miR in unique_mir_ls:

        read_hold = 0

        for x in range(0, len(microRNA_ls)):

            if miR != microRNA_ls[x]:

                continue

            else:

                read_hold += int(read_count_ls[x])

        outfile.write('{}\t{}\n'.format(miR, read_hold))
        # print(miR, read_hold)

#  -------------------------------------------------------------------------------------------------------------------

def mature_info(infile, miRBase_id_col, mature_file):

    # give any tab file with a mirbase MIMAT code column number and it will add species, name etc...

    outfile = open(infile.replace('.txt', '_with-info.txt'), 'w')
    infile = open(infile)
    mature_file = open(mature_file)

    mat_name_ls = []
    mat_name_ns_ls = []
    mimat_ls = []
    species_ls = []
    mir_ls = []
    arm_ls = []
    variant_ls = []

    for line in mature_file:

        if line.startswith('>'):

            line = line.strip()
            line = line.split(' ')

            mat_name_ls.append(line[0].replace('>', ''))
            mimat_ls.append(line[1])
            species_ls.append('{}_{}'.format(line[2], line[3]))

            mat_name_ns_ls.append(line[4])
            mir_temp = (line[4].split('-'))

            if '3p' in mir_temp or '5p' in mir_temp:

                arm_ls.append(mir_temp.pop())

                if len(mir_temp) > 2:

                    variant_ls.append(mir_temp.pop())

                else:

                    variant_ls.append('na')

                mir_ls.append('-'.join(mir_temp))

            else:

                arm_ls.append('na')

                if len(mir_temp) > 2:

                    variant_ls.append(mir_temp.pop())

                else:

                    variant_ls.append('na')

                mir_ls.append('-'.join(mir_temp))

    for line in infile:

        line = line.strip()
        line = line.split('\t')

        mimat = line[miRBase_id_col]
        mimat_dex = mimat_ls.index(mimat)

        outfile.write('{}\t{}\n'.format('\t'.join(line), '\t'.join([mat_name_ls[mimat_dex], species_ls[mimat_dex], mat_name_ns_ls[mimat_dex], mir_ls[mimat_dex], arm_ls[mimat_dex], variant_ls[mimat_dex]])))


#----------------------------------------------------------------------------------------------------------------------

def find_arm_ratio(mir, infile, output_dir):

    import pandas as pd
    #
    # requires output file from find_all_read + mature_info
    # input miroRNA i.e miR-9 checks ratio for each species

    outfile = open('{}\\{}_arm-ratio.txt'.format(output_dir, mir), 'w')
    outfile.write('species\tmicroRNA\t5p_reads\t3p_reads\t%5p\t%3p\n')
    mirframe = pd.DataFrame(pd.read_csv(infile, sep='\t', header=None))

    # print(mirframe)

    mirframe = mirframe.loc[mirframe[5] == mir]

    # print(mirframe)

    species_ls = list(set(list(mirframe[3])))

    for species in species_ls:

        subframe = mirframe.loc[mirframe[3] == species]

        p3 = sum(subframe.loc[subframe[6] == '3p'][1])

        p5 = sum(subframe.loc[subframe[6] == '5p'][1])

        if 'na' in subframe[6]:

            p0 = sum(subframe.loc[subframe[6] == 'na'][1])

        else:

            p0 = 'na'

        if p3 == 0 and p5 == 0:

            print(species, mir, p3, p5, p0, 0, 0)
            outfile.write('{}\n'.format('\t'.join(str(t) for t in [species, mir, p3, p5, p0, 0, 0])))

        else:

            print(species, mir, p3, p5, p0, p5/(p5+p3)*100, p3/(p5+p3)*100)
            outfile.write('{}\n'.format('\t'.join(str(t) for t in [species, mir, p3, p5, p0, p5/(p5+p3)*100, p3/(p5+p3)*100])))

        # print(subframe)

#----------------------------------------------------------------------------------------------------------------------

# find_all_reads('F:\\PhD - Nancy\\mirBase\\Database\\mature_read_count_by_experiment.txt')
# mature_info('F:\\PhD - Nancy\\mirBase\\Database\\mature_read_count_by_experimentsummed.txt', 0, 'F:\\PhD - Nancy\\mirBase\\Database\\mature.fa')

# find_arm_ratio('miR-9', 'F:\\PhD - Nancy\\mirBase\\Database\\mature_read_count_by_experimentsummed_with-info.txt', 'F:\\PhD - Nancy\\mirBase\\Database')


'''---------------------------------------------------------------------------------------------------------------------
This Section is for counting arm reads in miRbase files'''

def get_heads(sql_file):

    import pandas as pd

    infile = open(sql_file)

    table_name_ls = []
    table_heads_ls = []


    table_found = False
    head_set = False
    heads_listed = False
    temp_ls = []


    for line in infile:

        line = line.strip()
        line = line.replace(' ', '')

        if heads_listed and line.startswith('--'):

            head_set = False
            heads_listed = False
            # print(temp_ls)
            table_heads_ls.append('&'.join(temp_ls))
            temp_ls = []

        if head_set:

            if line.startswith('`'):

                heads_listed = True

                splitline = line.split('`')
                # print(splitline[1])
                temp_ls.append(splitline[1])

        if table_found:

            table_found = False

            if '`' in line:

                splitline = line.split('`')
                table_name_ls.append(splitline[1])
                head_set = True

        if line.startswith('--'):

            table_found = True

    table_heads_ls = pd.DataFrame(table_heads_ls, columns=['Heads'])
    table_name_ls = pd.DataFrame(table_name_ls, columns=['Table'])

    outframe = pd.DataFrame(pd.concat([table_name_ls, table_heads_ls],axis=1))

    return(outframe)

# print(get_heads('D:\\MisMatch\\tables.sql'))

def collect_files(direc, mismatch):

    import pandas as pd
    from os import listdir

    mismatch = int(mismatch)

    file_ls = ['mirna_pre_read.txt',
               'mirna.txt',
               'mirna_mature.txt',
               'mirna_pre_mature.txt',
               'experiment_pre_read.txt',
               'tables.sql',
               'mirna_read.txt']

    missing = False

    # make sure all files are present

    for file in file_ls:

        if file not in listdir(direc):
            missing = True
            print('{} is missing from selected directory'.format(file))

    if missing:

        exit()

    '''

    this needs a folder containing:
    - mirna_pre_read.txt
    - mirna.txt
    - mirna_pre_mature.txt
    - experiment_pre_read

    '''

    # get column headers for all files

    table_fame = pd.DataFrame(get_heads('{}/tables.sql'.format(direc)))

    mirna_pre_read_heads = table_fame.loc[table_fame['Table'] == 'mirna_pre_read'].iat[0, 1].split('&')
    mirna_heads = table_fame.loc[table_fame['Table'] == 'mirna'].iat[0, 1].split('&')
    mirna_mature_heads = table_fame.loc[table_fame['Table'] == 'mirna_mature'].iat[0, 1].split('&')
    mirna_pre_mature_heads = table_fame.loc[table_fame['Table'] == 'mirna_pre_mature'].iat[0, 1].split('&')
    experiment_pre_read_heads = table_fame.loc[table_fame['Table'] == 'experiment_pre_read'].iat[0, 1].split('&')
    mirna_read_heads = table_fame.loc[table_fame['Table'] == 'mirna_read'].iat[0, 1].split('&')

    print('mirna_pre_read', mirna_pre_read_heads)
    print('mirna', mirna_heads)
    print('mirna_pre_mautre', mirna_pre_mature_heads)
    print('experiment_pre_read', experiment_pre_read_heads)
    print('mirna_mature', mirna_mature_heads)
    print('mirna_read', mirna_read_heads)

    mirna_pre_read = pd.DataFrame(pd.read_csv('{}/mirna_pre_read.txt'.format(direc), sep='\t',
                                              header=None))

    mirna_read = pd.DataFrame(pd.read_csv('{}/mirna_read.txt'.format(direc), sep='\t',
                                              header=None))

    # test_load = open('{}/mirna_pre_read.txt'.format(direc))

    # for line in test_load:
    #
    #     print(line)


    mirna_pre_read.columns = mirna_pre_read_heads
    # print(mirna_pre_read.head())

    mirna_read.columns = mirna_read_heads

    # print(len(mirna_pre_read))

    mirna_pre_read = mirna_pre_read.loc[(mirna_pre_read['cost_5p'].isin([mismatch, 0])) & (mirna_pre_read['cost_3p'] == 0)]

    # print(len(mirna_pre_read))
    #
    # print(mirna_pre_read)

    mpr_auto_mirna_ls = list(mirna_pre_read['auto_mirna'])

    mirna_heads_drop = mirna_heads[1:]

    # print(mirna_heads)

    mirna_acc_ls = []
    mirna_id_ls = []

    mirna_frame = pd.DataFrame(pd.read_csv('{}/mirna.txt'.format(direc), sep='\t', header=None, index_col=0))
    mirna_frame.columns = mirna_heads_drop
    # print(mirna_frame.head())
    # exit()

    mirna_acc_dex = mirna_heads_drop.index('mirna_acc')
    mirna_id_dex = mirna_heads_drop.index('mirna_id')

    # print(sorted(list(mirna_frame.index), reverse=False))
    # print(sorted(list(mpr_auto_mirna_ls))[:20])
    # exit()

    no_match_ls = []

    '''--------------------------------------------------------------------------------------------------------------'''

    for x in mpr_auto_mirna_ls:

        if x not in mirna_frame.index:

            no_match_ls.append(x)

            mirna_acc_ls.append('no_match')
            mirna_id_ls.append('no_match')

            continue

        mirna_acc_ls.append(mirna_frame.loc[x].iat[mirna_acc_dex])
        mirna_id_ls.append(mirna_frame.loc[x].iat[mirna_id_dex])


    # print(len(mpr_auto_mirna_ls))
    # print(len(mirna_acc_ls))
    # print(len(mirna_id_ls))

    mirna_id_ls = pd.DataFrame(mirna_id_ls, columns=['mirna_id'])
    mirna_acc_ls = pd.DataFrame(mirna_acc_ls, columns=['mirna_acc'])

    mirna_id_ls.index = list(mirna_pre_read.index)
    mirna_acc_ls.index = list(mirna_pre_read.index)

    master_frame = pd.DataFrame(pd.concat([mirna_pre_read, mirna_id_ls, mirna_acc_ls], axis=1))

    # print(master_frame.head())

    # clear out some ram

    mirna_id_ls = 0
    mirna_acc_ls = 0
    mirna_pre_read = 0


    # print(master_frame.loc[master_frame['mirna_acc'] == 'no_match'])

    master_frame = master_frame.loc[master_frame['mirna_acc'] != 'no_match']

    mirna_pre_mature = pd.DataFrame(pd.read_csv('{}/mirna_pre_mature.txt'.format(direc), sep='\t'))
    mirna_pre_mature.columns = mirna_pre_mature_heads

    auto_mature1 = []
    mature_from1 = []
    mature_to1 = []

    auto_mature2 = []
    mature_from2 = []
    mature_to2 = []



    mirna_pre_mature_heads_drop = mirna_pre_mature_heads

    auto_mature_dex = mirna_pre_mature_heads_drop.index('auto_mature')
    mature_from_dex = mirna_pre_mature_heads_drop.index('mature_from')
    mature_to_dex = mirna_pre_mature_heads_drop.index('mature_to')

    mf_auto_read_ls = list(master_frame['auto_mirna'])

    mirna_pre_mature_ls = list(set(list(mirna_pre_mature['auto_mirna'])))

    '''--------------------------------------------------------------------------------------------------------------'''

    for x in mf_auto_read_ls:

        if x not in mirna_pre_mature_ls:

            print('---------------------------')
            print(x)
            print('---------------------------')

            auto_mature1.append('no_match')
            mature_from1.append('no_match')
            mature_to1.append('no_match')

            auto_mature2.append('no_match')
            mature_from2.append('no_match')
            mature_to2.append('no_match')

            continue

        mpm_loc = mirna_pre_mature.loc[mirna_pre_mature['auto_mirna'] == x]

        # print('**************************************************************')
        # print(mpm_loc)
        # print('**************************************************************')
        # print('\n')

        if len(mpm_loc) == 2:

            auto_mature1.append(mpm_loc.iat[0, auto_mature_dex])
            mature_from1.append(mpm_loc.iat[0, mature_from_dex])
            mature_to1.append(mpm_loc.iat[0, mature_to_dex])

            auto_mature2.append(mpm_loc.iat[1, auto_mature_dex])
            mature_from2.append(mpm_loc.iat[1, mature_from_dex])
            mature_to2.append(mpm_loc.iat[1, mature_to_dex])

        elif len(mpm_loc) == 1:

            auto_mature1.append(mpm_loc.iat[0, auto_mature_dex])
            mature_from1.append(mpm_loc.iat[0, mature_from_dex])
            mature_to1.append(mpm_loc.iat[0, mature_to_dex])

            auto_mature2.append('single')
            mature_from2.append('single')
            mature_to2.append('single')

        else:

            # print(mpm_loc)

            auto_mature1.append('something_strange')
            mature_from1.append('something_strange')
            mature_to1.append('something_strange')

            auto_mature2.append('something_strange')
            mature_from2.append('something_strange')
            mature_to2.append('something_strange')



    # print(len(auto_mature1))
    # print(len(mature_from1))
    # print(len(mature_to1))

    auto_mature1 = pd.DataFrame(auto_mature1, columns=['auto_mature1'])
    mature_from1 = pd.DataFrame(mature_from1, columns=['mature_from1'])
    mature_to1 = pd.DataFrame(mature_to1, columns=['mature_to1'])

    auto_mature1.index = list(master_frame.index)
    mature_from1.index = list(master_frame.index)
    mature_to1.index = list(master_frame.index)

    auto_mature2 = pd.DataFrame(auto_mature2, columns=['auto_mature2'])
    mature_from2 = pd.DataFrame(mature_from2, columns=['mature_from2'])
    mature_to2 = pd.DataFrame(mature_to2, columns=['mature_to2'])

    auto_mature2.index = list(master_frame.index)
    mature_from2.index = list(master_frame.index)
    mature_to2.index = list(master_frame.index)

    master_frame = pd.DataFrame(pd.concat([master_frame, auto_mature1, mature_from1,
                                           mature_to1, auto_mature2, mature_from2, mature_to2], axis=1))

    # print(master_frame.head())

    '''--------------------------------------------------------------------------------------------------------------'''

    #   collect and sum reads

    epr_df = pd.DataFrame(pd.read_csv('{}/experiment_pre_read.txt'.format(direc), sep='\t'))
    epr_df.columns = experiment_pre_read_heads

    mf_auto_read_ls = list(master_frame['auto_read'])

    read_sum_ls = []

    print('HEY!!!!')

    for y in mf_auto_read_ls:   # can probably drop reads here by filter for auto read length. Just need to incorporate

        # print(y)`

        read_mirna_sub = mirna_read.loc[mirna_read['auto_read'] == y]

        if len(list(read_mirna_sub['sequence'])[0]) < 17:

            print('short: ', y)

            read_sum_ls.append(0)

            continue


        read_sub = epr_df.loc[epr_df['auto_read'] == y]
        read_count_sum = sum(list(read_sub['count']))

        read_sum_ls.append(read_count_sum)

    read_sum_ls = pd.DataFrame(read_sum_ls, columns=['sum_count'])
    read_sum_ls.index = list(master_frame.index)

    master_frame = pd.DataFrame(pd.concat([master_frame, read_sum_ls], axis=1))

    print(master_frame)

    pd.DataFrame.to_csv(master_frame, '{}/miRBaseREADTABLE.txt'.format(direc), sep='\t', index=False)



# collect_files('D:\\MisMatch', 1)

def map_and_count(direc, collect_files_file, mismatch):

    import pandas as pd

    infile = open(collect_files_file)

    arm_ls = []
    match_ls = []

    overhang = 10

    for line in infile:

        if line.startswith('auto_read'):

            continue

        line = line.strip()
        line = line.split()

        if line[9] == 'something_strange':

            arm_ls.append('no_match')
            match_ls.append('no_match')

        elif int(line[2]) in range(int(line[9])-overhang, int(line[10])-overhang) and int(line[3]) == 0:

            arm_ls.append('5\'')
            match_ls.append(0)

        elif int(line[2]) in range(int(line[9])-overhang, int(line[10])-overhang) and int(line[3]) == 1:

                arm_ls.append('5\'')
                match_ls.append(1)

        elif line[12] != 'single':

            if int(line[2]) in range(int(line[12])-overhang, int(line[13])-overhang) and int(line[3]) == 0:

                arm_ls.append('3\'')
                match_ls.append(0)

            elif int(line[2]) in range(int(line[12])-overhang, int(line[13])-overhang) and int(line[3] == 1):

                arm_ls.append('3\'')
                match_ls.append(1)

            else:

                arm_ls.append('no_match')
                match_ls.append('no_match')
                # print(' '.join(line))

                continue

        else:

            arm_ls.append('no_match')
            match_ls.append('no_match')

            # print(' '.join(line))

            continue

    infile.close()

    master_frame = pd.DataFrame(pd.read_csv(collect_files_file, sep='\t', header=0, low_memory=False))

    arm_ls = pd.DataFrame(arm_ls, columns=['arm'])
    match_ls = pd.DataFrame(match_ls, columns=['mismatch'])

    mirna_mature = pd.DataFrame(pd.read_csv('{}/mirna_mature.txt'.format(direc), header=None, sep='\t', index_col=0))

    auto_mature1 = list(master_frame['auto_mature1'])
    auto_mature2 = list(master_frame['auto_mature2'])

    mimat1 = []
    mimat2 = []

    for auto_mat in auto_mature1:

        if auto_mat == 'something_strange' or auto_mat == 'single':

            mimat1.append('no_match')

            continue

        mimat1.append(mirna_mature.loc[int(auto_mat)].iat[2])

    for auto_mat in auto_mature2:

        if auto_mat == 'something_strange' or auto_mat == 'single':

            mimat2.append('no_match')

            continue

        mimat2.append(mirna_mature.loc[int(auto_mat)].iat[2])


    mimat1 = pd.DataFrame(mimat1, columns=['mature_acc_5\''])
    mimat2 = pd.DataFrame(mimat2, columns=['mature_acc_3\''])

    master_frame = pd.DataFrame(pd.concat([master_frame, arm_ls, match_ls, mimat1, mimat2], axis=1))

    # print(master_frame)
    print(len(master_frame))
    print(len(master_frame.loc[master_frame['mismatch'] == 0]))
    print(len(master_frame.loc[master_frame['mismatch'] == 1]))

    master_frame = master_frame.loc[master_frame['mismatch'] != 'no_match']

    if mismatch == 0:

        master_frame = master_frame.loc[master_frame['mismatch'] == 0]

    print(master_frame)

    mirna_id_ls = list(set(list(master_frame['mirna_id'])))

    mirna_acc_ls = []
    mirna_5_count = []
    mirna_3_count = []
    mature_5_acc = []
    mature_3_acc = []

    species = []

    for mirna in mirna_id_ls:

        mirna_sub = master_frame.loc[master_frame['mirna_id'] == mirna]
        mirna_acc_ls.append(mirna_sub.iat[0, 7])
        species.append(mirna_sub.iat[0, 6].split('-')[0])

        mirna_sub5 = mirna_sub.loc[mirna_sub['arm'] == '5\'']
        mirna_sub3 = mirna_sub.loc[mirna_sub['arm'] == '3\'']

        mirna_5_count.append(sum(int(x) for x in list(mirna_sub5['sum_count'])))
        mirna_3_count.append(sum(int(x) for x in list(mirna_sub3['sum_count'])))

        if len(mirna_sub5) > 0:

            mature_5_acc.append(mirna_sub5.iat[0, 17])

        else:

             mature_5_acc.append('no_match')

        if len(mirna_sub3) > 0:

            mature_3_acc.append(mirna_sub3.iat[0, 18])

        else:
            mature_3_acc.append('no_match')

    mirna_id_ls = pd.DataFrame(mirna_id_ls, columns=['mirna_id'])
    mirna_acc_ls = pd.DataFrame(mirna_acc_ls, columns=['mirna_acc'])
    mirna_5_count = pd.DataFrame(mirna_5_count, columns=['5_arm'])
    mirna_3_count = pd.DataFrame(mirna_3_count, columns=['3_arm'])

    mature_5_acc = pd.DataFrame(mature_5_acc, columns=['mature_acc_5'])
    mature_3_acc = pd.DataFrame(mature_3_acc, columns=['mature_acc_3'])
    species = pd.DataFrame(species, columns=['species'])

    outframe = pd.DataFrame(pd.concat([species, mirna_id_ls, mirna_acc_ls, mature_5_acc, mirna_5_count,
                                       mature_3_acc, mirna_3_count], axis=1))

    pd.DataFrame.to_csv(outframe, '{}/miRbase_assigned_reads_mis{}.txt'.format(direc, mismatch), sep = '\t', index=None)

    print(outframe)

    #make a MIMAT_count

    mimat_ls = list(set(list(outframe['mature_acc_5'])+list(['mature_acc_3'])))
    reads_ls = []

    for mimat in mimat_ls:

        mat5 = outframe.loc[outframe['mature_acc_5'] == mimat]
        mat3 = outframe.loc[outframe['mature_acc_3'] == mimat]

        reads_ls.append(sum(list(mat5['5_arm']))+sum(list(mat3['3_arm'])))

    mimat_ls = pd.DataFrame(mimat_ls, columns=['MIMAT'])
    reads_ls = pd.DataFrame(reads_ls, columns=['read_count'])

    outframe2 = pd.DataFrame(pd.concat([mimat_ls, reads_ls], axis=1))

    pd.DataFrame.to_csv(outframe2, '{}/MIMATreads_mis{}.txt'.format(direc, mismatch), sep='\t', index=None)



    print(outframe2)



# map_and_count('D:\\MisMatch\\DATA\\miRBaseREADTABLE.txt', 0)

'''-----------------------------------------------------------------------------------------------------------------'''

def count_mature(gff3):

    primary = 0
    arms = 0

    gff3 = open(gff3)

    for line in gff3:

        if line.startswith('#'):

            continue

        if "miRNA_primary_transcript" in line:

            primary+=1

        elif "MIMAT" in line:

            arms+=1

        else:

            continue

    print('pri: ', primary)
    print('arm: ', arms)

# count_mature('F:/PhD - Nancy/mirBase/hsa_miRs.gff3')

def fasta_get_species(species, file):


    fasta = open(file)
    outfile = open(file.replace('.fa', '-{}.fa'.format(species)), 'w')
    match = False
    for line in fasta:

        if match:

            outfile.write(line)
            match = False
            continue


        if species in line:

            outfile.write(line)
            match = True
            continue

        else:

            continue

# fasta_get_species('hsa', '/Users/mqbpktm3/Dropbox (The University of Manchester)/0-PhD/1-ChIP Network/2-REMAP/7-miRtar/input_files/mature.fa')

def get_list(fasta_file):

    # goes through fasta file such as high confidence microRNAs and gets a list of the microRNAs

    file = open(fasta_file)

    outfile = open(fasta_file.replace('.fa', '_listed.txt'),'w')

    outfile.write('microRNA_name\tmirna\n')

    for line in file:

        if line.startswith('>'):

            line = line.strip().split()

            outfile.write('{}\t{}\n'.format(line[0].replace('>', ''), line[1]))

        else:

            continue

# get_list('/Users/mqbpktm3/Dropbox (The University of Manchester)/0-PhD/1-ChIP Network/2-REMAP/7-miRtar/input_files/mature-hsa.fa')
