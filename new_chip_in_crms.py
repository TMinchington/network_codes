'''chip_in_crms2'''

def check_crm_file(crm):
    crm_ls = []
    current_pos = 0
    current_chr = 'lalala'
    with open(crm) as open_crm:
        
        for line in open_crm:
            
            split_line = line.strip().split('\t')
            chromosome, start, stop = split_line[0:3]
            start, stop = int(start), int(stop)
            if current_chr != chromosome:
                if chromosome in crm_ls:
                    exit(f'File out of order CHROMOSOME {line}')

                crm_ls.append(chromosome)
                current_chr = chromosome
                current_pos = 0
            
            if current_pos > start:
                exit(exit(f'File out of order START_pos {line}'))

            current_pos = start
    return crm_ls


def find_tfbs_in_crm(crm, tfbs):
    otfbs = open(tfbs)
    outfile = open(tfbs.replace('.', '-crm_filteredNew.'),'w')
    tf_pos = -100
    with open(crm) as ocrm:

        for line in ocrm:
            crm_split_line = line.strip().split('\t')
            current_crm_Chr, current_crm_start, current_crm_stop = crm_split_line[0:3]
            current_crm_start, current_crm_stop = int(current_crm_start), int(current_crm_stop)

            if  current_crm_start <= tf_pos <= current_crm_stop:
                    # print('in', current_crm_Chr, current_crm_start, current_crm_stop, tf_chr, tf_pos)
                outfile.write(line2)

            for line2 in otfbs:
                split_line2 = line2.strip().split('\t')
                tf_chr, tf_pos = split_line2[0], int(split_line2[6])
                # print(tf_pos, current_crm_start, current_crm_stop)
                if tf_chr != current_crm_Chr:
                    # print(tf_chr, current_crm_Chr)
                    break

                elif  current_crm_start <= tf_pos <= current_crm_stop:
                    # print('in', current_crm_Chr, current_crm_start, current_crm_stop, tf_chr, tf_pos)
                    outfile.write(line2)
                    continue

                elif tf_pos > current_crm_stop:
                    # print(tf_pos, current_crm_stop)
                    break

                else:
                    continue
                

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('crm')
    parser.add_argument('tfbs')

    args = parser.parse_args()

    # crm_chr = check_crm_file(args.crm)
    # tfbs_chr = check_crm_file(args.tfbs)

    # print(crm_chr)
    # print(tfbs_chr)
    # print(crm_chr == tfbs_chr)
    find_tfbs_in_crm(args.crm, args.tfbs)



