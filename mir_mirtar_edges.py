def assign_mir_mir(file, mir_host):
    
    exit()


def get_mir_codes(mirtarfile):
    """

    :param mirtarfile:
    :return:
    """
    o_mirtar = open(mirtarfile)
    first_line = True

    mirtardic = {}

    for line in o_mirtar:

        split_line = line.strip().split('\t')

        if first_line:

            id_dex = split_line.index('miRTarBase ID')
            miR_dex = split_line.index('miRNA')

            first_line = False

            continue

        mirtardic[split_line[miR_dex]] = split_line[id_dex]

    return mirtardic


def import_mir_hosts(mirhost_file, mirtardic):

    """

    :param mirhost_file:
    :param mirtardic:
    :return:

    """

    mirhost_dic = {}

    o_host = open(mirhost_file)

    first_line = True

    for line in o_host:

        split_line = line.strip().split('\t')

        if first_line:

            gene_dex = split_line.index('gene')

            micro_dex = split_line.index('mir')

    return mirhost_dic


def get_file_paths(control_file):

    """

    :param control_file:
    :return control_dic:

    parses control_file into a dictionary providing the code access to all the paths it needs to run

    control file format is @\tkey=path

    key is used as a dictionary key and the path is assigned to that key

    the dictionary is then returned for use by other functions

    """

    o_control = open(control_file)

    control_dic = {}

    for line in o_control:

        if line.startswith('@'):

            split_line = line.strip().split('\t')

            split1 = split_line[1].split('=')

            control_dic[split1[0]] = split1[1]

    return control_dic


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('mirtar_file')
    parser.add_argument('mir_host_file')

    args = parser.parse_args()

    # control_dic = get_file_paths(args.mirtar_file)

    mirtardic = get_mir_codes(args.mirtar_file)
    import_mir_hosts(args.mir_host_file, mirtardic)
