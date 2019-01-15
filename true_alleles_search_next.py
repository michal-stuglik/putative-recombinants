#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# License: GPL v3

from Bio import SeqIO
import os
import sys
import argparse

seqs = {}

diff_threshold = 1


def parse_fasta_to_dic(input_file_hlr):
    seq_dic = dict()

    for sequence in SeqIO.parse(input_file_hlr, "fasta"):
        seq_dic[sequence.id] = str(sequence.seq)

    input_file_hlr.close()
    return seq_dic


def recombinant(seq_dic: dict) -> dict:

    target_set = set()

    recombinant_dict = dict()
    recombinant_counter = 0

    for seq_id, seq in sorted(seq_dic.items()):
        recombinant_counter += 1

        log_str = 'searching recombinant for {}'.format(seq_id)
        print(log_str)

        seq = str(seq)

        # i_found_it = False
        for i in range(len(seq)):

            sub_seq_l = seq[0:i + 1]
            sub_seq_r = seq[i + 1:len(seq)]

            for id_parent_l, seq_parent_l in seq_dic.items():
                seq_parent_l = str(seq_parent_l)

                if seq == seq_parent_l:
                    continue

                if str(seq_parent_l).startswith(sub_seq_l):
                    for id_parent_r, seq_parent_r in seq_dic.items():

                        key = "{}_{}_{}".format(seq_id, id_parent_l, id_parent_r)

                        if key in recombinant_dict:
                            continue

                        seq_parent_r = str(seq_parent_r)

                        if seq == seq_parent_r:
                            continue

                        if seq == seq_parent_l:
                            continue

                        if seq_parent_l == seq_parent_r:
                            continue

                        if str(seq_parent_r).endswith(sub_seq_r):
                            log_str = 'putative parents were discovered for allele {}'.format(seq_id)
                            print(log_str)

                            print("target: {} | sources: {}, {}".format(seq_id, id_parent_l, id_parent_r))
                            recombinant_dict[key] = [seq_id, id_parent_l, id_parent_r]

                            target_set.add(seq_id)
                            # set_chimera_allele_for_individual(ind.idx, seq_id)
                            i_found_it = True

                        # if i_found_it:
                        #     break

                # if i_found_it:
                #     break

            # if i_found_it:
            #     break

        if seq_id not in target_set:
            recombinant_dict[seq_id] = [seq_id, 0, 0]

    return recombinant_dict


def parsing_fasta(input_fasta) -> dict:
    fasta_dict = dict()
    for sequence in SeqIO.parse(open(input_fasta, 'r'), "fasta"):
        fasta_dict[sequence.id] = sequence.seq
    return fasta_dict


def write_output(file_name: str, recombinant_dict: dict):
    with open(file_name, 'w') as f:
        f.write("{}\t{}\t{}\n".format("target", "source1", "source2"))
        for  rec in recombinant_dict.values():
            f.write("{}\t{}\t{}\n".format(rec[0], rec[1], rec[2]))


def main(args: list):
    parser = argparse.ArgumentParser(
        description='Program searches putative chimeric sequences in specific range of alleles frequency')

    parser.add_argument("-f", "--input_fasta", dest="INPUT_FASTA",
                        help="path to input fasta file with allel sequences", action='store')

    parser.add_argument("-o", "--output", dest="OUTPUT",
                        help="output file, no chimeras, no substitutions", action='store')

    options = parser.parse_args(args)
    print("\nEntering program: " + parser.prog + \
          '\nDesc: ' + parser.description)

    os.chdir(os.curdir)

    print('Parsing fasta file...')
    fasta_dict = parsing_fasta(options.INPUT_FASTA)

    print('Checking chimeras...')
    recombinant_dict = recombinant(seq_dic=fasta_dict)

    write_output(options.OUTPUT, recombinant_dict)

    print("Leaving program: " + parser.prog)


if __name__ == "__main__":
    main(sys.argv[1:])
