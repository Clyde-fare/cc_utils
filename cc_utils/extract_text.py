#! /usr/bin/env python
__author__ = 'clyde'

import sys


def extract_txt(pattern_1, pattern_2, in_file, out_file):
    """extracts all text found in in_file that is between pattern_1 and pattern_2 and writes it to out_file"""

    with open(in_file) as inp_f:
        in_section = False

        section = []
        sub_section = []

        for line in inp_f:
            if pattern_1 in line:
                in_section = True

            elif pattern_2 in line and in_section:
                section += sub_section

                in_section = False
                sub_section = []

            elif in_section:
                sub_section.append(line)


    with open(out_file, 'w') as outp_f:
        outp_f.write('\n'.join(section))


if __name__ == '__main__':
    extract_txt(*sys.argv[1:])