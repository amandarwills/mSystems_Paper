#!/usr/bin/env python3

import re

fin = open("Stress_Gene_List.csv")
for line in fin:
   # line = re.sub(r"-", r" ", line)
    split_line = line.split()
    fother = open("MC_T5_FC_protvgene_correlation_results_Annotation.csv")
    for other_line in fother:
        column_f = other_line.split(",")[5]
        split_other_line = column_f.split()

        num_matches = 0
        for word in set(split_line):
            for other_word in set(split_other_line):
                if word.lower() == other_word.lower():
                    num_matches += 1
            

        if num_matches >= min(3, len(split_other_line)):
            print(line.rstrip())
            print(column_f.rstrip())
            print(other_line.rstrip())
            print()

    fother.close()
