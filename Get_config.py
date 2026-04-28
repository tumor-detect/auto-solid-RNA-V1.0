#!/public/software/anaconda3/envs/python27/bin/python
# -*- coding: utf-8 -*-

import configparser
import os
import io


script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def get_disease_dict():
    disease_file_path = os.path.join(script_dir, "pipeline", "disease_info.txt")
    with open(disease_file_path, "rt") as inf:
        disease_dict = {}
        lines = [line.strip().split("\t") for line in inf.readlines()]
        for line in lines:
            disease_dict[line[0]] = line[1]
    return disease_dict


