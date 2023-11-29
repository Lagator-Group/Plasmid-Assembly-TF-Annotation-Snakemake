import re
import logging
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# import basic python packages
import numpy as np

from Bio import SeqIO

# import torch packages
import torch
import torch.nn as nn

# import scikit learn packages
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.metrics import f1_score, precision_score, recall_score


def argument_parser(version=None):
    def boolean_string(s):
        if s not in {'False', 'True'}:
            raise ValueError('Not a valid boolean string')
        return s == 'True'
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_dir', required=True, 
                        help='Output directory')
    parser.add_argument('-g', '--gpu', required=False, 
                        default='cpu', help='Specify gpu')
    parser.add_argument('-b', '--batch_size', required=False, type=int,
                        default=128, help='Batch size')
    parser.add_argument('-ckpt', '--checkpoint', required=False, 
                        default='./trained_model/DeepTFactor_ckpt.pt', help='Checkpoint file')
    parser.add_argument('-i', '--seq_file', required=False, 
                        default='./Dataset/ec_seq.fa', help='Sequence data')
    parser.add_argument('-cpu', '--cpu_num', required=False, type=int,
                        default=1, help='Number of cpus to use')    
    return parser