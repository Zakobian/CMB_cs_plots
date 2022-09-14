import argparse
import numpy as np

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# Training settings
parser = argparse.ArgumentParser(description='CT comparison framework')
parser.add_argument('--cs2', type=float, default=1, metavar='CS2',
                    help='Plot for given cs')
                    # [0.2,0.5,0.7,1]
parser.add_argument('--nts', default=[0.01,0.03,0.1,1], metavar='NT',
                    help='Transition times as percentage of max time')
