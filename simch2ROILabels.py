import os
import argparse
import numpy as np
import glob
from time import time as tm
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("--path", "-p", default="/nfs/public/romarco/datasets/20201124/val/evts",
                    type=str, help='Dataset directory path to convert')
args = vars(parser.parse_args())

napas = 6
istep = 800
cstep = 960
apastep = 2*istep + cstep
numch = apastep*napas
tdcs = 6000

def main():
    fnames = os.path.join(args["path"], "p*_simch_ev*")
    for fname in tqdm(glob.glob(fnames)):
        hits = np.load(fname)
        labels = np.zeros([numch, tdcs])
        assert hits[:,1].min() >= 0
        assert hits[:,1].max() < 15360
        mask = np.logical_and(hits[:,2] >= 0, hits[:,2] < 6000)

        masked = hits[mask]
        labels[masked[:,1].astype(int), masked[:,2].astype(int)] += hits[mask][:,4]
        nevt = np.ones_like(labels[:,:1])*hits[0,0] # event number
        chs = np.arange(numch)[:,None] # channels
        labels = np.concatenate([nevt, chs, labels], axis=1)

        sname = fname.split("_")
        sname.insert(-1, "labels")
        sname = "_".join(sname)
        np.save(sname, labels)


if __name__ == "__main__":
    start = tm()
    main()
    print(f"Program done in {tm()-start}")

