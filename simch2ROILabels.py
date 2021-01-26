import os
import numpy as np
import glob
from time import time as tm

folder = "/nfs/public/romarco/datasets/20201124"

napas = 6
istep = 800
cstep = 960
apastep = 2*istep + cstep
numch = apastep*napas
tdcs = 6000

def main():
    for fname in glob.glob(f"{folder}/val/evts/p*_simch_*"):
        hits = np.load(fname)
        labels = np.zeros(numch, tdcs)
        assert hits[:,1].min() >= 0
        assert hits[:,1].max() < 15360
        mask = np.logical_and(hits[:,2] >= 0, hits[:,2] < 6000)

        labels[hits[mask][:,1], hits[mask][:,2]] += hits[mask][:,4]
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
    