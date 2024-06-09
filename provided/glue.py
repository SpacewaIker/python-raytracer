# Thibaut Baguette
# 261001513
# COMP557/ECSE532 Assignment 4
# 2024-04-08

import numpy as np
import argparse
from PIL import Image
import os

parse = argparse.ArgumentParser()
parse.add_argument("--outfile", type=str, default="out.png")
parse.add_argument("--dir", type=str)
args = parse.parse_args()

# get all files in the directory
files = os.listdir(args.dir)
files.sort(key=lambda x: int(x.split(".")[0]))
arrays = [np.load(args.dir + "/" + file) for file in files]

# concatenade all arrays
image = np.concatenate(arrays, axis=0)
image = np.rot90(image, k=1, axes=(0, 1))

im = Image.fromarray((image * 255).astype(np.uint8))
im.save(args.outfile)
im.show()
