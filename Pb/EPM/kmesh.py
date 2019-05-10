
from argparse import ArgumentParser
import sys
parser = ArgumentParser()
parser.add_argument('grid_size', type=int, nargs=3)
parser.add_argument('weighted', type=int, default=0, nargs='?')
args = parser.parse_args()

n1 = args.grid_size[0]
n2 = args.grid_size[1]
n3 = args.grid_size[2]
weighted = args.weighted

if n1 <= 0:
    sys.exit("n1 needs to be > 0")
if n2 <= 0:
    sys.exit("n2 needs to be > 0")
if n3 <= 0:
    sys.exit("n3 needs to be > 0")
  
totpts = n1 * n2 * n3

if not weighted:
    print("K_POINTS crystal")
    print("%d" % (totpts))
    for i in range(0, n1):
        for j in range(0, n2):
            for k in range(0, n3):
                print("%12.8f%12.8f%12.8f" % (i/n1,j/n2,k/n3))

else:
    print("K_POINTS crystal")
    print("%d" % (totpts))
    for i in range(0, n1):
        for j in range(0, n2):
            for k in range(0, n3):
                print("%12.8f%12.8f%12.8f%14.6e" % (i/n1,j/n2,k/n3,1/totpts))
