 
from argparse import ArgumentParser
import sys

#read command line
parser = ArgumentParser()
parser.add_argument('nbndsub', type=int)
parser.add_argument('bands', type=int)
parser.add_argument('points', type=int)
parser.add_argument('Ef', type=float)
parser.add_argument('bands_file', type=str)
args = parser.parse_args()

nbndsub = args.nbndsub
bands = args.bands
points = args.points
Ef = args.Ef
bands_file = args.bands_file

#read the band energies
e_vec = []

ifs = open(bands_file, "r")
for e in ifs:
    e_vec.append(float(e))

#get the index (nbndskip) of the first band that gets wannierized i.e. crosses the fermi energy (indices starting at 0)
nbndskip = 0
for i in range(0, bands):
    for j in range(0, points):
        index = i*points + j
        if e_vec[index] >= Ef:
            nbndskip = i
            break
    if not nbndskip == 0:
        break

if nbndskip + nbndsub > bands :
    sys.exit("You chose too many bands to wannierize")
    
#get the outer window by finding the smallest and largest values of the specified bands
smallest_outer = e_vec[nbndskip*points]
largest_outer = smallest_outer
for i in range(nbndskip*points, (nbndskip+nbndsub)*points):
    if e_vec[i] < smallest_outer:
        smallest_outer = e_vec[i]
    elif e_vec[i] > largest_outer:
        largest_outer = e_vec[i]

#find the "disturbing" bands that fall into the outer window and do not belong to the specified bands
dist_bands = []
for i in range(0, bands):
    if i >= nbndskip and i < nbndskip + nbndsub:
        continue
    else:
        for j in range(0, points):
            index = i*points + j
            if e_vec[index] > smallest_outer and e_vec[index] < largest_outer:
                dist_bands.append(i)
                break

#find the inner window using the disturbing bands
smallest_inner = smallest_outer;
largest_inner = largest_outer;
for i in range(0, len(dist_bands)):
    smallest_local = e_vec[dist_bands[i]*points]
    largest_local = smallest_local
    for j in range(0, points):
        index = dist_bands[i]*points + j
        #get the largest and smallest values of the current band
        if e_vec[index] < smallest_local:
            smallest_local = e_vec[index]
        elif e_vec[index] > largest_local:
            largest_local = e_vec[index]

    #adjust the inner band using these maximum points of the disturbing band
    if largest_local > smallest_inner and largest_local < largest_inner:
        smallest_inner = largest_local
    if smallest_local < largest_inner and  smallest_local > smallest_inner:
        largest_inner = smallest_local

#band index starts at 1 in q-e making the the nbndskip of this script the proper nbndskip for EPW (i.e. the last unwannierized band)
print("%12.8f    %12.8f          %12.8f          %12.8f          %d" %(smallest_inner, largest_inner, smallest_outer, largest_outer, nbndskip))

