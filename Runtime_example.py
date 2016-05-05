import makenc
import sys
# cliff change this to where the FRF_natneighbor.py lives if not in same folder
sys.path.append('/home/spike/repos/frf_bathy_interp')
import FRF_natneighbor as nn
import getopt
# look for files = flist

for fname in filelist:
    gridDict = nn.frf_grid_product(fname, dxdy=10)
    ofname = fname[:-4]+'_grid.txt'
    nn.write_grid(ofname=ofname, grid_dict=gridDict)
    # put plotting functing that saves file here

if __name__ == "__main__":
    opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
    #.
    #
    fname = args[0]