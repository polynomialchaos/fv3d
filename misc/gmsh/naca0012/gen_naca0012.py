####################################################################################################################################
# Imports and definitions
#-----------------------------------------------------------------------------------------------------------------------------------
import numpy as np

fileName = 'xy_naca0012.geo'
fileNameXY = 'xy_naca0012.dat'

####################################################################################################################################
# Main function
#-----------------------------------------------------------------------------------------------------------------------------------
with open(fileName,'w') as f:
    # definitions
    f.write('lc_naca = 0.005;\n')

    # write the xyz points of naca profile
    with open(fileNameXY,'r') as xyf:
        for idx,line in enumerate(xyf.readlines()):
            xyData = [float(x) for x in line.split()]
            f.write('Point({:}) = {{{:},{:},{:},{:}}};\n'.format(idx,xyData[0],xyData[1],0.0,'lc_naca'))
