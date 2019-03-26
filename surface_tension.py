#!/usr/bin/env python
"""surface_tension

 calculate surface tension from pressure tensor
 							                                    
 -Andrew P. Santos					                            
 
"""
import sys, argparse
import numpy as np
import math as ma        
from scipy import stats
import scipy.special
from scipy.optimize import curve_fit

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

class SurfaceTension(object):
    """Callable SurfaceTension class

       Reads pressure tensor data
       Calculates and writes surface tension
    """
    def __init__(self):
        self.sqrtpi = ma.sqrt( ma.pi ) # sqrt of pi

    def addParser(self, parser):
        """
        Get relevant values from an argparse parser
        """
        if (parser.parse_args().start_line != None):
            self.start_line = parser.parse_args().start_line
        else:
            print 'assuming startline is 1'
            self.start_line = 1

        self.ifilename = parser.parse_args().filename
        if parser.parse_args().outfilename:
            self.ofilename = parser.parse_args().outfilename
        else:
            self.ofilename = 'surface_tension.dat'

        if parser.parse_args().perp_direction == 'x':
            self.direction = 0
        elif parser.parse_args().perp_direction == 'y':
            self.direction = 1
        elif parser.parse_args().perp_direction == 'z':
            self.direction = 2

        if parser.parse_args().L_perp:
            self.L_perp = parser.parse_args().L_perp
            self.half_L_perp = parser.parse_args().L_perp / 2.0

    def readPressureTensor(self, filename=None):
        """
        Read density profile
        """
        if filename == None:
            filename = self.ifilename

        try:
            ifile = open(filename, 'r')
        except IOError:
            raise IOError('cannot find: %s. Either put the file in the working directory or fix the input file.' % filename)

        ntime = file_len(filename) - self.start_line

        self.pressure_tensor = np.zeros( (3,3,ntime), dtype=np.float )

        i_line = -1
        data_line = 0
        for line in ifile:
            i_line += 1

            if (i_line < self.start_line):
                continue

            data = line.strip().split()
            if len(data) < 7:
                raise IOError('Not enough collumns in %s based on input information.' % filename)
        
            self.pressure_tensor[0, 0, data_line] = float( data[1] )    # xx
            self.pressure_tensor[1, 1, data_line] = float( data[2] )    # yy
            self.pressure_tensor[2, 2, data_line] = float( data[3] )    # zz
            self.pressure_tensor[0, 1, data_line] = float( data[4] )    # xy
            self.pressure_tensor[0, 2, data_line] = float( data[5] )    # xz
            self.pressure_tensor[1, 2, data_line] = float( data[6] )    # yz

            data_line += 1

        ifile.close()

    def calculateSurfaceTension(self, filename=None):
        """
        Calculate the surface tension
        """
        P_per_mean = np.mean(self.pressure_tensor[self.direction,self.direction,:])
        para_directions = [0,1,2]
        para_directions.pop(self.direction)
        P_para_mean = np.mean(self.pressure_tensor[para_directions[0],para_directions[0],:])
        P_para2_mean = np.mean(self.pressure_tensor[para_directions[1],para_directions[1],:])
        self.surface_tension = self.half_L_perp * (P_per_mean - 0.5*(P_para_mean + P_para2_mean))

    def write(self, filename=None):
        """
        Write the surface tension
        """
        if filename == None:
            filename = self.ofilename

        ofile = open(filename, 'w')

        ofile.write("%10.8f\n" % ( self.surface_tension ) )
        ofile.close()

def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Fit a density profile and extract valuble information')
    parser.add_argument('-f', "--filename", type=str, required=True,
                   help='file with pressure tensor data. Must be in (LAMMPS) form: "time P_xx P_yy P_zz P_xy P_xz Pyz"')
    parser.add_argument("-o", "--outfilename", type=str, default='surface_tension.dat',
                   help='output file name, assumed to be surface_tension.dat')
    parser.add_argument('-l', "--start_line", type=int,
                   help='Set the starting line number.')
    parser.add_argument("--perp_direction", type=str, default='z', choices=['x','y','z'],
                   help='Dimension perpendicular to the surface.')
    parser.add_argument("--L_perp", type=float, required=True,
                   help='box length perpendicular to the interface')

    # Initialize the Activity class
    surface_tension = SurfaceTension()

    # Tell the class everything specified in files and command line
    err = surface_tension.addParser(parser)

    if err != None:
        print 'error parsing info'
        return

    # read in data 
    surface_tension.readPressureTensor()
    calc_error = surface_tension.calculateSurfaceTension()

    if not calc_error:
        surface_tension.write()

if __name__ == '__main__':
    sys.exit(main())

# "The greatest happiness is to know the source of unhappiness." -Fyodor Dostoevsky
