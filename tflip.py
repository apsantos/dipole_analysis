#!/usr/bin/env python
"""dipole_orientation

 Calculate the dipole-orientaiton parameter, cos(theta)
 							                                    
 -Andrew P. Santos					                            
 
"""
import sys, argparse
import numpy as np
import math as ma        
from scipy import stats
import scipy.special
from scipy.optimize import curve_fit

class DipoleOrientation(object):
    """Callable Density Profile class

       Reads  density profile and fits
    """
    def __init__(self):
        self.tflip = []

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
            self.ofilename = 'density_profile_fit.dat'

        if parser.parse_args().t_collumn >= 0:
            self.ic_t = parser.parse_args().t_collumn
            self.t_conversion = parser.parse_args().t_conversion

        if parser.parse_args().dipole_collumn >= 0:
            self.ic_dmag = parser.parse_args().dipole_collumn

        if parser.parse_args().field_collumn >= 0:
            self.ic_dfield = parser.parse_args().field_collumn

    def readDipoleProfile(self, filename=None):
        """
        Read the dipole profile
        """
        if filename == None:
            filename = self.ifilename

        try:
            ifile = open(filename, 'r')
        except IOError:
            raise IOError('cannot find: %s. Either put the file in the working directory or fix the input file.' % filename)

        self.z = []
        self.cosTheta = []
        self.P2 = []

        max_ic = max( self.ic_t, self.ic_dfield, self.ic_dmag )
        i_line = 0
        previous_dip_sign = 1.0
        for line in ifile:
            i_line += 1
            if (i_line < self.start_line):
                continue

            data = line.strip().split()
            if len(data) <= max_ic:
                raise IOError('Not enough collumns in %s based on input information.' % filename)

            if np.sign( previous_dip_sign * float( data[self.ic_dmag] ) ) == np.sign( float( data[self.ic_dfield] ) ):
                self.tflip.append( float( data[self.ic_t] ) * self.t_conversion )
                previous_dip_sign *= -1.0

        ifile.close()

    def write(self, filename=None):
        """
        Write the fit parameters and the fit values
        """
        if filename == None:
            filename = self.ofilename

        ofile = open(filename, 'w')


        ofile.write('# flip times\n')

        for i in range(len(self.tflip)):
            ofile.write("%10.8f 0\n" % ( self.tflip[i] ) )

        ofile.close()

def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Fit a density profile and extract valuble information')
    parser.add_argument('-f', "--filename", type=str, required=True,
                   help='file with density data.')
    parser.add_argument("-o", "--outfilename", type=str, default='dipole_flip_times.dat',
                   help='output file name, assumed to be dipole_flip_times.dat')
    parser.add_argument('-l', "--start_line", type=int,
                   help='Set the starting line number.')
    parser.add_argument("--t_collumn", type=int, 
                   help='Collumn with the simulation time')
    parser.add_argument("--t_conversion", type=float, default=1.0,
                   help='convert the times by this conversion')
    parser.add_argument("--field_collumn", type=int, 
                   help='Collumn with the total dipole magnitude')
    parser.add_argument("--dipole_collumn", type=int, 
                   help='Collumn with the dipole magnitude in the applied field direction')

    # Initialize the Activity class
    dipole = DipoleOrientation()

    # Tell the class everything specified in files and command line
    err = dipole.addParser(parser)

    if err != None:
        print 'error parsing info'
        return

    # read in data 
    dipole.readDipoleProfile()

    dipole.write()

if __name__ == '__main__':
    sys.exit(main())

# "The greatest happiness is to know the source of unhappiness." -Fyodor Dostoevsky
