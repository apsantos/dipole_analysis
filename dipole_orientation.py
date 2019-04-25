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
    # def __init__(self):

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

        if parser.parse_args().z_collumn >= 0:
            """distance collumn"""
            self.ic_z = parser.parse_args().z_collumn
            self.z_conversion = parser.parse_args().z_conversion

        if parser.parse_args().dipole_magnitude_collumn >= 0:
            self.ic_dmag = parser.parse_args().dipole_magnitude_collumn

        if parser.parse_args().dipole_field_direction_collumn >= 0:
            self.ic_dfield = parser.parse_args().dipole_field_direction_collumn

        if parser.parse_args().L_z:
            self.L_z = parser.parse_args().L_z

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

        max_ic = max( self.ic_z, self.ic_dfield, self.ic_dmag )
        i_line = 0
        for line in ifile:
            i_line += 1
            if (i_line < self.start_line):
                continue

            data = line.strip().split()
            if len(data) <= max_ic:
                raise IOError('Not enough collumns in %s based on input information.' % filename)

            self.z.append( float( data[self.ic_z] ) * self.z_conversion)
            if  float( data[self.ic_dmag] ) == 0:
                icosTheta = 0.0
            else:
                icosTheta = float( data[self.ic_dfield] ) / float( data[self.ic_dmag] ) # * self.L_z

            self.cosTheta.append( icosTheta )
            self.P2.append( 0.5 * ( (3.0 * icosTheta**2.0) - 1.0) )

        ifile.close()

    def write(self, filename=None):
        """
        Write the fit parameters and the fit values
        """
        if filename == None:
            filename = self.ofilename

        ofile = open(filename, 'w')


        ofile.write('# z cosTheta(z) P_2(z)\n')

        for i in range(len(self.cosTheta)):
            ofile.write("%10.8f %10.8f %10.8f\n" % ( self.z[i], self.cosTheta[i], self.P2[i] ) )

        ofile.close()

def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Fit a density profile and extract valuble information')
    parser.add_argument('-f', "--filename", type=str, required=True,
                   help='file with density data.')
    parser.add_argument("-o", "--outfilename", type=str, default='density_profile_fit.dat',
                   help='output file name, assumed to be density_profile_fit.dat')
    parser.add_argument('-l', "--start_line", type=int,
                   help='Set the starting line number.')
    parser.add_argument("--dipole_magnitude_collumn", type=int, 
                   help='Collumn with the total dipole magnitude')
    parser.add_argument("--dipole_field_direction_collumn", type=int, 
                   help='Collumn with the dipole magnitude in the applied field direction')
    parser.add_argument("--z_collumn", type=int, 
                   help='Collumn with the bin locations')
    parser.add_argument("--z_conversion", type=float, default=1.0,
                   help='convert the distances by this conversion')
    parser.add_argument("--L_z", type=float,
                   help='box length in the profile direction')

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
