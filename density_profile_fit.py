#!/usr/bin/env python
"""density_profile_fit

 fit a VLE density profile 
 							                                    
 -Andrew P. Santos					                            
 
"""
import sys, argparse
import numpy as np
import math as ma        
from scipy import stats
import scipy.special
from scipy.optimize import curve_fit

class DensityProfileFit(object):
    """Callable Density Profile class

       Reads  density profile and fits
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
            self.ofilename = 'density_profile_fit.dat'

        if parser.parse_args().z_collumn >= 0:
            """distance collumn"""
            self.ic_z = parser.parse_args().z_collumn
            self.z_conversion = parser.parse_args().z_conversion

        if parser.parse_args().ndensity_collumn >= 0:
            """denisty collumn in N/V"""
            self.ic_rho = parser.parse_args().ndensity_collumn

        elif parser.parse_args().N_collumn >= 0:
            self.ic_rho = parser.parse_args().N_collumn

        else:
            print 'you need either the collumn for N or rho in the file'
            return 1

        self.density_conversion = parser.parse_args().density_conversion

        if parser.parse_args().L_z:
            self.L_z = parser.parse_args().L_z
            self.half_L_z = parser.parse_args().L_z / 2.0

        if parser.parse_args().fit == "IGS":
            self.fit_style = parser.parse_args().fit
            
    def readDensity(self, filename=None):
        """
        Read density profile
        """
        if filename == None:
            filename = self.ifilename

        try:
            ifile = open(filename, 'r')
        except IOError:
            raise IOError('cannot find: %s. Either put the file in the working directory or fix the input file.' % filename)

        self.z = []
        self.rho = []

        i_line = 0
        shift = 0
        for line in ifile:
            i_line += 1
            if (i_line < self.start_line):
                continue

            data = line.strip().split()
            if len(data) <= self.ic_rho:
                raise IOError('Not enough collumns in %s based on input information.' % filename)

            if shift == 0:
                if float( data[self.ic_z] ) < 0:
                    # the positions need to be from 0 to L_z
                    shift = self.half_L_z
            self.z.append( float( data[self.ic_z] ) * self.z_conversion)
            self.rho.append( float( data[self.ic_rho] ) * self.density_conversion)
        self.z = [x+shift for x in self.z]

        ifile.close()

    def IGS_fit(self, z, rho_v, rho_l, w, z_0):
        sqrtpi_w = self.sqrtpi / w
        return (rho_v + (0.5 * (rho_l - rho_v) * (
                           scipy.special.erf( sqrtpi_w * (z - self.half_L_z + z_0) ) 
                         - scipy.special.erf( sqrtpi_w * (z - self.half_L_z - z_0) ) )))

    def fit_profile(self, filename=None):
        """
        Fit the density profile
        """

        if (self.fit_style == 'IGS'):
            popt, pcov = curve_fit(self.IGS_fit, np.array(self.z), np.array(self.rho) )
        self.rho_vap = popt[0]
        self.rho_liq = popt[1]
        self.w = popt[2]
        self.z_0 = popt[3]
        self.delta = self.w**2.0 / 2.0 / ma.pi

        if (self.fit_style == 'IGS'):
            self.rho_fit = self.IGS_fit( np.array(self.z), self.rho_vap, self.rho_liq, self.w, self.z_0 )


    def write(self, filename=None):
        """
        Write the fit parameters and the fit values
        """
        if filename == None:
            filename = self.ofilename

        ofile = open(filename, 'w')


        ofile.write('# z rho(z) rho_fit(z) fitted with:')

        if self.fit_style == 'IGS':
            ofile.write('# rho(z) = rho_vap + (rho_liq-rho_vap)/2 * '
                         '{erf[sqrt(\pi)/w*(z-L_z/2+z_0)] - erf[sqrt(\pi)/w*(z-L_z/2-z_0)]}\n')
            ofile.write('# rho_vap = %f\n' % self.rho_vap)
            ofile.write('# rho_liq = %f\n' % self.rho_liq)
            ofile.write('# w = %f\n' % self.w)
            ofile.write('# z_0 = %f\n' % self.z_0)
            ofile.write('# D = %f = w^2/2pi\n' % self.delta)

        for i in range(len(self.rho)):
            ofile.write("%10.8f %10.8f %10.8f\n" % ( self.z[i], self.rho[i], self.rho_fit[i]) )
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
    parser.add_argument("--ndensity_collumn", type=int, 
                   help='Collumn with the number density')
    parser.add_argument("--N_collumn", type=int, 
                   help='Collumn with the number of particles in a bin')
    parser.add_argument("--density_conversion", type=float, default=1.0,
                   help='convert the density or n collumn by this value')
    parser.add_argument("--z_collumn", type=int, 
                   help='Collumn with the bin locations')
    parser.add_argument("--z_conversion", type=float, default=1.0,
                   help='convert the distances by this conversion')
    parser.add_argument("--L_z", type=float,
                   help='box length in the profile direction')
    parser.add_argument("--fit", type=str, choices=["IGS"],
                   help='fit to the Ismail-Grest-Stevens functional: '
                         'rho(z) = rho_vap + (rho_liq+rho_vap)/2 * '
                         '{erf[sqrt(\pi)/w*(z-L_z/2+z_0)] - erf[sqrt(\pi)/w*(z-L_z/2-z_0)]}')

    # Initialize the Activity class
    rhofit = DensityProfileFit()

    # Tell the class everything specified in files and command line
    err = rhofit.addParser(parser)

    if err != None:
        print 'error parsing info'
        return

    # read in data 
    rhofit.readDensity()
    calc_error = rhofit.fit_profile()

    if not calc_error:
        rhofit.write()

if __name__ == '__main__':
    sys.exit(main())

# "The greatest happiness is to know the source of unhappiness." -Fyodor Dostoevsky
