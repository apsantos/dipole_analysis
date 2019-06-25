#!/usr/bin/env python
"""dipole_orientation

 Calculate the AC susceptibilities(chi_i, chi_ii)
 							                                    
 -Andrew P. Santos					                            
 
"""
import sys, argparse
import numpy as np
import math as ma        
from scipy import stats
import scipy.special
from scipy.optimize import curve_fit

class ACSusceptibility(object):
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

        if parser.parse_args().time_collumn >= 0:
            self.ic_t = parser.parse_args().time_collumn
            self.t_conversion = parser.parse_args().time_conversion

        if parser.parse_args().B_collumn >= 0:
            self.ic_b = parser.parse_args().B_collumn

        if parser.parse_args().M_collumn >= 0:
            self.ic_m = parser.parse_args().M_collumn

    def readBMdata(self, filename=None):
        """
        Read the dipole profile
        """
        if filename == None:
            filename = self.ifilename

        try:
            ifile = open(filename, 'r')
        except IOError:
            raise IOError('cannot find: %s. Either put the file in the working directory or fix the input file.' % filename)

        self.t = []
        self.b = []
        self.m = []

        max_ic = max( self.ic_t, self.ic_b, self.ic_m )
        i_line = 0
        for line in ifile:
            i_line += 1
            if (i_line < self.start_line):
                continue

            data = line.strip().split()
            if len(data) <= max_ic:
                raise IOError('Not enough collumns in %s based on input information.' % filename)

            self.t.append( float( data[self.ic_t] ) * self.t_conversion)
            self.b.append( float( data[self.ic_b] ))
            self.m.append( float( data[self.ic_m] ))

        ifile.close()

    def sine_fit_fixfrequency(self, t, mag, shift):
        return mag * np.sin( self.omega * (t - shift) )

    def sine_fit(self, t, mag, frequency, shift):
        return mag * np.sin( frequency * (t - shift) )

    def fit(self):
        """
        Fit sine functions to B and M curves
        """
        self.t = np.array(self.t)
        self.b = np.array(self.b)
        self.m = np.array(self.m)

        # get guess for magnitude and frequency
        magguess = np.max(self.b)
        mean = np.mean(self.b)
        fguess = 1
        i_1 = 0
        over_hump = False
        for i in range( 15, len(self.t)):
            # have we past the region near 0?
            if i_1 != 0:
                if (self.b[i] - mean)/magguess >= 0.01:
                    over_hump = True

            # is B near zero?
            if (self.b[i] - mean)/magguess < 0.01:
                # If this is the 1st time, keep track
                if i_1 == 0:
                    print 'hi', self.b[i]
                    i_1 = i
                # If it is the second time, we have the separation from 0
                elif over_hump:
                    print 'ho', (self.t[i] - self.t[i_1]), self.b[i]
                    # the frequency
                    fguess = 3.0 / (self.t[i] - self.t[i_1])
                    break
        print fguess
        shiftguess = 0
        Bopt, Bcov = curve_fit(self.sine_fit, self.t, self.b, p0 = [ magguess, fguess, shiftguess] )
        self.b0 = Bopt[0]
        self.omega = Bopt[1]
        self.tau = Bopt[2]
        self.b_fit = self.sine_fit( self.t, self.b0, self.omega, self.tau )

        magguess = np.max(self.b)
        # use the same shiftguess and fguess as the m
        Mopt, Mcov = curve_fit(self.sine_fit_fixfrequency, self.t, self.m, p0 = [ magguess, shiftguess] )
        self.m0 = Mopt[0]
        #self.omega_m = Mopt[1]
        self.tau_ref = Mopt[1]
        self.m0_std = np.sqrt(np.diag(Mcov))[0]
        self.tau_ref_std = np.sqrt(np.diag(Mcov))[1]
        self.m0_tau_ref_std = Mcov[0,1]
        self.m_fit = self.sine_fit_fixfrequency( self.t, self.m0, self.tau_ref )
        #self.m_fit = self.sine_fit_fixfrequency( self.t, self.m0, self.omega_m, self.tau_ref )

        self.chi_i = (self.m0/self.b0) * ma.cos(self.omega * self.tau_ref )
        self.chi_i_std = abs(self.chi_i) * ma.sqrt( 
                            (self.m0_std/ self.m0)**2.0
                          + (self.omega * ma.tan(self.omega * self.tau_ref) * self.tau_ref_std)**2.0
                          + (2.0 * self.m0_tau_ref_std / self.chi_i) )
                                                    
        self.chi_ii = (self.m0/self.b0) * ma.sin(self.omega * self.tau_ref )
        self.chi_ii_std = abs(self.chi_i) * ma.sqrt( 
                            (self.m0_std/ self.m0)**2.0
                          + (self.omega / ma.tan(self.omega * self.tau_ref) * self.tau_ref_std)**2.0
                          + (2.0 * self.m0_tau_ref_std / self.chi_i) )
        self.chi = ma.sqrt( self.chi_i**2.0 + self.chi_ii**2.0)
        self.chi_std = ma.sqrt( ((self.chi_i/self.chi)**2.0 * self.chi_i_std**2.0) 
                              + ((self.chi_ii/self.chi)**2.0 * self.chi_ii_std**2.0) )

    def write(self, filename=None):
        """
        Write the fit parameters and the fit values
        """
        if filename == None:
            filename = self.ofilename

        ofile = open(filename, 'w')


        ofile.write('# B0: %e omega: %e T: %e\n' % (self.b0, self.omega, self.tau) )
        ofile.write('# M0: %e Tref: %e\n' % (self.m0, self.tau_ref) )
        ofile.write('# d(M0): %e d(Tref): %e\n' % (self.m0_std, self.tau_ref_std) )
        ofile.write('# Chi: %e Chi_i: %e Chi_ii: %e\n' % (self.chi, self.chi_i, self.chi_ii) )
        ofile.write('# d(Chi): %e d(Chi_i): %e d(Chi_ii): %e\n' % (self.chi_std, self.chi_i_std, self.chi_ii_std) )
        ofile.write('# t b m m(fit)\n')

        for i in range(len(self.t)):
            ofile.write(" %12.10f %12.10f %12.10f %12.10f %12.10f\n" % ( self.t[i], self.b[i], self.b_fit[i], self.m[i], self.m_fit[i] ) )

        ofile.close()

def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Fit a density profile and extract valuble information')
    parser.add_argument('-f', "--filename", type=str, required=True,
                   help='file with applied field(B) and magnetization(M) data.')
    parser.add_argument("-o", "--outfilename", type=str, default='ac_susceptibility_fit.dat',
                   help='output file name, assumed to be density_profile_fit.dat')
    parser.add_argument('-l', "--start_line", type=int,
                   help='Set the starting line number.')
    parser.add_argument("--time_collumn", type=int,  required=True,
                   help='Collumn with the simulation time step')
    parser.add_argument("--time_conversion", type=float, default=1.0,
                   help='convert the timegt by this value')
    parser.add_argument("--B_collumn", type=int,  required=True,
                   help='Collumn with the applied field')
    parser.add_argument("--M_collumn", type=int,  required=True,
                   help='Collumn with the magnetization')

    # Initialize the Activity class
    susc = ACSusceptibility()

    # Tell the class everything specified in files and command line
    err = susc.addParser(parser)

    if err != None:
        print 'error parsing info'
        return

    # read in data 
    susc.readBMdata()

    # fit to sine waves
    susc.fit()

    susc.write()

if __name__ == '__main__':
    sys.exit(main())

# "The greatest happiness is to know the source of unhappiness." -Fyodor Dostoevsky
