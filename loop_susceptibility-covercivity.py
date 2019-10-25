#!/usr/bin/env python
"""loop_susceptibility-coercivity

 Calculate the susceptibility and coercivity
 							                                    
 -Andrew P. Santos					                            
 
"""
import sys, argparse
import numpy as np
import math as ma        
from scipy import stats
import scipy.special
from scipy.optimize import curve_fit

class MH(object):
    """Callable Hysteresis loop class

       Reads M(H) and calculates
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

        self.fittype = parser.parse_args().fittype
        if parser.parse_args().fitrange:
            self.fitrange = True
            self.fitmin = parser.parse_args().fitrange[0]
            self.fitmax = parser.parse_args().fitrange[1]
        elif self.fittype == 'linear':
            print 'WARNING: you should probably set --fitrange if using a linear fit'
        else:
            self.fitrange = False

        if parser.parse_args().time_collumn >= 0:
            self.ic_t = parser.parse_args().time_collumn
            self.t_conversion = parser.parse_args().time_conversion
        else:
            self.ic_t = -1

        if parser.parse_args().B_collumn >= 0:
            self.ic_b = parser.parse_args().B_collumn
            self.BtoH_conversion = parser.parse_args().BtoH_conversion
        else:
            if parser.parse_args().H_collumn < 0:
                print "H or B collumn needs to be defined"
                return -1
            self.ic_b = -1

        if parser.parse_args().B_max >= 0:
            self.bmag_max = parser.parse_args().B_max
        else:
            self.bmag_max = 0

        if parser.parse_args().H_collumn >= 0:
            self.ic_h = parser.parse_args().H_collumn
        else:
            self.ic_h = -1

        if parser.parse_args().M_collumn >= 0:
            self.ic_m = parser.parse_args().M_collumn
        else:
            self.ic_m = -1

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
        self.h = []
        self.m = []

        max_ic = max( self.ic_t, self.ic_b, self.ic_h, self.ic_m )
        i_line = 0
        for line in ifile:
            i_line += 1
            if (i_line < self.start_line):
                continue

            data = line.strip().split()
            if len(data) <= max_ic:
                raise IOError('Not enough collumns in %s based on input information.' % filename)

            if self.ic_t >= 0:
                self.t.append( float( data[self.ic_t] ) * self.t_conversion)
            if self.ic_b >= 0:
                self.b.append( float( data[self.ic_b] ))
            if self.ic_h >= 0:
                self.h.append( float( data[self.ic_h] ))
            if self.ic_m >= 0:
                self.m.append( float( data[self.ic_m] ))

        ifile.close()

    def fit_linear(self, h, m):
        """
        Fit a linear equation to M(H)
        M = a*H+c
        susceptbility = a
        c is like the coercivity
        """
        susceptibility, coercivity, r, p, se = scipy.stats.linregress(h, m)
        return coercivity/susceptibility, susceptibility
    
    def langevin_derivative(self, h, a, b, c):
        x = b*(h+c)
        return -(a*b) * (np.sinh(x)**-2.0 - x**-2.0)
        #return a * (np.linalg.matrix_power(np.sinh(x),-2.0) - np.linalg.matrix_power(x, -2.0))

    def langevin(self, h, a, b, c):
        x = b*(h+c)
        return a * (1/np.tanh(x) - 1/(x))
        #return a * (np.linalg.matrix_power(np.tanh(x), -1.0) - np.linalg.matrix_power(x, -1.0))

    def fit_langevin(self, h, m):
        """
        Fit the langevin equation to M(B)
        M = A[coth(x)-1/x]
        dMdB = A/(B(H +- C))^2 - ABcsch^2(B(H +- C))

        x = B(H +- C)
        C is like the coercivity
        A should be Msat
        B should be mu_0*mu*beta 
        """
        # get guess for magnitude and frequency
        msatguess = max(m)
        bguess = 4.0 * (ma.pi * 1439.5315)**2.0 * 151 / ( 1.9872e-3 *298 )
        Bopt, Bcov = curve_fit(self.langevin, h, m, p0 = [ msatguess, bguess, 0] )
        A = Bopt[0]
        B = Bopt[1]
        C = Bopt[2]
        A_std = np.sqrt(np.diag(Bcov))[0]
        B_std = np.sqrt(np.diag(Bcov))[1]
        C_std = np.sqrt(np.diag(Bcov))[2]
        AB_std = Bcov[0,1]
        m_fit = self.langevin(h, A, B, C )
        susceptibility = self.langevin_derivative( 0, A, B, C )
        return A, B, C, susceptibility
    
    def fit_segments(self):
        if self.ic_h < 0:
           self.h = np.array(self.b) * self.BtoH_conversion 
        else:
            self.h = np.array(self.h)

        self.m = np.array(self.m)

        Hmax = np.max(self.h)
        Mmax = np.max(self.m)
        hs = [[]]
        ms = [[]]
        seg = 0
        pdiffh = self.h[1] - self.h[0]

        for i in range(1,len(self.h)):
            # if you know how high B goes, turn around once it hits the edge of a sweep
            # otherwise when the sweep starts to go the opposite direction start recording another segment
            if self.bmag_max:
                if abs(self.b[i]) > (self.bmag_max):
                    continue
            diffh = self.h[i] - self.h[i-1]
            if np.sign(diffh) != np.sign(pdiffh):
                #if (self.bmag_max and (self.b[i]) < self.bmag_max)) or (not self.bmag_max) :
                    seg += 1
                    hs.append( [] )
                    ms.append( [] )
                    pdiffh = diffh

            if self.fitrange:
                if (self.b[i] < self.fitmin) or (self.b[i] > self.fitmax):
                    continue

            hs[seg].append(self.h[i])
            ms[seg].append(self.m[i])

        coercivity = []
        susceptibility = []
        for iseg in range(1,seg):
            h_seg = np.array(hs[iseg])
            m_seg = np.array(ms[iseg])
            if self.fittype == 'langevin':
                ia, ib, ic, iss = self.fit_langevin( h_seg, m_seg)
            elif self.fittype == 'linear':
                ic, iss = self.fit_linear( h_seg, m_seg)

            coercivity.append(ic)
            susceptibility.append(iss)

        self.coercivity_mean = np.mean(np.abs(coercivity))
        self.susceptibility_mean = np.mean(susceptibility)
        self.coercivity_std = np.std(np.abs(coercivity))
        self.susceptibility_std = np.std(susceptibility)

    def write(self, filename=None):
        """
        Write the fit parameters and the fit values
        """
        if filename == None:
            filename = self.ofilename

        ofile = open(filename, 'w')

        ofile.write('# Susceptibility: %E d(susc): %E Coercivity: %E d(coer): %E\n' % (self.susceptibility_mean, self.susceptibility_std, self.coercivity_mean, self.coercivity_std) )
        ofile.write('# H[] M[] Mfit[]\n')

        #for i in range(len(self.h)):
        #    ofile.write(" %12.10f %12.10f %12.10f\n" % ( self.h[i], self.m[i], self.m_fit[i] ) )

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
    parser.add_argument("--time_collumn", type=int, 
                   help='Collumn with the simulation time step')
    parser.add_argument("--time_conversion", type=float, default=1.0,
                   help='convert the timegt by this value')
    parser.add_argument("--B_collumn", type=int, 
                   help='Collumn with the applied field')
    parser.add_argument("--B_max", type=float,
                   help='maximum b magnitude')
    parser.add_argument("--BtoH_conversion", type=float, default=18089.46015,
                   help='convert the timegt by this value')
    parser.add_argument("--H_collumn", type=int, 
                   help='Collumn with the applied field')
    parser.add_argument("--M_collumn", type=int,  required=True,
                   help='Collumn with the magnetization')
    parser.add_argument("--fittype", type=str, default='langevin', choices=['langevin','linear'],
                   help='either fit the entire curve to the Langevin equation or a specified linear region near 0')
    parser.add_argument("--fitrange", nargs=2, type=float,
                   help='set the min and max B for the fit (argparse will not accept negative scientific notation values)')

    # Initialize the MH
    loop = MH()

    # Tell the class everything specified in files and command line
    err = loop.addParser(parser)

    if err != None:
        print 'error parsing info'
        return

    # read in data 
    loop.readBMdata()

    # fit to the langevin curve
    loop.fit_segments()

    loop.write()

if __name__ == '__main__':
    sys.exit(main())

# "The greatest happiness is to know the source of unhappiness." -Fyodor Dostoevsky
