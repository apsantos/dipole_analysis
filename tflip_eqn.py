#!/usr/bin/env python
"""tflip_eqn

    Find the analytical time to flip a dipole

 -Andrew P. Santos					                            
 
"""
import sys, argparse
import numpy as np
import math as ma        
from scipy import stats
import scipy.special
from scipy.optimize import curve_fit

class DipoleFlipTime(object):
    """Callable Fliptime class

    """
    def __init__(self):
        self.sqrtpi = ma.sqrt( ma.pi ) # sqrt of pi

    def addParser(self, parser):
        """
        Get relevant values from an argparse parser
        """
        if parser.parse_args().mass:
            self.mass = parser.parse_args().mass
        else:
            self.mass = -1

        if parser.parse_args().diameter:
            self.radius = parser.parse_args().diameter / 2.0
        else:
            self.radius = -1

        if parser.parse_args().field:
            self.field = parser.parse_args().field
        else:
            self.field = -1

        if parser.parse_args().moment:
            self.moment_i = parser.parse_args().moment
        else:
            self.moment_i = -1

        if parser.parse_args().moment2:
            self.moment_j = parser.parse_args().moment_j
        else:
            self.moment_j = -1

        if parser.parse_args().separation:
            self.separation = parser.parse_args().separation
        else:
            self.separation = -1

        if parser.parse_args().timescale:
            self.timescale = parser.parse_args().timescale
        else:
            self.timescale = -1

        self.solve_mass = parser.parse_args().solve_mass
        self.solve_separation = parser.parse_args().solve_separation

        if self.solve_mass:
            if self.mass > 0:
                print "Cannot set the mass and solve for the mass"
                return -1
            if self.timescale <= 0:
                print "You need the timescale to calculate the mass"
                return -1
        elif self.mass <= 0:
            print "You need to set the mass"
            return -1

        if self.solve_separation:
            if self.separation > 0:
                print "Cannot set the separation and solve for the separation"
                return -1
        elif self.separation <= 0:
            if self.moment_j > 0:
                print "You need to set the separation"
                return -1

        if self.field <= 0:
            if self.moment_j <= 0:
                print "You need to set at least the field or moment of dipole 2"
                return -1
        else:
            if self.moment_j > 0:
                print "WARNING: you have set both the field and the moment of dipole 2, assuming you want the Field-dipole time"
            
    def MomentInertia(self, m, r):
        return (2.0 / 5.0) * m * r**2.0 

    def quadratic(self, a, b, c):
        return (-b + ma.sqrt(b**2 - 4*a*c)) / (2 * a)
            
    def Mass(self):
        """
        calculate the mass needed to match a timescale of a dipole in an applied field
        """
        self.mass = minimize(self.DipoleField, 1, method='BFGS')
        print ("Mass: %10.8f [g/mol]" % self.mass)

    def Separation(self):
        """
        calculate the separation needed to match the timescale of the dipole-dipole and dipole-field time scales
        """
        I = self.MomentInertia(self.mass, self.radius)
        a = (4.0 * self.self.moment_i**2.0) - (8.0 * ma.pi * self.moment_j * self.self.moment_i**2.0 / I)
        b = 4 * ma.pi * self.moment_i * self.moment_j * self.field**2.0 / I
        c = -self.field**2.0
        self.separation = self.quadratic(a, b, c)
        print ("Separation: %10.8f [Angstrom]" % self.separation)

    def DipoleField(self, mass=None):
        """
        calculate the time for a dipole to go from anti-parallel to perpendicular to an applied field
        """
        if not mass:
            mass = self.mass
            
        I = self.MomentInertia(self.mass, self.radius)
        alpha = self.moment_i * self.field / I
        a = alpha / 2.0
        b = -alpha * (ma.pi)
        c = -ma.pi / 2.0
        self.timescale = self.quadratic(a, b, c)
        print ("t_flip: %10.8f [fs]" % self.timescale)

    def DipoleDipole(self):
        """
        calculate the time for a dipole to go from anti-parallel to perpendicular to another dipole
        """
        I = self.MomentInertia(self.mass, self.radius)
        a = 2.0 * self.moment_i * self.moment_j / (2.0 * I * self.separation**3.0 )
        b = -1.0
        c = -ma.pi / 2.0
        self.timescale = self.quadratic(a, b, c)
        print ("t_flip: %10.8f [fs]" % self.timescale)

def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Fit a density profile and extract valuble information')
    parser.add_argument("--field", type=float, 
                   help='value of applied field (likely kcal/mol-Ampere-Angstrom^2)')
    parser.add_argument("--mass", type=float, 
                   help='value of the mass (likely g/mol)')
    parser.add_argument("--moment", type=float, required=True,
                   help='value of moment (likely Ampere - Angstrom)')
    parser.add_argument("--moment2", type=float, 
                   help='value of moment of the other dipole (Ampere - Angstrom)')
    parser.add_argument("--separation", type=float, 
                   help='value of separation (likely Angstrom)')
    parser.add_argument("--diameter", type=float, required=True,
                   help='value of the particle diameter (likely Angstrom)')
    parser.add_argument("--timescale", type=float,
                   help='value of the timescale which the mass should match (likely femptoseconds)')
    parser.add_argument("--solve_mass", action='store_true', default=False,
                   help='Solve for the mass that would satisify the dipole-Field timescale')
    parser.add_argument("--solve_separation", action='store_true', default=False,
                   help='solve for the separation value which would lead to t_field-dipole = t_dipole-dipole')

    # Initialize the Activity class
    mutflip = DipoleFlipTime()

    # Tell the class everything specified in files and command line
    err = mutflip.addParser(parser)

    if err != None:
        print 'error parsing info'
        return

    if mutflip.solve_mass:
        calc_error = mutflip.Mass()
    elif mutflip.solve_separation:
        calc_error = mutflip.Separation()
    elif mutflip.field > 0:
        calc_error = mutflip.DipoleField()
    else:
        calc_error = mutflip.DipoleDipole()

if __name__ == '__main__':
    sys.exit(main())

# "The greatest happiness is to know the source of unhappiness." -Fyodor Dostoevsky
