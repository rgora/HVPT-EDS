#!/usr/bin/env python
"""
This is a simple EDS helper script. By default it prepares a set of
input files for EDS calculations of each xyz structure provided as an
input, in which case an input template is required. If none is found a standard
one is prepared. Invoking the script with -c option allows to retrieve the
properties, provided that the calculations are completed, in which case paths
to data dir(s) are expected on input.

Usage: xyz2eds.py [options] xyz file(s) or data dir(s)

Options:
  -h, --help       show this help


"""

#     Copyright (C) 2011, Robert W. Gora (robert.gora@pwr.wroc.pl)
#  
#     This program is free software; you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation; either version 2 of the License, or (at your
#     option) any later version.
#  
#     This program is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#     Public License for more details.
#  
#     You should have received a copy of the GNU General Public License along
#     with this program; if not, write to the Free Software Foundation, Inc.,
#     675 Mass Ave, Cambridge, MA 02139, USA.

__author__  = "Robert Gora (robert.gora@pwr.wroc.pl)"
__version__ = list(filter(str.isdigit, "$Revision$"))

# Import necessary modules
import os, sys, getopt, re

from string import Template
from numpy import *

# Regular expressions
reflags = re.DOTALL

#----------------------------------------------------------------------------
# Usage
#----------------------------------------------------------------------------
def Usage():
    """Print usage information and exit."""
    print(__doc__)
    print("Machine epsilon is: ",finfo(float64).eps,"for float64 type\n")

    sys.exit()

#----------------------------------------------------------------------------
# Main
#----------------------------------------------------------------------------
def Main(argv):
    '''Parse commandline and loop throught the logs'''

    # Parse commandline
    try:
        opts, args = getopt.getopt(argv, "h",
                                        ["help"])
    except getopt.GetoptError as error:
        print(error)
        Usage()
    if not argv:
        Usage()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            Usage()

    # Parse each data file (with xyz coords) or dir (with the results)
    data_files = args

    for f in data_files:
        GAMESS_INPUTS(f)

#----------------------------------------------------------------------------
# Common input template routines
#----------------------------------------------------------------------------
class INPUT_TEMPLATE(Template):
    delimiter = '@'

class INPUTS:
    """Common input routines"""
    def __init__(self, data):
        self.data = data
        self.ReadTemplate()

    def ReadTemplate(self):
        """Read or punch standard template"""
        try:
            self.tmpl = open(self.pkg+'.tmpl','r').read()
            self.tmpl = INPUT_TEMPLATE(self.tmpl)
            self.WriteInputs()
        except IOError:
            print("There's no " + self.pkg + " template. I'm punching one - please check")
            open(self.pkg+'.tmpl','w').write(self.tmpl)
            sys.exit()
        except AttributeError:
            pass

    def WriteInputs(self):
        pass


#----------------------------------------------------------------------------
# Gamess (US) routines
#----------------------------------------------------------------------------

class GAMESS_INPUTS(INPUTS):
    """Gamess US input routines"""

    def __init__(self, data):
        # template name
        self.pkg = "gamess"
        # template content
        self.tmpl="""\
 $system mwords=10 memddi=100 parall=.t. $end
 $contrl scftyp=rhf runtyp=eds icharg=0 mult=1 units=angs
         maxit=100 nprint=7 exetyp=run mplevl=2 ispher=1
         icut=20 itol=30 aimpac=.f. cctyp=none $end
 $eds    mch(1)=0,0 mmul(1)=1,1 mnr(1)=3 $end
 $mp2    mp2prp=.t. code=ddi cutoff=1d-20 $end
 $ccinp  iconv=14 maxcc=100 $end
 $trans  cuttrf=1d-15 $end
 $scf    dirscf=.t. fdiff=.f. diis=.t. soscf=.f.
         conv=1d-11 swdiis=0.0001 npreo(1)=0,-1,1,9999 $end
 $basis  gbasis=ccd $end
@data
"""
        INPUTS.__init__(self, data)

    def WriteInputs(self):

        # initialize periodic table
        p=Periodic(0)
        
        # read xyz file
        try:
            xyz=open(self.data,'r').readlines()
            xyz=xyz[2:int(xyz[0])+2]
        except ValueError:
            print("Problem with *.xyz file?")
            sys.exit(1)
        
        for i in range(len(xyz)):
            xyz[i]=xyz[i].split()
            xyz[i].insert(1, str( Atomn(xyz[i][0], p) ))
            xyz[i]='%-5s %5s %15s %15s %15s\n' % tuple([xyz[i][j] for j in range(5)])
        
        xyz.insert(0, ' $data\nEDS\nc1 0\n')
        xyz.append(' $end')
        xyz=''.join(xyz)
        
        # write input files
        filename = self.data.replace('.xyz','')
        filename = filename.replace(' ','_')
        finput = self.tmpl.substitute(data=xyz)
        # write inputs
        open(filename+'.inp','w').write(finput)


def Periodic(mendeleiev):
    '''Returns the mendeleiev table as a python list of tuples. Each cell
    contains either None or a tuple (symbol, atomic number), or a list of pairs
    for the cells * and **. Requires: "import re". Source: Gribouillis at
    www.daniweb.com - 2008 '''

    # L is a consecutive list of tuples ('Symbol', atomic number)
    L = [ (e,i+1) for (i,e) in enumerate( re.compile ("[A-Z][a-z]*").findall('''
    HHeLiBeBCNOFNeNaMgAlSiPSClArKCaScTiVCrMnFeCoNiCuZnGaGeAsSeBrKr
    RbSrYZrNbMoTcRuRhPdAgCdInSnSbTeIXeCsBaLaCePrNdPmSmEuGdTbDyHoEr
    TmYbLuHfTaWReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaUNpPuAmCmBkCfEsFm
    MdNoLrRfDbSgBhHsMtDsRgUubUutUuqUupUuhUusUuo'''))]

    # The following fills the void with nones and returns the list of lists
    mendeleiev = 0

    if mendeleiev:
        for i,j in ( (88,103), (56,71) ):
            L[i] = L[i:j]
            L[i+1:] = L[j:]
        for i,j in ( (12,10 ), (4,10), (1,16) ):
            L[i:i]=[None]*j 

        return [ L[18*i:18*(i+1)] for i in range (7) ]

    # Return a plain list of tuples
    else:
        return L

def Atomn(s,ptable):
    '''Returns the atomic number based on atomic symbol string
    ptable is a list of consecutive (symbol, atomic number) tuples.'''

    for n,a in enumerate(ptable):
        if a[0].lower().find(s.strip().lower()) !=-1 :
            return float(n+1)


#----------------------------------------------------------------------------
# Main routine
#----------------------------------------------------------------------------
if __name__ == "__main__":
    Main(sys.argv[1:])

