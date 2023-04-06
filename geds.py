#!/usr/bin/env python

"""
Parser for EDS GAMESS module

Usage: ggms.py [options] [output file(s)]

Options:
  -h, --help            show this help

  -o, --out             Output format: txt, csv, tex

  -t                    Grep total energies

  -e, --energy-units=   set energy units; chose from: kcal, kJ, meV,
                        mH (mili Hartree) or au which is the default

  -p, --property-units= set energy units; chose from: esu, si,
                        asi (si multiplied by electric permittivity of free
                        space, mau (mili au) or au which is the default

  -s, --sort=           select alternative mode of sorting: float
                        it looks for '\d+\.\d+' substrings and compares
                        the values of first occurences; if such strings are
                        not present in the filename the normal sorting is resumed

  -r, --relative=       the absolute values are converted to relative
                        if relative=first the first value is used as a reference
                        if relative=last its the opposite

  -d                    show debugging information while parsing
"""

__author__  = "Robert Gora (robert.gora@pwr.wroc.pl)"
__version__ = "$Revision: 0.2 $"

from numpy import *

import sys
import getopt
import re

# Regular expressions
reflags = re.DOTALL

#----------------------------------------------------------------------------
# Usage
#----------------------------------------------------------------------------
def Usage():
    """Print usage information."""
    print(__doc__)

#----------------------------------------------------------------------------
# Main
#----------------------------------------------------------------------------
def Main(argv):
    """Parse commandline and loop throught the logs"""

    global Monomers, Systems, SortMode, Relative, OutFormat
    global _TotEn_, TotOutFile


    Energies      = {}
    TotEnergies   = {}
    Properties    = {}
    EnergyTables  = {}
    ClusterTables = {}
    MbodyTables   = {}
    FieldTables   = {}
    PropTables    = {}
    TitleLen      = [25]
    SortMode      = False
    Relative      = ''
    OutFormat     = SetOutFormat('txt')

    # Dictionary of sorted labels
    OrdLabel = SetLabels()
    OldLabel = SetLabels()

    # Set units
    SetEnUnits('au')
    SetPrUnits('au')

    _TotEn_ = 0

    # Parse commandline
    try:
        opts, args = getopt.getopt(argv, "ho:e:p:s:r:dt", 
                                        ["help",
                                         "out=",
                                         "energy-units=",
                                         "property-units=",
                                         "sort=",
                                         "relative="])
    except getopt.GetoptError :
        Usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            Usage()
            sys.exit()
        elif opt == '-d':
            global _debug
            _debug = 1
        elif opt == '-t':
            _TotEn_ = 1
            TotOutFile=open('toten.txt','w')
        elif opt in ("-o", "--out"):
            OutFormat = SetOutFormat(arg)
        elif opt in ("-e", "--energy-units"):
            EnUnits = SetEnUnits(arg)
        elif opt in ("-p", "--property-units"):
            PrUnits = SetPrUnits(arg)
        elif opt in ("-s", "--sort"):
            SortMode=arg
        elif opt in ("-r", "--relative"):
            Relative=arg

    if not args:
        Usage()
        sys.exit()

    # Parse each log file
    for LogFile in args:
        OrdLabel = SetLabels()
        ParseFile(LogFile,OrdLabel,TitleLen,Energies,Properties,TotEnergies)
        Labels = SaveLabels(OrdLabel,OldLabel)

    # Format results
    FormatSubEnergies(EnergyTables,ClusterTables,Energies,Labels)
    
    if ManyBody:
        FormatMnbEnergies(MbodyTables,Energies,Labels)

    if FiniteField:
        FormatFieldEnergies(FieldTables,Energies,Labels)
        FormatProperties(PropTables,Properties,Labels)

    # Write results
    WriteEnergies(max(TitleLen),EnergyTables,ClusterTables,MbodyTables,FieldTables)
    if FiniteField:
        WriteProperties(max(TitleLen),PropTables,Properties)

#----------------------------------------------------------------------------
# Parse
#----------------------------------------------------------------------------
def ParseFile(LogFile,OrdLabel,TitleLen,Energies,Properties,TotEnergies):
    """Parse current log file"""

    # open the current log file ...
    File=open(LogFile,'r')

    # ... parse ...
    Title = ReadPreamble(File,LogFile)
    TitleLen.append(len(Title.split()[1].replace('.log',''))+1)

    # ... read field free energies of subsystems ...
    Energies[Title] = {}
    Field=(0,0,0)
    Energies[Title][Field]={}

    for i in range(Systems-Monomers): 
        ReadSubEnergies(File,'',OrdLabel['SubLabel'],Energies[Title][Field])
        if MpLevel == 2 or (CcLevel and CcLevel.upper() != 'NONE'):
            ReadSubEnergies(File,'(CORR)',OrdLabel['SubLabel'],Energies[Title][Field])

    # ... read many-body partitioning ...
    if ManyBody:
        File.seek(0)
        if MpLevel == 2 or (CcLevel and CcLevel.upper() != 'NONE'):
            SeekIter=3
        else:
            SeekIter=2
        ReadMnbEnergies(File,OrdLabel['MnbLabel'],Energies[Title][Field],SeekIter)

    # ... read properties ...
    if FiniteField:

        File.seek(0)

        # ... energies in field ...
        while 1:
            line = File.readline()
            if line == '': break
            if line.find('APPLIED FIELD') !=-1:
                FirstLine = line
                line = line.split()

                # ... set field label ...
                Field = (float(line[-3]), float(line[-2]), float(line[-1]))
                if (Field not in Energies[Title]) and OrdLabel['FieldLabel'].count(Field) == 0 :
                    OrdLabel['FieldLabel'].append(Field)
                    Energies[Title][Field]={}

                # ... read energies of subsystems ...
                #for i in range(Systems-Monomers):
                #    print Energies
                #    ReadSubEnergies(File,'',OrdLabel['SubLabel'],Energies[Title][Field])
                #    if MpLevel == 2 or (CcLevel and CcLevel.upper() != 'NONE'):
                #        ReadSubEnergies(File,'(CORR)',OrdLabel['SubLabel'],Energies[Title][Field])
               
                # ... read many-body partitioning ...
                File.seek(0)
                FindLine(File, FirstLine)                    

                if ManyBody:
                    ReadMnbEnergies(File,OrdLabel['MnbLabel'],Energies[Title][Field],1)
                else:
                    ReadFldEnergies(File,OrdLabel['MnbLabel'],Energies[Title][Field],1)

        # ... FEDS properties ...
        Properties[Title] = {}
        ReadProperties(File,OrdLabel,Properties[Title])

    # ... Total energies
    if _TotEn_:
        TotEnergies[Title] = {}
        ReadTotEnergies(File,TotOutFile,TotEnergies[Title])

    # ... and close
    File.close()

#----------------------------------------------------------------------------
# Set Label
#----------------------------------------------------------------------------
def SetLabels():
    """Reset labels"""

    OrdLabel = { 'SubLabel': [],
                 'MnbLabel': [],
                 'FieldLabel': [],
                 'TotPropLabel': [],
                 'IntPropLabel': [],
                 'SumPropLabel': [],
                 'ExcPropLabel': [] }

    return OrdLabel

def SaveLabels(OrdLabel,OldLabel):
    """Set labels"""

    for LabelType in list(OrdLabel.keys()):
        if len(OrdLabel[LabelType]) > len(OldLabel[LabelType]):
            OldLabel[LabelType] = OrdLabel[LabelType]
        else:
            OrdLabel[LabelType] = OldLabel[LabelType]

    return OrdLabel

#----------------------------------------------------------------------------
# Read Properties
#----------------------------------------------------------------------------
def ReadProperties(File,OrdLabel,Properties):
    """Read selected properties for this system."""

    File.seek(0)

    # Read total properties
    Properties['Total'] = {}

    line = FindLine(File,'ELECTRIC PROPERTIES OF SUBSYSTEMS')
    line = SkipLines(File,1)

    #for i in range(4*Monomers+5*(Systems-Monomers-1)):
    #    line = FindLine(File,'BASED RESULTS')
    #for i in range(5):
    #    ReadTotalProperty(File,OrdLabel['TotPropLabel'],Properties['Total'])

    while 1:
        line = File.readline()
        if line == '': break
        if line.find(30*'=') !=-1: break
        if line.find('BASED RESULTS') !=-1:
            ReadTotalProperty(File,line,OrdLabel['TotPropLabel'],Properties['Total'])

    # Read interaction induced properties
    Properties['Interaction'] = {}

    line = FindLine(File,'  INTERACTION INDUCED PROPERTIES')
    line = SkipLines(File,1)

    while 1:
        line = File.readline()
        if line == '': break
        if line.find(30*'=') !=-1: break
        if line.find('BASED RESULTS') !=-1:
            ReadInteractionProperty(File,line,OrdLabel['IntPropLabel'],Properties['Interaction'])

    # Read sum of interaction induced properties
    Properties['SumInteraction'] = {}

    line = FindLine(File,'SUM OF INTERACTION INDUCED PROPERTIES')
    line = SkipLines(File,1)

    while 1:
        line = File.readline()
        if line == '': break
        if line.find(30*'=') !=-1: break
        if line.find('BASED RESULTS') !=-1:
            ReadInteractionProperty(File,line,OrdLabel['SumPropLabel'],Properties['SumInteraction'])

    # Read excess properties
    Properties['Excess'] = {}

    line = FindLine(File,'EXCESS PROPERTIES')
    line = SkipLines(File,1)

    #for i in range(5*(Systems-Monomers-1)):
    #    line = FindLine(File,'BASED RESULTS')
    #for i in range(5):
    #    ReadTotalProperty(File,OrdLabel['ExcPropLabel'],Properties['Excess'])

    while 1:
        line = File.readline()
        if line == '': break
        if line.find(30*'=') !=-1: break
        if line.find('BASED RESULTS') !=-1:
            ReadTotalProperty(File,line,OrdLabel['ExcPropLabel'],Properties['Excess'])

#----------------------------------------------------------------------------
# Read Interaction Induced Properties
#----------------------------------------------------------------------------
def ReadInteractionProperty(File,line,OrdLabel,Properties):
    """Read selected interaction property."""

    Label = line.split()

    if len(Label) == 3:
        Label = Label[0]
    elif Label[1].count('CC') > 0:
        Label = ','.join([''.join(Label[:3]), Label[-2].split('-')[0]])
    else:
        Label = ','.join([Label[0], Label[-2].split('-')[0]])

    OrdLabel.append(Label)
    PrUnits['LabLen'].append(len(Label))
    if OutFormat == 'tex':
        TexLabel(Label)

    ReadProperty(File,Label,Properties)

#----------------------------------------------------------------------------
# Read Excess and Total Properties
#----------------------------------------------------------------------------
def ReadTotalProperty(File,line,OrdLabel,Properties):
    """Read selected total property."""

    #line  = FindLine(File,'BASED RESULTS')

    Label = line.split()

    if Label[1].count('CC') > 0:
        Label = ''.join([''.join(Label[:3]), ',', Label[-2], Label[-1]])
    else:
        Label = ''.join([Label[0], ',', Label[-2], Label[-1]])

    OrdLabel.append(Label)
    PrUnits['LabLen'].append(len(Label))
    if OutFormat == 'tex':
        TexLabel(Label)

    ReadProperty(File,Label,Properties)

#----------------------------------------------------------------------------
# Read property common routine
#----------------------------------------------------------------------------
def ReadProperty(File,Label,Properties):
    """Read selected property."""

    Properties[Label] = {}

    # Read dipole moment
    line = SkipLines(File,4)
    Mu = array(line.split(),dtype=float64)
    line = SkipLines(File,2)
    Dipole = float(line.split()[1])
    if abs(sqrt(dot(Mu,Mu))-Dipole) > 1e-5:
        print('Warning! Large difference in dipole moment'+Label)
    else:
        Dipole = sqrt(dot(Mu,Mu))

    # Read polarizability tensor
    line = SkipLines(File,2)
    X = SkipLines(File,1).split()[1:]
    Y = SkipLines(File,1).split()[1:]
    Z = SkipLines(File,1).split()[1:]
    Alpha = array([X,Y,Z],dtype=float64)

    # Read average polarizability
    line = SkipLines(File,2)
    AvgPolar = float(line.split()[1])
    if abs(trace(Alpha)/3.0-AvgPolar) > 1e-4:
        print('Warning! Large difference in <alpha>'+Label)
    else:
        AvgPolar = trace(Alpha)/3.0

    # Read polarizability anisotropy (Z is the rotation axis)
    line = SkipLines(File,1)
    AnzPolar = float(line.split()[3])
    if abs((Alpha[2,2]-Alpha[0,0])-AnzPolar) > 1e-4:
        print('Warning! Large difference in <alpha,ani>'+Label)
    else:
        AnzPolar = Alpha[2,2]-Alpha[0,0]

    # Read first hyperpolarizability tensor
    line = SkipLines(File,2)
    XX = SkipLines(File,1).split()[1:]
    YY = SkipLines(File,1).split()[1:]
    ZZ = SkipLines(File,1).split()[1:]
    Beta = array([XX,YY,ZZ],dtype=float64)

    # Read vector component of hyperpolarizability tensor
    # (Z is the permament dipole moment direction)
    line = SkipLines(File,2)
    VecFirstHyper = float(line.split()[3])
    if abs(VecFirstHyper-(3.0/5.0)*Beta.sum(axis=0)[2]) > 1e-4:
        print('Warning! Large difference in <beta,vec>'+Label)
    else:
        VecFirstHyper=(3.0/5.0)*Beta.sum(axis=0)[2]

    # Read second hyperpolarizability tensor
    line = SkipLines(File,2)
    XX = SkipLines(File,1).split()[1:]
    YY = SkipLines(File,1).split()[1:]
    ZZ = SkipLines(File,1).split()[1:]
    Gamma = array([XX,YY,ZZ],dtype=float64)

    # Read scalar component of second hyperpolarizability tensor
    # given by the isotropic average
    line = SkipLines(File,2)
    AvgSecHyper = float(line.split()[1])
    if abs(AvgSecHyper-(Gamma.trace()+2.0*(Gamma[0,1]+Gamma[1,2]+Gamma[0,2]))/5.0) > 1e-2:
        print('Warning! Large difference in <Gamma>'+Label)
        print((AvgSecHyper-(Gamma.trace()+2.0*(Gamma[0,1]+Gamma[1,2]+Gamma[0,2]))/5.0))
    else:
        AvgSecHyper = (Gamma.trace()+2.0*(Gamma[0,1]+Gamma[1,2]+Gamma[0,2]))/5.0

    Properties[Label]['Mu']    = around(Mu            * PropertyConFac['Mu'],    decimals=PrUnits['Round']['m'])
    Properties[Label]['|D|']   =  round(Dipole        * PropertyConFac['|D|'],            PrUnits['Round']['m'])
    Properties[Label]['Alpha'] = around(Alpha         * PropertyConFac['Alpha'], decimals=PrUnits['Round']['a'])
    Properties[Label]['<A>']   =  round(AvgPolar      * PropertyConFac['<A>'],            PrUnits['Round']['a'])
    Properties[Label]['<B>']   =  round(AnzPolar      * PropertyConFac['<B>'],            PrUnits['Round']['a'])
    Properties[Label]['Beta']  = around(Beta          * PropertyConFac['Beta'],  decimals=PrUnits['Round']['b'])
    Properties[Label]['B(Z)']  =  round(VecFirstHyper * PropertyConFac['B(Z)'],           PrUnits['Round']['b'])
    Properties[Label]['Gamma'] = around(Gamma         * PropertyConFac['Gamma'], decimals=PrUnits['Round']['g'])
    Properties[Label]['<G>']   =  round(AvgSecHyper   * PropertyConFac['<G>'],            PrUnits['Round']['g'])

#----------------------------------------------------------------------------
# Format property tables
#----------------------------------------------------------------------------
def FormatProperties(PropTables,Properties,OrdLabel):
    """
    Form tables of interaction induced properties
    
        Properties are stuck in a following dictionary:

        Properties{'File': 
                  {'Property type':
                  {'Energy term':
                  {'Property': Value }}}}
    """

    # List of files
    RunFiles = list(Properties.keys())
    if SortMode == 'float':
        RunFiles.sort(cmp=SortFiles)
    else:
        RunFiles.sort()

    # Property labels
    for Property in PropertyLabels:
        PropTables[Property] = {}

    # Labels of Total, Sum and Excess energy terms
    TotalTerms  = OrdLabel['TotPropLabel']
    SumTerms    = OrdLabel['SumPropLabel']
    ExcessTerms = OrdLabel['ExcPropLabel']
    IntTerms    = OrdLabel['IntPropLabel']

    if OutFormat == 'tex':
        TableHeader = '%'
    else:
        TableHeader = '#'

    # Comparison of Total, Sum and Excess properties
    for Property in PropertyLabels:

        # Set header
        Table = [[TableHeader]]
        Table[0].extend(TotalTerms)
        Table[0].extend(ExcessTerms)
        Table[0].extend(SumTerms)

        # Put the proper terms for all files
        for Column, RunFile in enumerate(RunFiles):

            Table.append([RunFile])

            for EnergyTerm in TotalTerms:
                Table[Column+1].append(Properties[RunFile]['Total'][EnergyTerm][Property])

            for EnergyTerm in ExcessTerms:
                Table[Column+1].append(Properties[RunFile]['Excess'][EnergyTerm][Property])

            for EnergyTerm in SumTerms:
                Table[Column+1].append(Properties[RunFile]['SumInteraction'][EnergyTerm][Property])

        # Form the property table for printout
        #Table = array(Table)
        PropTables[Property]['Excess'] = Table

    # Comparison of interaction induced properties
    for Property in PropertyLabels:

        # Set header
        Table = [[TableHeader]]
        Table[0].extend(IntTerms)

        # Put the proper terms for all files
        for Column, RunFile in enumerate(RunFiles):

            Table.append([RunFile])

            for EnergyTerm in IntTerms:
                Table[Column+1].append(Properties[RunFile]['Interaction'][EnergyTerm][Property])

        # Form the property table for printout
        #Table = array(Table)
        PropTables[Property]['Interaction'] = Table

#----------------------------------------------------------------------------
# Write Property Tables
#----------------------------------------------------------------------------
def WriteProperties(TitleLen,PropTables,Properties):
    """Save data to file."""

    if OutFormat == 'csv':
        Separator = ' ;'
        C         = '#'
        EndRow    = '\n'
        DataFileN = 'properties.csv'
        DataFileT = 'troperties.csv'
    elif OutFormat == 'tex':
        Separator = ' &'
        C         = '%'
        EndRow    = '\\\\\n'
        DataFileN = 'properties.tex'
        DataFileT = 'troperties.tex'
    else:
        Separator = '; '
        C         = '#'
        EndRow    = '\n'
        DataFileN = 'properties.txt'
        DataFileT = 'troperties.txt'

    DataFile=open(DataFileN,'w')
    TataFile=open(DataFileT,'w')

    Files = list(Properties.keys())
    Files.sort()

    if Relative:
        Format = '%14.1f' + Separator
        RoundP = 1

    LabLen = max(PrUnits['LabLen'])

    if OutFormat == 'tex':
        TitleFormat = '%-'+ str(TitleLen) + 's' + Separator
        LabelFormat = ' %' + str(LabLen) + 's ' + Separator
    else:
        TitleFormat = '"%-'+ str(TitleLen) + 's"' + Separator
        LabelFormat = '"%' + str(LabLen) + 's"' + Separator

    for l in list(PropertyFormats.keys()):
        PropertyFormats[l] = PropertyFormats[l].replace('%ln','%'+str(LabLen+2))+Separator

    # Write excess and total properties as well as sum of interaction
    # induced properties
    DataFile.write(C+' Total, Excess and Sum of Interaction Induced Properties\n'+C+'\n')
    TataFile.write(C+' Total, Excess and Sum of Interaction Induced Properties\n'+C+'\n')

    if OutFormat == 'csv':
        TableHeader = ''
    elif OutFormat == 'tex':
        TableHeader = '\\begin{tabular}{@{\extracolsep{\\fill}}l' + \
                      len(list(Properties[Files[0]]['Total'].keys()))  * ' r' + \
                      len(list(Properties[Files[0]]['Excess'].keys())) * ' r' + \
                      len(list(Properties[Files[0]]['SumInteraction'].keys())) * ' r' + '}\hline' + EndRow
        TableFooter = '\\end{tabular}' + EndRow
    else:
        TableHeader = '"' + C + 'Property'.center(TitleLen,'_') + '|' + \
                      'Total Properties'.center(len(list(Properties[Files[0]]['Total'].keys()))*(LabLen+4)-1,'_') + '|' + \
                      'Excess Properties'.center(len(list(Properties[Files[0]]['Excess'].keys()))*(LabLen+4)-1,'_') + '|' + \
                      'Sum of Interaction Properties'.center(len(list(Properties[Files[0]]['SumInteraction'].keys()))*(LabLen+4)-1,'_') + '|' + '"' + EndRow

    for Property in PropertyLabels:

        Table = PropTables[Property]['Excess']

        # Write table header
        DataFile.write(PropertyDescription[Property] + '\n' + C +'\n')
        DataFile.write(TableHeader)
        DataFile.write(TitleFormat % Table[0][0])

        for Column in range(len(Table[0])-1):
            Label = Table[0][Column+1]
            DataFile.write(LabelFormat % Label.rjust(LabLen))

        DataFile.write(EndRow)

        if OutFormat == 'tex':
            DataFile.write(TitleFormat % Table[0][0])
            for Column in range(len(Table[0])-1):
                Label = Table[0][Column+1]
                Label = TexLabel(Label)
                DataFile.write(LabelFormat % Label.rjust(LabLen))

            DataFile.write(EndRow)

        for Row in range(len(Table)-1):
            DataFile.write(TitleFormat % Table[Row+1][0].split()[1].replace('.log',''))

            for Column in range(len(Table[Row])-1):
                PrValue = Table[Row+1][Column+1]

                if type(PrValue) == ndarray:
                    PrValue = PrValue.flatten()[PropertyIndex[Property][0]]

                if str(type(PrValue)).count('float') > 0:
                    if Relative == 'first':
                        PrValue = round(100.0*PrValue / Table[Row+1][1],RoundP)
                        DataFile.write(Format % PrValue)
                    elif Relative == 'last':
                        PrValue = round(100.0*PrValue / Table[Row+1][-1],RoundP)
                        DataFile.write(Format % PrValue)
                    else:
                        DataFile.write(PropertyFormats[Property] % PrValue)
                else:
                    DataFile.write('%s' % PrValue.rjust(LabLen))
       
            DataFile.write(EndRow)
       
        if OutFormat == 'tex':
            DataFile.write(TableFooter)

        DataFile.write('\n')

        # Write Transposed properties
        TataFile.write(PropertyDescription[Property] + '\n' + C +'\n')
        TataFile.write(TableHeader)
        TataFile.write(TitleFormat % Table[0][0])

        for Row in range(len(Table)-1):
            TataFile.write(TitleFormat % Table[Row+1][0].split()[1].replace('.log',''))

        TataFile.write(EndRow)

        for Column in range(len(Table[0])-1):
            if OutFormat == 'tex':
                Label = Table[0][Column+1]
                Label = TexLabel(Label)
            else:
                Label = Table[0][Column+1]

            TataFile.write(LabelFormat % Label.rjust(LabLen))

            for Row in range(len(Table)-1):

                PrValue = Table[Row+1][Column+1]

                if type(PrValue) == ndarray:
                    PrValue = PrValue.flatten()[PropertyIndex[Property][0]]

                if str(type(PrValue)).count('float') > 0:
                    if Relative == 'first':
                        PrValue = round(100.0*PrValue / Table[Row+1][1],RoundP)
                        TataFile.write(Format % PrValue)
                    elif Relative == 'last':
                        PrValue = round(100.0*PrValue / Table[Row+1][-1],RoundP)
                        TataFile.write(Format % PrValue)
                    else:
                        TataFile.write(PropertyFormats[Property] % PrValue)
                else:
                    TataFile.write('%s' % PrValue.rjust(LabLen))

            if OutFormat == 'tex':
                Label = Table[0][Column+1]
                TataFile.write(' \\\\ ')
                TataFile.write(LabelFormat % ('% '+Label).rjust(LabLen))
                TataFile.write('\n')
            else:
                TataFile.write(EndRow)
       
        if OutFormat == 'tex':
            TataFile.write(TableFooter)

        TataFile.write('\n')

    # Write excess and total properties as well as sum of interaction
    # induced properties
    DataFile.write(C+' Interaction Induced Properties\n'+C+'\n')
    TataFile.write(C+' Interaction Induced Properties\n'+C+'\n')

    if OutFormat == 'csv':
        TableHeader = ''
    elif OutFormat == 'tex':
        TableHeader = '\\begin{tabular}{@{\extracolsep{\\fill}}l' + \
                      len(list(Properties[Files[0]]['Interaction'].keys())) * ' r' + EndRow
    else:
        TableHeader = '"' + C + 'Property'.center(TitleLen,'_') + '|' + \
                      'Finite Field Estimates of Interaction Induced Properties'.center(len(list(Properties[Files[0]]['Interaction'].keys()))*(LabLen+4)-1,'_') + '|' +'"' + EndRow

    for Property in PropertyLabels:

        Table = PropTables[Property]['Interaction']

        # Write table header
        DataFile.write(PropertyDescription[Property] + '\n'+C+'\n')
        DataFile.write(TableHeader)
        DataFile.write(TitleFormat % Table[0][0])

        for Column in range(len(Table[0])-1):
            Label = Table[0][Column+1]
            DataFile.write(LabelFormat % Label.rjust(LabLen))
       
        DataFile.write(EndRow)
       
        if OutFormat == 'tex':
            DataFile.write(TitleFormat % Table[0][0])
            for Column in range(len(Table[0])-1):
                Label = Table[0][Column+1]
                Label = TexLabel(Label)
                DataFile.write(LabelFormat % Label.rjust(LabLen))

            DataFile.write(EndRow)

        for Row in range(len(Table)-1):
            DataFile.write(TitleFormat % Table[Row+1][0].split()[1].replace('.log',''))

            for Column in range(len(Table[Row])-1):
                PrValue = Table[Row+1][Column+1]

                if type(PrValue) == ndarray:
                    PrValue = PrValue.flatten()[PropertyIndex[Property][0]]

                if str(type(PrValue)).count('float') > 0:
                    if Relative == 'first':
                        PrValue = round(100.0*PrValue / Table[Row+1][1],RoundP)
                        DataFile.write(Format % PrValue)
                    elif Relative == 'last':
                        PrValue = round(100.0*PrValue / Table[Row+1][-1],RoundP)
                        DataFile.write(Format % PrValue)
                    else:
                        DataFile.write(PropertyFormats[Property] % PrValue)
                else:
                    DataFile.write('%s' % PrValue.rjust(LabLen))

            DataFile.write(EndRow)

        if OutFormat == 'tex':
            DataFile.write(TableFooter)

        DataFile.write('\n')

        # Write Transposed properties
        TataFile.write(PropertyDescription[Property] + '\n' + C +'\n')
        TataFile.write(TableHeader)
        TataFile.write(TitleFormat % Table[0][0])

        for Row in range(len(Table)-1):
            TataFile.write(TitleFormat % Table[Row+1][0].split()[1].replace('.log',''))

        TataFile.write(EndRow)

        for Column in range(len(Table[0])-1):
            if OutFormat == 'tex':
                Label = Table[0][Column+1]
                Label = TexLabel(Label)
            else:
                Label = Table[0][Column+1]

            TataFile.write(LabelFormat % Label.rjust(LabLen))

            for Row in range(len(Table)-1):

                PrValue = Table[Row+1][Column+1]

                if type(PrValue) == ndarray:
                    PrValue = PrValue.flatten()[PropertyIndex[Property][0]]

                if str(type(PrValue)).count('float') > 0:
                    if Relative == 'first':
                        PrValue = round(100.0*PrValue / Table[Row+1][1],RoundP)
                        TataFile.write(Format % PrValue)
                    elif Relative == 'last':
                        PrValue = round(100.0*PrValue / Table[Row+1][-1],RoundP)
                        TataFile.write(Format % PrValue)
                    else:
                        TataFile.write(PropertyFormats[Property] % PrValue)
                else:
                    TataFile.write('%s' % PrValue.rjust(LabLen))

            if OutFormat == 'tex':
                Label = Table[0][Column+1]
                TataFile.write(' \\\\ ')
                TataFile.write(LabelFormat % ('% '+Label).rjust(LabLen))
                TataFile.write('\n')
            else:
                TataFile.write(EndRow)
       
        if OutFormat == 'tex':
            TataFile.write(TableFooter)

        TataFile.write('\n')

    # Close data file
    DataFile.close()
    TataFile.close()

#----------------------------------------------------------------------------
# Read Energies
#----------------------------------------------------------------------------
def ReadSubEnergies(File,CorrLabel,OrdLabel,Energies):
    """Read energies for this system."""

    line   = FindLine(File,'  INTERACTION ENERGY TERMS')
    Mer    = re.split('\D+',line)[1]
    ConfNo = re.split('\D+',line)[2]
    line   = SkipLines(File,4)

    while 1:
        line = File.readline()
        if line == '': break
        if line.find(20*'-') !=-1: break
        if line.find('(') !=-1:
            line=line.split()

            if len(line) == 3: 
                EnLabel = line[0]
                EnValue = line[1]
            if len(line) == 4:
                EnLabel = ' '.join(line[:2])
                EnValue = line[2]
            if MpLevel == 2 or (CcLevel and CcLevel.upper() != 'NONE'):
                EnLabel += CorrLabel
            if (EnLabel not in Energies) and OrdLabel.count(EnLabel) == 0 :
                OrdLabel.append(EnLabel)
                EnUnits['LabLen'].append(len(EnLabel))
                if OutFormat == 'tex':
                    TexLabel(EnLabel)
                Energies[EnLabel] = {}

            Energies[EnLabel][ConfNo] = EnValue

#----------------------------------------------------------------------------
# Read Total Energies
#----------------------------------------------------------------------------
def ReadTotEn(File,Energies,Field):
    """Read energies for this system."""

    line = FindLine(File,'TOTAL SCF ENERGIES')
    Energies[Field]['SCF'] = {}
    if Field==(0,0,0):
        line = SkipLines(File,2)
    else:
        line = FindLine(File,'FREE ENERGIES')
        line = SkipLines(File,1)

    while 1:
        line = File.readline()
        if line == '': break
        if line == '\n': break
        if line.find(10*'-') !=-1: break
        if line.find('(') !=-1:
            line   = re.split('\(|\)',line)
            Mer    = line[0].split('-')[0]
            ConfNo = int(line[1])
            Energies[Field]['SCF'][ConfNo] = float(line[2])

    if MpLevel == 2 or (CcLevel and CcLevel.upper() != 'NONE'):

        line = FindLine(File,'MP2 E(2) CORRECTIONS')
        Energies[Field]['MP2'] = {}
        line = SkipLines(File,2)

        while 1:
            line = File.readline()
            if line == '': break
            if line == '\n': break
            if line.find(10*'-') !=-1: break
            if line.find('(') !=-1:
                line   = re.split('\(|\)',line)
                Mer    = line[0].split('-')[0]
                ConfNo = int(line[1])
                Energies[Field]['MP2'][ConfNo] = Energies[Field]['SCF'][ConfNo] + float(line[2])

    if CcLevel and CcLevel.upper().count('CCSD(TQ') >= 1:

        line = FindLine(File,'CC CORRELATION ENERGY E(  CCSD(TQ))')
        Energies[Field]['CCSDTQ'] = {}
        line = SkipLines(File,2)

        while 1:
            line = File.readline()
            if line == '': break
            if line == '\n': break
            if line.find(10*'-') !=-1: break
            if line.find('(') !=-1:
                line   = re.split('\(|\)',line)
                Mer    = line[0].split('-')[0]
                ConfNo = int(line[1])
                Energies[Field]['CCSDTQ'][ConfNo] = Energies[Field]['SCF'][ConfNo] + float(line[2])

def ReadTotEnergies(File,out,Energies):
    """Read energies for this system."""

    File.seek(0)

    Field=(0,0,0)
    Energies[Field] = {}

    ReadTotEn(File,Energies,Field)

    while FiniteField:
        line = File.readline()
        if line == '': break
        if line.find('APPLIED FIELD') !=-1:
            line = line.split()
            # ... set field label ...
            Field = (float(line[-3]), float(line[-2]), float(line[-1]))
            if (Field not in Energies):
                Energies[Field]={}

            ReadTotEn(File,Energies,Field)

    for S in list(Energies[(0,0,0)]['SCF'].keys()):
        out.write( '# Subsystem %d \n' % S )

        row='# %21s' % 'Field'
        for L in list(Energies[(0,0,0)].keys()):
            row+='%26s' % L
        out.write( row+'\n' )

        for F in list(Energies.keys()):
            row='%7.4f %7.4f %7.4f' % F
            for L in list(Energies[(0,0,0)].keys()):
                row+='%26.15f' % Energies[F][L][S]
            out.write( row+'\n' )

#----------------------------------------------------------------------------
# Read Many Body Energy Terms
#----------------------------------------------------------------------------
def ReadMnbEnergies(File,OrdLabel,Energies,SeekIter):
    """Read energies for this system."""

    for i in range(SeekIter): 
        line = FindLine(File,'MANY BODY INTERACTION ENERGY TERMS')

    MnbLabel = line.split()[6]
    line     = SkipLines(File,4)

    while 1:
        line = File.readline()
        if line == '': break
        if line.find(20*'-') !=-1: break
        if line.find('(') !=-1:
            line=line.split()

            if len(line) == 3: 
                EnLabel = line[0]
                EnValue = line[1]
            if len(line) == 4:
                EnLabel = ' '.join(line[:2])
                EnValue = line[2]

            EnLabel += '(MNB)'

            if (EnLabel not in Energies) and OrdLabel.count(EnLabel) == 0 :
                OrdLabel.append(EnLabel)
                EnUnits['LabLen'].append(len(EnLabel))
                if OutFormat == 'tex':
                    TexLabel(EnLabel)
  
            Energies[EnLabel] = EnValue

def ReadFldEnergies(File,OrdLabel,Energies,SeekIter):
    """Read energies for this system."""

    for i in range(SeekIter): 
        line   = FindLine(File,'  INTERACTION ENERGY TERMS')

    Mer    = re.split('\D+',line)[1]
    ConfNo = re.split('\D+',line)[2]
    line   = SkipLines(File,4)

    while 1:
        line = File.readline()
        if line == '': break
        if line.find(20*'-') !=-1: break
        if line.find('(') !=-1:
            line=line.split()

            if len(line) == 3: 
                EnLabel = line[0]
                EnValue = line[1]
            if len(line) == 4:
                EnLabel = ' '.join(line[:2])
                EnValue = line[2]

            EnLabel += '(MNB)'

            if (EnLabel not in Energies) and OrdLabel.count(EnLabel) == 0 :
                OrdLabel.append(EnLabel)
                EnUnits['LabLen'].append(len(EnLabel))
                if OutFormat == 'tex':
                    TexLabel(EnLabel)
  
            Energies[EnLabel] = EnValue

#----------------------------------------------------------------------------
# Read Preamble
#----------------------------------------------------------------------------
def ReadPreamble(File,LogFile):
    """Read run title, and basic informations concerning the system."""

    global MpLevel, CcLevel, Monomers, Systems, ManyBody, FiniteField

    # Read filename and run title
    line  = FindLine(File,'RUN TITLE')
    line  = SkipLines(File,2)
    Title = "File: " + LogFile + " Run Title: " + line.strip()

    # Read filename and run title
    line    = FindLine(File,'CONTRL OPTIONS')
    line    = FindLine(File,'MPLEVL= ')
    MpLevel = int(line.split()[1])
    CcLevel = line.split()[5][1:]

    # Read number of monomers and subsystems
    line     = FindLine(File,'BODY COMPLEX')
    line     = re.split('\D+',line)[1:3]
    Monomers = int(line[1])
    Systems  = int(line[0])
    ManyBody = Systems > 3

    # Check for interaction induced properties run
    line        = FindLine(File,'FFEDS =')
    FiniteField = line.split()[11] == 'T'

    return Title

#----------------------------------------------------------------------------
# Format Energies
#----------------------------------------------------------------------------
def FormatSubEnergies(EnergyTables,ClusterTables,Energies,Labels):
    """
    Form tabularized energies.
    
        Energies are stuck in a following dictionary:
        Energies{'File':
                {(Field):
                {'Component':
                {'Subsystem': Value }}}}
    """

    # Labels of interaction energy components
    EnergyTerms = Labels['SubLabel']
    Field=(0,0,0)

    if OutFormat == 'tex':
        TableHeader = '%'
    else:
        TableHeader = '#'

    # List of files
    RunFiles = list(Energies.keys())
    RunFiles.sort()
    if len(RunFiles) > 1:
        Compare = True
    else:
        Compare = False

    # List of subsystems
    ComplexList = list(Energies[RunFiles[0]][Field]['DE(HF)'].keys())
    ComplexList.sort(key=lambda x: int(x))

    # Stack energies for comparison among subsystems
    for RunFile in RunFiles:

        # It's a bit difficult to sort therefore use OrdLabel
        #EnergyTerms = Energies[RunFile].keys()
        #EnergyTerms.sort()

        EnergyTable = [[TableHeader]]

        TempList = list(Energies[RunFile][Field]['DE(HF)'].keys())
        TempList.sort(key=lambda x: int(x))

        if ComplexList != TempList:
            ComplexList = TempList
            Compare = False

        EnergyTable[0].extend(ComplexList)

        for Column, EnergyTerm in enumerate(EnergyTerms):

            EnergyTable.append([EnergyTerm])

            for Cluster in ComplexList:
                if Cluster in Energies[RunFile][Field][EnergyTerm]:
                    EnergyTable[Column+1].append(Energies[RunFile][Field][EnergyTerm][Cluster])
                else:
                    EnergyTable[Column+1].append('-')

        #EnergyTable = transpose(array(EnergyTable,PyObject))
        EnergyTable = transpose(array(EnergyTable))
        EnergyTables[RunFile] = EnergyTable

    # Stack energies for comparison among files
    if Compare:

        for Cluster in ComplexList:
            
            ClusterTable = [[TableHeader]]
            ClusterTable[0].extend(RunFiles)

            for Column, EnergyTerm in enumerate(EnergyTerms):

                ClusterTable.append([EnergyTerm])

                for RunFile in RunFiles:
                    if Cluster in Energies[RunFile][Field][EnergyTerm]:
                        ClusterTable[Column+1].append(Energies[RunFile][Field][EnergyTerm][Cluster])
                    else:
                        ClusterTable[Column+1].append('-')

            #ClusterTable = transpose(array(ClusterTable,PyObject))
            ClusterTable = transpose(array(ClusterTable))
            ClusterTables[Cluster] = ClusterTable

#----------------------------------------------------------------------------
# Format Many Body Energies
#----------------------------------------------------------------------------
def FormatMnbEnergies(MbodyTables,Energies,Labels):
    """
    Form tabularized many body energies.
    
        Energies are stuck in a following dictionary:
        Energies{'File':
                {'Component': Value }}
    """

    # Labels of interaction energy components
    EnergyTerms = Labels['MnbLabel']
    Field=(0,0,0)

    # List of files
    RunFiles = list(Energies.keys())
    RunFiles.sort()

    if OutFormat == 'tex':
        TableHeader = '%'
    else:
        TableHeader = '#'

    # Stack energies for comparison among subsystems
    MbodyTable = [[TableHeader]]
    MbodyTable[0].extend(RunFiles)

    for Column, EnergyTerm in enumerate(EnergyTerms):

        MbodyTable.append([EnergyTerm])

        for RunFile in RunFiles:
            if EnergyTerm in Energies[RunFile][Field]:
                MbodyTable[Column+1].append(Energies[RunFile][Field][EnergyTerm])
            else:
                MbodyTable[Column+1].append('-')

    MbodyTable = transpose(array(MbodyTable))
    MbodyTables[Field] = MbodyTable

#----------------------------------------------------------------------------
# Format Field Energies
#----------------------------------------------------------------------------
def FormatFieldEnergies(FieldTables,Energies,Labels):
    """
    Form tabularized field energies.
        Energies are stuck in a following dictionary:
        Energies{'File':
                {'Component': Value }}
    """

    RunFiles = sorted(Energies.keys())

    FieldLabels = []
    RunFieldLabels = []
    EnLabels = []

    for RunFile in RunFiles:
        Fields = list(Energies[RunFile].keys())
        Fields = sorted(sorted(sorted(Fields,SortY),SortZ),SortX)
        #Fields = sorted( sorted( sorted(Fields, lambda a,b: cmp(abs(a[0]), abs(b[0]))), lambda a,b: cmp(abs(a[1]), abs(b[1]))), lambda a,b: cmp(abs(a[2]), abs(b[2])))
        for Field in Fields:
            RunFieldLabel=RunFile+'FIELD='+str(Field)
            if FieldLabels.count(Field) == 0:
                FieldLabels.append(Field)
            if RunFieldLabels.count(RunFieldLabel) == 0:
                RunFieldLabels.append(RunFieldLabel)
            
            for EnLabel in list(Energies[RunFile][Field].keys()):
                if EnLabels.count(EnLabel) == 0:
                    EnLabels.append(EnLabel)

    EnergyTerms = Labels['MnbLabel']

    if OutFormat == 'tex':
        TableHeader = '%'
    else:
        TableHeader = '#'
    
    # Stack energies for comparison among subsystems
    FieldTable = [[TableHeader]]
    FieldTable[0].extend(RunFieldLabels)
    
    for Column, EnergyTerm in enumerate(EnergyTerms):
    
        FieldTable.append([EnergyTerm])
    
        for RunFile in RunFiles:
            Fields = list(Energies[RunFile].keys())
            Fields = sorted(sorted(sorted(Fields,SortY),SortZ),SortX)
            for Field in Fields:
                if EnergyTerm in Energies[RunFile][Field]:
                    FieldTable[Column+1].append(Energies[RunFile][Field][EnergyTerm])
                else:
                    FieldTable[Column+1].append('-')
    
    FieldTable = transpose(array(FieldTable))
    FieldTables['MnbEn'] = FieldTable

#----------------------------------------------------------------------------
# Write Energies
#----------------------------------------------------------------------------
def WriteEnergies(TitleLen,EnergyTables,ClusterTables,MbodyTables,FieldTables):
    """Save data to file."""

    if OutFormat == 'csv':
        Separator = ' ;'
        C         = '#'
        EndRow    = '\n'
        DataFileN = 'energies.csv'
    elif OutFormat == 'tex':
        Separator = ' &'
        C         = '%'
        EndRow    = '\\\\\n'
        DataFileN = 'energies.tex'
    else:
        Separator = '; '
        C         = '#'
        EndRow    = '\n'
        DataFileN = 'energies.txt'

    DataFile=open(DataFileN,'w')

    Files = list(EnergyTables.keys())
    Files.sort()

    Clusters = list(ClusterTables.keys())
    Clusters.sort(key=lambda x: int(x))

    ConFac = EnUnits['ConFac']
    Format = EnUnits['Format']
    RoundE = EnUnits['Round']
    LabLen = max(EnUnits['LabLen'])

    if OutFormat == 'tex':
        TitleFormat = '%-' + str(TitleLen) + 's' + Separator
        LabelFormat = ' %' + str(LabLen) + 's ' + Separator
    else:
        TitleFormat = '"%-' + str(TitleLen) + 's"' + Separator
        LabelFormat = '"%' + str(LabLen) + 's"' + Separator

    ValueFormat = Format.replace('%ln','%'+str(LabLen+2))+Separator

    # Write interaction energies for all files
    for File in Files:

        Table = EnergyTables[File]

        # Write table header
        DataFile.write(C+' %s\n\n' % File)

        if OutFormat == 'tex':
            TableHeader = '\\begin{tabular}{@{\extracolsep{\\fill}}l' + (len(Table)-1) * ' r' + '}\hline' + EndRow
            TableFooter = '\\end{tabular}' + EndRow
            DataFile.write(TableHeader)

        DataFile.write(TitleFormat % Table[0][0])
        for Column in range(len(Table[0])-1):
            Label = Table[0][Column+1].replace('(CORR)','')
            DataFile.write(LabelFormat % Label.rjust(LabLen))

        DataFile.write(EndRow)

        if OutFormat == 'tex':
            DataFile.write(TitleFormat % Table[0][0])
            for Column in range(len(Table[0])-1):
                Label = Table[0][Column+1].replace('(CORR)','')
                Label = TexLabel(Label)
                DataFile.write(LabelFormat % Label.rjust(LabLen))

            DataFile.write(EndRow)

        for Row in range(len(Table)-1):
            DataFile.write(TitleFormat % Table[Row+1][0])

            for Column in range(len(Table[Row])-1):
                EnValue = Table[Row+1][Column+1]

                try:
                    EnValue = round(float(EnValue)*ConFac,RoundE)
                    DataFile.write(ValueFormat % EnValue)
                except TypeError:
                    EnValue = float(EnValue)*ConFac
                    DataFile.write(ValueFormat % EnValue)
                except ValueError:
                    DataFile.write(LabelFormat % EnValue.rjust(LabLen))

            DataFile.write(EndRow)

        if OutFormat == 'tex':
            DataFile.write(TableFooter)
        DataFile.write('\n')

    # This time compare respective components in files
    for Cluster in Clusters:
        Table = ClusterTables[Cluster]

        # Write table header
        DataFile.write(C+' Subsystem No: %s\n\n' % Cluster)

        if OutFormat == 'tex':
            TableHeader = '\\begin{tabular}{@{\extracolsep{\\fill}}l' + (len(Table)-1) * ' r' + '}\hline' + EndRow
            TableFooter = '\\end{tabular}' + EndRow
            DataFile.write(TableHeader)

        DataFile.write(TitleFormat % Table[0][0])

        for Column in range(len(Table[0])-1):
            Label = Table[0][Column+1].replace('(CORR)','')
            DataFile.write(LabelFormat % Label.rjust(LabLen))

        DataFile.write(EndRow)

        if OutFormat == 'tex':
            DataFile.write(TitleFormat % Table[0][0])
            for Column in range(len(Table[0])-1):
                Label = Table[0][Column+1].replace('(CORR)','')
                Label = TexLabel(Label)
                DataFile.write(LabelFormat % Label.rjust(LabLen))

            DataFile.write(EndRow)

        for Row in range(len(Table)-1):
            DataFile.write(TitleFormat % Table[Row+1][0].split()[1].replace('.log',''))

            for Column in range(len(Table[Row])-1):
                EnValue = Table[Row+1][Column+1]
                try:
                    EnValue = round(float(EnValue)*ConFac,RoundE)
                    DataFile.write(ValueFormat % EnValue)
                except TypeError:
                    EnValue = float(EnValue)*ConFac
                    DataFile.write(ValueFormat % EnValue)
                except ValueError:
                    DataFile.write(LabelFormat % EnValue.rjust(LabLen))

            DataFile.write(EndRow)

        if OutFormat == 'tex':
            DataFile.write(TableFooter)
        DataFile.write('\n')

    # Write many body energy components
    if OutFormat == 'tex':
        DataFile.write(TableHeader)

    if ManyBody:

        Table = MbodyTables[(0,0,0)]

        # Write table header
        DataFile.write(C+' Many-body energy terms for selected systems\n\n')

        if OutFormat == 'tex':
            TableHeader = '\\begin{tabular}{@{\extracolsep{\\fill}}l' + (len(Table)-1) * ' r' + '}\hline' + EndRow
            TableFooter = '\\end{tabular}' + EndRow
            DataFile.write(TableHeader)

        DataFile.write(TitleFormat % Table[0][0])
       
        for Column in range(len(Table[0])-1):
            Label = Table[0][Column+1].replace('(MNB)','')
            Label = Label.replace(' ','')
            DataFile.write(LabelFormat % Label.rjust(LabLen))
        
        DataFile.write(EndRow)
       
        if OutFormat == 'tex':
            DataFile.write(TitleFormat % Table[0][0])
            for Column in range(len(Table[0])-1):
                Label = Table[0][Column+1].replace('(MNB)','')
                Label = Label.replace(' ','')
                Label = TexLabel(Label)
                DataFile.write(LabelFormat % Label.rjust(LabLen))

            DataFile.write(EndRow)

        for Row in range(len(Table)-1):
            DataFile.write(TitleFormat % Table[Row+1][0].split()[1].replace('.log',''))
        
            for Column in range(len(Table[Row])-1):
                EnValue = Table[Row+1][Column+1]
                try:
                    EnValue = round(float(EnValue)*ConFac,RoundE)
                    DataFile.write(ValueFormat % EnValue)
                except TypeError:
                    EnValue = float(EnValue)*ConFac
                    DataFile.write(ValueFormat % EnValue)
                except ValueError:
                    DataFile.write(LabelFormat % EnValue.rjust(LabLen))
        
            DataFile.write(EndRow)
        
        if OutFormat == 'tex':
            DataFile.write(TableFooter)
        DataFile.write('\n')

    # Write many body field energy components
    if OutFormat == 'tex':
        DataFile.write(TableHeader)

    if FiniteField:

        Table = FieldTables['MnbEn']
        #savetxt('test.out', Table, fmt='%'+str(TitleLen)+'s'+(len(Table[0])-1)*'%23s')

        # Write table header
        DataFile.write(C+' Many-body energy terms for selected fields\n\n')

        if OutFormat == 'tex':
            TableHeader = '\\begin{tabular}{@{\extracolsep{\\fill}}l' + (len(Table)-1) * ' r' + '}\hline' + EndRow
            TableFooter = '\\end{tabular}' + EndRow
            DataFile.write(TableHeader)

        DataFile.write(TitleFormat % Table[0][0])
       
        for Column in range(len(Table[0])-1):
            Label = Table[0][Column+1].replace('(MNB)','')
            Label = Label.replace(' ','')
            DataFile.write(LabelFormat % Label.rjust(LabLen))
        
        DataFile.write(EndRow)
       
        if OutFormat == 'tex':
            DataFile.write(TitleFormat % Table[0][0])
            for Column in range(len(Table[0])-1):
                Label = Table[0][Column+1].replace('(MNB)','')
                Label = Label.replace(' ','')
                Label = TexLabel(Label)
                DataFile.write(LabelFormat % Label.rjust(LabLen))
        
            DataFile.write(EndRow)

        for Row in range(len(Table)-1):
            FieldString='%7.4f, %7.4f, %7.4f' % tuple(asarray(mat( Table[Row+1][0].split('FIELD=')[1] )).tolist()[0])
            DataFile.write(TitleFormat % FieldString)
        
            for Column in range(len(Table[Row])-1):
                EnValue = Table[Row+1][Column+1]
                try:
                    EnValue = round(float(EnValue)*ConFac,RoundE)
                    DataFile.write(ValueFormat % EnValue)
                except TypeError:
                    EnValue = float(EnValue)*ConFac
                    DataFile.write(ValueFormat % EnValue)
                except ValueError:
                    DataFile.write(LabelFormat % EnValue.rjust(LabLen))
        
            DataFile.write(TitleFormat % C+' '+Table[Row+1][0].split('FIELD=')[0])
            DataFile.write(EndRow)
        
        if OutFormat == 'tex':
            DataFile.write(TableFooter)
        DataFile.write('\n')

    # Close data file
    DataFile.close()

#----------------------------------------------------------------------------
# Utilities
#----------------------------------------------------------------------------
def SkipLines(File,n):
    """Read n lines from file f."""

    for i in range(n):
        line = File.readline()
        if line == '' : break

    return line

def FindLine(File,pattern):
    """Read lines until pattern matches."""

    while 1:
        line = File.readline()
        if line.find(pattern) !=-1 : break
        if line == '' : break

    return line

def SortFiles(x, y):

    a=re.compile('\d+\.\d+').findall(x)
    b=re.compile('\d+\.\d+').findall(y)

    if a and b:
        for i in range(len(a)):
            a[i]=float(a[i])
        for i in range(len(b)):
            b[i]=float(b[i])
        #a=float(a.group())
        #b=float(b.group())
    else:
        a=x
        b=y

    if a>b:
        return 1
    elif a==b:
        return 0
    else: # a<b
        return -1

def SortX(x, y):

    a=x[0]
    b=y[0]

    if x[0] == 0 and x[1] == 0 and x[2] == 0:
        return -1
    else:
        if a == 0:
            return -1
        elif b == 0:
            return 1
        else:
            if a<0 and b<0:
                if a>b:
                    return 1
                elif a==b:
                    return 0
                else: # a<b
                    return -1
            else:
                if a>b:
                    return -1
                elif a==b:
                    return 0
                else: # a<b
                    return 1

def SortY(x, y):

    a=x[1]
    b=y[1]

    if x[0] == 0 and x[1] == 0 and x[2] == 0:
        return 0
    else:
        if a == 0:
            return -1
        elif b == 0:
            return 1
        else:
            if a<0 and b<0:
                if a>b:
                    return 1
                elif a==b:
                    return 0
                else: # a<b
                    return -1
            else:
                if a>b:
                    return -1
                elif a==b:
                    return 0
                else: # a<b
                    return 1

def SortZ(x, y):

    a=x[2]
    b=y[2]

    if x[0] == 0 and x[1] == 0 and x[2] == 0:
        return 0
    else:
        if a == 0:
            return -1
        elif b == 0:
            return 1
        else:
            if a<0 and b<0:
                if a>b:
                    return -1
                elif a==b:
                    return 0
                else: # a<b
                    return 1
            else:
                if a>b:
                    return 1
                elif a==b:
                    return 0
                else: # a<b
                    return -1

#----------------------------------------------------------------------------
# Units, conversion factors and formats
#----------------------------------------------------------------------------
def TexLabel(label):

    label = label.replace('(CORR)','')
    label = label.replace('(MNB)','')

    if re.search('[A-Z]+\((\w+|\w-\w+)\(\w+\)\)', label):
        label = re.split('\(|\)',label)
        if label[0] == 'E' or label[0] == 'G':
            label[0] = label[0].replace('E','\\epsilon')
            label[0] = label[0].replace('G','\\epsilon')
            label    = ''.join([label[0], '^{\\rm ', label[1], '(', label[2], ')', '}'])
        if label[0] == 'DE' or label[0] == 'DG':
            label[0] = label[0].replace('DE','\\Delta E')
            label[0] = label[0].replace('DG','\\Delta E')
            label    = ''.join([label[0], '^{\\rm ', label[1], '(', label[2], ')', '}'])
    elif re.search('[A-Z]+\(\w+\)', label):
        label = re.split('\(|\)',label)
        if label[0] == 'E' or label[0] == 'G':
            label[0] = label[0].replace('E','\\epsilon')
            label[0] = label[0].replace('G','\\epsilon')
        if label[0] == 'DE' or label[0] == 'DG':
            label[0] = label[0].replace('DE','\\Delta E')
            label[0] = label[0].replace('DG','\\Delta E')
        label = ''.join([label[0], '^{\\rm ', label[1], '}'])
    elif re.search('[A-Z]+\(\w+,\w+\)', label):
        label = re.split('\(|\)',label)
        if label[0] == 'E' or label[0] == 'G':
            label[0] = label[0].replace('E','\\epsilon')
            label[0] = label[0].replace('G','\\epsilon')
        if label[0] == 'DE' or label[0] == 'DG':
            label[0] = label[0].replace('DE','\\Delta E')
            label[0] = label[0].replace('DG','\\Delta E')
        label = ''.join([label[0], '_{\\rm ', label[1].split(',')[0], '}', '^{\\rm ', label[1].split(',')[1], '}'])
    else:
        print("Warning! Unknown label")

    EnUnits['LabLen'].append(len(label)+2)
    PrUnits['LabLen'].append(len(label)+2)

    return '$'+label+'$'

def SetOutFormat(arg):
    """Set output format."""
    if arg.lower() == 'csv' or arg.lower() == 'tex':
        return arg.lower()
    else:
        return 'txt'

def SetEnUnits(Units):
    """Set energy units and appropriate formats."""

    global EnUnits

    EnUnits = {}

    if Units.lower() == 'kcal':
        EnUnits['ConFac'] = 627.509541
        EnUnits['Format'] = '%ln.3f'
        EnUnits['Round']  = 3
        EnUnits['LabLen'] = [14]
    elif Units.lower() == 'kj':
        EnUnits['ConFac'] = 627.509541 * 4.1840
        EnUnits['Format'] = '%ln.3f'
        EnUnits['Round']  = 3
        EnUnits['LabLen'] = [14]
    elif Units.lower() == 'au':
        EnUnits['ConFac'] = 1.0
        EnUnits['Format'] = '%ln.12e'
        EnUnits['Round']  = ''
        EnUnits['LabLen'] = [23]
    elif Units.lower() == 'mh':
        EnUnits['ConFac'] = 1.0e3
        EnUnits['Format'] = '%ln.3f'
        EnUnits['Round']  = 3
        EnUnits['LabLen'] = [14]
    elif Units.lower() == 'mev':
        EnUnits['ConFac'] = 27211.3845
        EnUnits['Format'] = '%ln.3f'
        EnUnits['Round']  = 3
        EnUnits['LabLen'] = [14]
    else:
        Usage()
        sys.exit(2)


def SetPrUnits(Units):
    """Set units of electric properties and the respective formats."""

    global PrUnits, PropertyLabels, PropertyConFac, PropertyFormats, PropertyIndex, PropertyDescription

    PrUnits = {}

    if Units.lower() == 'au':
        PrUnits['Mu']    = 1.0
        PrUnits['Alpha'] = 1.0
        PrUnits['Beta']  = 1.0
        PrUnits['Gamma'] = 1.0
        PrUnits['Format'] = {'m':'%ln.4f', 'a':'%ln.3f', 'b':'%ln.2f', 'g':'%ln.1f'}
        PrUnits['Round']  = {'m':4,        'a':3,        'b':2,        'g':1}
        PrUnits['LabLen'] = [10]
    elif Units.lower() == 'mau':
        PrUnits['Mu']    = 1.0e3
        PrUnits['Alpha'] = 1.0e3
        PrUnits['Beta']  = 1.0e3
        PrUnits['Gamma'] = 1.0e3
        PrUnits['Format'] = {'m':'%ln.1f', 'a':'%ln.1f', 'b':'%ln.1f', 'g':'%ln.1f'}
        PrUnits['Round']  = {'m':1,        'a':1,        'b':1,        'g':1}
        PrUnits['LabLen'] = [10]
    elif Units.lower() == 'si':
        PrUnits['Mu']    = 8.478358e-30 # C m
        PrUnits['Alpha'] = 1.648778e-41 # C^2 m^2 J^-1
        PrUnits['Beta']  = 3.206361e-53 # C^3 m^3 J^-2
        PrUnits['Gamma'] = 6.235377e-65 # C^4 m^4 J^-3
        PrUnits['Format'] = {'m':'%ln.5e', 'a':'%ln.5e', 'b':'%ln.5e', 'g':'%ln.5e'}
        PrUnits['Round']  = {'m':3,        'a':3,        'b':3,        'g':3}
        PrUnits['LabLen'] = [10]
    elif Units.lower() == 'asi':
        PrUnits['Mu']    = 8.4784e-30 # C m
        PrUnits['Alpha'] = 1.8621e-30 # m^3
        PrUnits['Beta']  = 3.6213e-42 # m^4 V^-1
        PrUnits['Gamma'] = 7.0423e-54 # m^5 V^-2
        PrUnits['Format'] = {'m':'%ln.4e', 'a':'%ln.4e', 'b':'%ln.4e', 'g':'%ln.4e'}
        PrUnits['Round']  = {'m':3,        'a':3,        'b':3,        'g':3}
        PrUnits['LabLen'] = [10]
    elif Units.lower() == 'esu':
        PrUnits['Mu']    = 2.5418e-18 # statvolt cm^2
        PrUnits['Alpha'] = 1.4817e-25 # cm^3
        PrUnits['Beta']  = 8.6392e-33 # statvolt^-1 cm^4
        PrUnits['Gamma'] = 5.0367e-40 # statvolt^-2 cm^5
        PrUnits['Format'] = {'m':'%ln.4e', 'a':'%ln.4e', 'b':'%ln.4e', 'g':'%ln.4e'}
        PrUnits['Round']  = {'m':40,        'a':40,        'b':40,        'g':40}
        PrUnits['LabLen'] = [10]
    else:
        Usage()
        sys.exit(2)

    # Property labels - make sure they are the same as in the
    # ReadProperty() routine
    PropertyLabels = ['Mu', '|D|', 'Alpha', '<A>', '<B>', 'Beta', 'B(Z)', 'Gamma', '<G>']

    # Property conversion factors
    PropertyConFac = {
        'Mu'    : PrUnits['Mu'],
        '|D|'   : PrUnits['Mu'],
        'Alpha' : PrUnits['Alpha'],
        '<A>'   : PrUnits['Alpha'],
        '<B>'   : PrUnits['Alpha'],
        'Beta'  : PrUnits['Beta'],
        'B(Z)'  : PrUnits['Beta'],
        'Gamma' : PrUnits['Gamma'],
        '<G>'   : PrUnits['Gamma'] }

    # Property formats
    PropertyFormats = {
        'Mu'    : PrUnits['Format']['m'],
        '|D|'   : PrUnits['Format']['m'],
        'Alpha' : PrUnits['Format']['a'],
        '<A>'   : PrUnits['Format']['a'],
        '<B>'   : PrUnits['Format']['a'],
        'Beta'  : PrUnits['Format']['b'],
        'B(Z)'  : PrUnits['Format']['b'],
        'Gamma' : PrUnits['Format']['g'],
        '<G>'   : PrUnits['Format']['g'] }

    # Selected tensor elements
    #               Mu:     Alpha:     Beta:         Gamma:
    # [[0, 1, 2],   x y z   xx xy xz   xxx xxy xxz   xxxx xxyy xxzz
    #  [3, 4, 5],           yx yy yz   yyx yyy yyz   yyxx yyyy yyzz
    #  [6, 7, 8]]           zx zy zz   zzx zzy zzz   zzxx zzyy zzzz
    # 
    PropertyIndex = {
        'Mu'    : [2, 'x'   ,'y'   ,'z'  ],
        'Alpha' : [8, 'xx'  ,'xy'  ,'xz'  , 'yx'  ,'yy'  ,'yz'  , 'zx'  ,'zy'  ,'zz'  ],
        'Beta'  : [8, 'xxx' ,'xxy' ,'xxz' , 'yyx' ,'yyy' ,'yyz' , 'zzx' ,'zzy' ,'zzz' ],
        'Gamma' : [8, 'xxxx','xxyy','xxzz', 'yyxx','yyyy','yyzz', 'zzxx','zzyy','zzzz'] }

    # Descriptions of properties
    PropertyDescription = {
        'Mu'    : '# Dipole Moment Vector i=' + \
                     PropertyIndex['Mu'][PropertyIndex['Mu'][0]+1],
        '|D|'   : '# Total Dipole Moment',
        'Alpha' : '# Polarizability Tensor i=' + \
                     PropertyIndex['Alpha'][PropertyIndex['Alpha'][0]+1],
        '<A>'   : '# Isotropic Polatizability',
        '<B>'   : '# Anisotropy of Polarizability (Z-axis is the rotation axis)',
        'Beta'  : '# First hyperpolarizability tensor i=' + \
                     PropertyIndex['Beta'][PropertyIndex['Beta'][0]+1],
        'B(Z)'  : '# Vector component of hyperpolarizability tensor (Z is the permament dipole moment direction)',
        'Gamma' : '# Second hyperpolarizability tensor i=' + \
                     PropertyIndex['Gamma'][PropertyIndex['Gamma'][0]+1],
        '<G>'   : '# Scalar component of second hyperpolarizability tensor given by the isotropic average' }

#----------------------------------------------------------------------------
# Main routine
#----------------------------------------------------------------------------
if __name__ == "__main__": Main(sys.argv[1:])

