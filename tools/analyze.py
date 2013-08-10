#!/usr/local/bin/env python

#=============================================================================================
# Analyze datafiles produced by YANK.
#=============================================================================================

#=============================================================================================
# REQUIREMENTS
#
# The netcdf4-python module is now used to provide netCDF v4 support:
# http://code.google.com/p/netcdf4-python/
#
# This requires NetCDF with version 4 and multithreading support, as well as HDF5.
#=============================================================================================

#=============================================================================================
# TODO
#=============================================================================================

#=============================================================================================
# CHAGELOG
#=============================================================================================

#=============================================================================================
# VERSION CONTROL INFORMATION
#=============================================================================================

#=============================================================================================
# IMPORTS
#=============================================================================================

import os
import os.path
import sys
import math

import numpy

import netCDF4 as netcdf # netcdf4-python

from pymbar import MBAR # multistate Bennett acceptance ratio
import timeseries # for statistical inefficiency analysis

import simtk.unit as units

from yank.analysis import *

#=============================================================================================
# MAIN
#=============================================================================================

if __name__ == '__main__':    

    #molecules = ['1-methylpyrrole', '12-dichlorobenzene', '2-fluorobenzaldehyde', '23-benzofuran', 'benzene', 'benzenedithiol', 'ethylbenzene', 'indene', 'indole', 'isobutylbenzene', 'n-butylbenzene', 'n-methylaniline', 'n-propylbenzene', 'o-xylene', 'p-xylene', 'phenol', 'thieno_23c_pyridine', 'toluene']
    #base_directory = '/Users/yank/amber-gbsa'


    # DEBUG: ANALYSIS PATH IS HARD-CODED FOR NOW

    source_directory = "examples/p-xylene"
    #source_directory = "examples/benzene-toluene"
    #source_directory = "../test-systems/T4-lysozyme-L99A/amber-gbsa/amber-gbsa/p-xylene"
    #source_directory = "../test-systems/T4-lysozyme-L99A/amber-gbsa/amber-gbsa/1-methylpyrrole"

    # Storage for different phases.
    data = dict()

    phases = ['vacuum', 'solvent', 'complex']

    # Process each netcdf file.
    for phase in phases:
        # Construct full path to NetCDF file.
        fullpath = os.path.join(source_directory, phase + '.nc')

        # Skip if the file doesn't exist.
        if (not os.path.exists(fullpath)): continue

        # Open NetCDF file for reading.
        print "Opening NetCDF trajectory file '%(fullpath)s' for reading..." % vars()
        ncfile = netcdf.Dataset(fullpath, 'r')

        # DEBUG
        print "dimensions:"
        for dimension_name in ncfile.dimensions.keys():
            print "%16s %8d" % (dimension_name, len(ncfile.dimensions[dimension_name]))
    
        # Read dimensions.
        niterations = ncfile.variables['positions'].shape[0]
        nstates = ncfile.variables['positions'].shape[1]
        natoms = ncfile.variables['positions'].shape[2]
        print "Read %(niterations)d iterations, %(nstates)d states" % vars()

#        # Compute torsion trajectories
#        if phase in ['complex', 'receptor']:
#            print "Computing torsions..."
#            compute_torsion_trajectories(ncfile, os.path.join(source_directory, phase + ".val111"))

#        # Write out replica trajectories
#        print "Writing replica trajectories...\n"
#        title = 'Source %(source_directory)s phase %(phase)s' % vars()        
#        write_netcdf_replica_trajectories(source_directory, phase, title, ncfile)

        # Read reference PDB file.
        if phase in ['vacuum', 'solvent']:
            reference_pdb_filename = os.path.join(source_directory, "ligand.pdb")
        else:
            reference_pdb_filename = os.path.join(source_directory, "complex.pdb")
        atoms = read_pdb(reference_pdb_filename)

        # Write replica trajectories.
        title = 'title'
        write_pdb_replica_trajectories(reference_pdb_filename, source_directory, phase, title, ncfile, trajectory_by_state=False)

        # Check to make sure no self-energies go nan.
        check_energies(ncfile, atoms)

        # Check to make sure no positions are nan
        check_positions(ncfile)

        # Choose number of samples to discard to equilibration
        #nequil = 50
        #if phase == 'complex':
        #    nequil = 2000 # discard 2 ns of complex simulations
        u_n = extract_u_n(ncfile)
        [nequil, g_t, Neff_max] = detect_equilibration(u_n)
        print [nequil, Neff_max]
 
        # Examine acceptance probabilities.
        show_mixing_statistics(ncfile, cutoff=0.05, nequil=nequil)

        # Estimate free energies.
        (Deltaf_ij, dDeltaf_ij) = estimate_free_energies(ncfile, ndiscard = nequil)
    
        # Estimate average enthalpies
        (DeltaH_i, dDeltaH_i) = estimate_enthalpies(ncfile, ndiscard = nequil)
    
        # Accumulate free energy differences
        entry = dict()
        entry['DeltaF'] = Deltaf_ij[0,nstates-1] 
        entry['dDeltaF'] = dDeltaf_ij[0,nstates-1]
        entry['DeltaH'] = DeltaH_i[nstates-1] - DeltaH_i[0]
        entry['dDeltaH'] = numpy.sqrt(dDeltaH_i[0]**2 + dDeltaH_i[nstates-1]**2)
        data[phase] = entry

        # Get temperatures.
        ncvar = ncfile.groups['thermodynamic_states'].variables['temperatures']
        temperature = ncvar[0] * units.kelvin
        kT = kB * temperature

        # Close input NetCDF file.
        ncfile.close()
    
    # Compute hydration free energy (free energy of transfer from vacuum to water)
    DeltaF = data['vacuum']['DeltaF'] - data['solvent']['DeltaF']
    dDeltaF = numpy.sqrt(data['vacuum']['dDeltaF']**2 + data['solvent']['dDeltaF']**2)
    print "Hydration free energy: %.3f +- %.3f kT (%.3f +- %.3f kcal/mol)" % (DeltaF, dDeltaF, DeltaF * kT / units.kilocalories_per_mole, dDeltaF * kT / units.kilocalories_per_mole)

    # Compute enthalpy of transfer from vacuum to water
    DeltaH = data['vacuum']['DeltaH'] - data['solvent']['DeltaH']
    dDeltaH = numpy.sqrt(data['vacuum']['dDeltaH']**2 + data['solvent']['dDeltaH']**2)
    print "Enthalpy of hydration: %.3f +- %.3f kT (%.3f +- %.3f kcal/mol)" % (DeltaH, dDeltaH, DeltaH * kT / units.kilocalories_per_mole, dDeltaH * kT / units.kilocalories_per_mole)

    # Read standard state correction free energy.
    DeltaF_restraints = 0.0
    phase = 'complex'
    fullpath = os.path.join(source_directory, phase + '.nc')
    ncfile = netcdf.Dataset(fullpath, 'r')
    DeltaF_restraints = ncfile.groups['metadata'].variables['standard_state_correction'][0]
    ncfile.close()

    # Compute binding free energy (free energy of transfer from vacuum to water)
    DeltaF = data['solvent']['DeltaF'] - DeltaF_restraints - data['complex']['DeltaF']
    dDeltaF = numpy.sqrt(data['solvent']['dDeltaF']**2 + data['complex']['dDeltaF']**2)
    print ""
    print "Binding free energy : %16.3f +- %.3f kT (%16.3f +- %.3f kcal/mol)" % (DeltaF, dDeltaF, DeltaF * kT / units.kilocalories_per_mole, dDeltaF * kT / units.kilocalories_per_mole)
    print ""
    print "DeltaG vacuum       : %16.3f +- %.3f kT" % (data['vacuum']['DeltaF'], data['vacuum']['dDeltaF'])
    print "DeltaG solvent      : %16.3f +- %.3f kT" % (data['solvent']['DeltaF'], data['solvent']['dDeltaF'])
    print "DeltaG complex      : %16.3f +- %.3f kT" % (data['complex']['DeltaF'], data['complex']['dDeltaF'])
    print "DeltaG restraint    : %16.3f          kT" % DeltaF_restraints
    print ""

    # Compute binding enthalpy
    DeltaH = data['solvent']['DeltaH'] - DeltaF_restraints - data['complex']['DeltaH'] 
    dDeltaH = numpy.sqrt(data['solvent']['dDeltaH']**2 + data['complex']['dDeltaH']**2)
    print "Binding enthalpy    : %16.3f +- %.3f kT (%16.3f +- %.3f kcal/mol)" % (DeltaH, dDeltaH, DeltaH * kT / units.kilocalories_per_mole, dDeltaH * kT / units.kilocalories_per_mole)

