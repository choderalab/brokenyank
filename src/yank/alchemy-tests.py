#!/usr/bin/python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Alchemical factory for free energy calculations that operates directly on OpenMM swig System objects.

DESCRIPTION

This module contains enumerative factories for generating alchemically-modified System objects
usable for the calculation of free energy differences of hydration or ligand binding.

The code in this module operates directly on OpenMM Swig-wrapped System objects for efficiency.

EXAMPLES

COPYRIGHT

@author John D. Chodera <jchodera@gmail.com>

All code in this repository is released under the GNU General Public License.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

TODO

* Can we store serialized form of Force objects so that we can save time in reconstituting
  Force objects when we make copies?  We can even manipulate the XML representation directly.
* Allow protocols to automatically be resized to arbitrary number of states, to 
  allow number of states to be enlarged to be an integral multiple of number of GPUs.
* Add GBVI support to AlchemicalFactory.
* Add analytical dispersion correction to softcore Lennard-Jones, or find some other
  way to deal with it (such as simply omitting it from lambda < 1 states).
* Deep copy Force objects that don't need to be modified instead of using explicit 
  handling routines to copy data.  Eventually replace with removeForce once implemented?
* Can alchemically-modified System objects share unmodified Force objects to avoid overhead
  of duplicating Forces that are not modified?

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
import numpy
import copy
import time

import simtk.openmm as openmm

from sets import Set

from alchemy import AlchemicalState, AbsoluteAlchemicalFactory

#=============================================================================================
# MAIN AND UNIT TESTS
#=============================================================================================

def testAlchemicalFactory(reference_system, coordinates, receptor_atoms, ligand_atoms, platform_name='CUDA', annihilateElectrostatics=True, annihilateLennardJones=False):
    """
    Compare energies of reference system and fully-interacting alchemically modified system.

    ARGUMENTS
    
    reference_system (simtk.openmm.System) - the reference System object to compare with
    coordinates - the coordinates to assess energetics for
    receptor_atoms (list of int) - the list of receptor atoms 
    ligand_atoms (list of int) - the list of ligand atoms to alchemically modify

    """

    import simtk.unit as units
    import simtk.openmm as openmm
    import time

    # Create a factory to produce alchemical intermediates.
    print "Creating alchemical factory..."
    initial_time = time.time()
    factory = AbsoluteAlchemicalFactory(reference_system, ligand_atoms=ligand_atoms)
    final_time = time.time()
    elapsed_time = final_time - initial_time
    print "AbsoluteAlchemicalFactory initialization took %.3f s" % elapsed_time

    # Create an alchemically-perturbed state corresponding to nearly fully-interacting.
    # NOTE: We use a lambda slightly smaller than 1.0 because the AlchemicalFactory does not use Custom*Force softcore versions if lambda = 1.0 identically.
    lambda_value = 1.0 - 1.0e-6
    alchemical_state = AlchemicalState(0.00, lambda_value, lambda_value, lambda_value)
    alchemical_state.annihilateElectrostatics = annihilateElectrostatics
    alchemical_state.annihilateLennardJones = annihilateLennardJones

    #platform_name = 'Reference' # DEBUG
    platform = openmm.Platform.getPlatformByName(platform_name)
    
    # Create the perturbed system.
    print "Creating alchemically-modified state..."
    initial_time = time.time()
    alchemical_system = factory.createPerturbedSystem(alchemical_state)    
    final_time = time.time()
    elapsed_time = final_time - initial_time
    # Compare energies.
    timestep = 1.0 * units.femtosecond
    print "Computing reference energies..."
    reference_integrator = openmm.VerletIntegrator(timestep)
    reference_context = openmm.Context(reference_system, reference_integrator, platform)
    reference_context.setPositions(coordinates)
    reference_state = reference_context.getState(getEnergy=True)
    reference_potential = reference_state.getPotentialEnergy()    
    print "Computing alchemical energies..."
    alchemical_integrator = openmm.VerletIntegrator(timestep)
    alchemical_context = openmm.Context(alchemical_system, alchemical_integrator, platform)
    alchemical_context.setPositions(coordinates)
    alchemical_state = alchemical_context.getState(getEnergy=True)
    alchemical_potential = alchemical_state.getPotentialEnergy()
    delta = alchemical_potential - reference_potential 
    print "reference system       : %24.8f kcal/mol" % (reference_potential / units.kilocalories_per_mole)
    print "alchemically modified  : %24.8f kcal/mol" % (alchemical_potential / units.kilocalories_per_mole)
    print "ERROR                  : %24.8f kcal/mol" % ((alchemical_potential - reference_potential) / units.kilocalories_per_mole)
    print "elapsed alchemical time  %.3f s" % elapsed_time

    return delta

def test_overlap():
    """
    BUGS TO REPORT:
    * Even if epsilon = 0, energy of two overlapping atoms is 'nan'.
    * Periodicity in 'nan' if dr = 0.1 even in nonperiodic system
    """

    # Create a reference system.    
    import testsystems

    print "Creating Lennard-Jones cluster system..."
    #[reference_system, coordinates] = testsystems.LennardJonesFluid()
    #receptor_atoms = [0]
    #ligand_atoms = [1]

    [reference_system, coordinates] = testsystems.LysozymeImplicit()
    receptor_atoms = range(0,2603) # T4 lysozyme L99A
    ligand_atoms = range(2603,2621) # p-xylene

    import simtk.unit as units
    unit = coordinates.unit
    coordinates = units.Quantity(numpy.array(coordinates / unit), unit)

    factory = AbsoluteAlchemicalFactory(reference_system, ligand_atoms=ligand_atoms)
    alchemical_state = AlchemicalState(0.00, 0.00, 0.00, 1.0)

    # Create the perturbed system.
    print "Creating alchemically-modified state..."
    alchemical_system = factory.createPerturbedSystem(alchemical_state)    
    # Compare energies.
    import simtk.unit as units
    import simtk.openmm as openmm
    timestep = 1.0 * units.femtosecond
    print "Computing reference energies..."
    integrator = openmm.VerletIntegrator(timestep)
    context = openmm.Context(reference_system, integrator)
    context.setPositions(coordinates)
    state = context.getState(getEnergy=True)
    reference_potential = state.getPotentialEnergy()    
    del state, context, integrator
    print reference_potential
    print "Computing alchemical energies..."
    integrator = openmm.VerletIntegrator(timestep)
    context = openmm.Context(alchemical_system, integrator)
    dr = 0.1 * units.angstroms # TODO: Why does 0.1 cause periodic 'nan's?
    a = receptor_atoms[-1]
    b = ligand_atoms[-1]
    delta = coordinates[a,:] - coordinates[b,:]
    for k in range(3):
        coordinates[ligand_atoms,k] += delta[k]
    for i in range(30):
        r = dr * i
        coordinates[ligand_atoms,0] += dr
          
        context.setPositions(coordinates)
        state = context.getState(getEnergy=True)
        alchemical_potential = state.getPotentialEnergy()    
        print "%8.3f A : %f " % (r / units.angstroms, alchemical_potential / units.kilocalories_per_mole)
    del state, context, integrator

    return

def test_intermediates():
    # Run tests on individual systems.
    import testsystems

    print "Creating Lennard-Jones fluid system without dispersion correction..."
    [reference_system, coordinates] = testsystems.LennardJonesFluid(dispersion_correction=False)
    ligand_atoms = range(0,1) # first atom
    receptor_atoms = range(2,3) # second atom
    testAlchemicalFactory(reference_system, coordinates, receptor_atoms, ligand_atoms)
    print ""

    print "Creating Lennard-Jones fluid system with dispersion correction..."
    [reference_system, coordinates] = testsystems.LennardJonesFluid(dispersion_correction=True)
    ligand_atoms = range(0,1) # first atom
    receptor_atoms = range(2,3) # second atom
    testAlchemicalFactory(reference_system, coordinates, receptor_atoms, ligand_atoms)
    print ""

    print "Creating T4 lysozyme system..."
    [reference_system, coordinates] = testsystems.LysozymeImplicit()
    receptor_atoms = range(0,2603) # T4 lysozyme L99A
    ligand_atoms = range(2603,2621) # p-xylene
    testAlchemicalFactory(reference_system, coordinates, receptor_atoms, ligand_atoms)    
    print ""

    print "Creating Lennard-Jones cluster..."
    [reference_system, coordinates] = testsystems.LennardJonesCluster()
    ligand_atoms = range(0,1) # first atom
    receptor_atoms = range(1,2) # second atom
    testAlchemicalFactory(reference_system, coordinates, receptor_atoms, ligand_atoms)
    print ""

    print "Creating alanine dipeptide implicit system..."
    [reference_system, coordinates] = testsystems.AlanineDipeptideImplicit()
    ligand_atoms = range(0,4) # methyl group
    receptor_atoms = range(4,22) # rest of system
    testAlchemicalFactory(reference_system, coordinates, receptor_atoms, ligand_atoms)
    print ""

    print "Creating alanine dipeptide explicit system..."
    [reference_system, coordinates] = testsystems.AlanineDipeptideExplicit()
    ligand_atoms = range(0,22) # alanine residue
    receptor_atoms = range(22,25) # one water
    testAlchemicalFactory(reference_system, coordinates, receptor_atoms, ligand_atoms)
    print ""

    print "Creating alanine dipeptide explicit system without dispersion correction..."
    #forces = { reference_system.getForce(index).__class__.__name__ : reference_system.getForce(index)) for index in range(reference_system.getNumForces()) } # requires Python 2.7 features
    forces = dict( (reference_system.getForce(index).__class__.__name__,  reference_system.getForce(index)) for index in range(reference_system.getNumForces()) ) # python 2.6 compatible
    forces['NonbondedForce'].setUseDispersionCorrection(False) # turn off dispersion correction
    ligand_atoms = range(0,22) # alanine residue
    receptor_atoms = range(22,25) # one water
    testAlchemicalFactory(reference_system, coordinates, receptor_atoms, ligand_atoms)
    print ""

#=============================================================================================
# MAIN
#=============================================================================================

if __name__ == "__main__":
    # Run overlap tests.
    test_overlap()
    
    # Test energy accuracy of intermediates near lambda = 1.
    test_intermediates()
