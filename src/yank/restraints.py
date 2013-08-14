#!/usr/local/bin/env python
#
#=============================================================================================
# FILE DOCSTRING
#=============================================================================================

"""

Automated selection and imposition of receptor-ligand restraints for absolute alchemical binding
free energy calculations.

DESCRIPTION

The restraints are chosen in such a way as to be able to provide a standard state binding free
energy.

@author John D. Chodera <jchodera@gmail.com>

Portions of this code copyright (c) 2009 University of California, Berkeley, Vertex
Pharmaceuticals, Stanford University, and the Authors.

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

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import copy
import math
import numpy

import scipy.integrate

import simtk.unit as units
import simtk.openmm as openmm

#=============================================================================================
# MODULE CONSTANTS
#=============================================================================================

kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA # Boltzmann constant

#=============================================================================================
# Base class for receptor-ligand restraints.
#=============================================================================================

class ReceptorLigandRestraint(object):
    """
    Impose a single restraint between ligand and protein to prevent ligand from drifting too far
    from protein in implicit solvent calculations.

    This restraint strength is controlled by a global context parameter called 'restraint_lambda'.

    NOTE

    These restraints should not be used with explicit solvent calculations, since the CustomBondForce
    does not respect periodic boundary conditions, and the analytical correction term does not include
    truncation due to a finite simulation box.

    EXAMPLE
        
    >>> # Create a test system.
    >>> import testsystems
    >>> [system, coordinates] = testsystems.LysozymeImplicit()
    >>> # Identify receptor and ligand atoms.
    >>> receptor_atoms = range(0,2603)
    >>> ligand_atoms = range(2603,2621)
    >>> # Construct a reference thermodynamic state.
    >>> from thermodynamics import ThermodynamicState
    >>> temperature = 298.0 * units.kelvin
    >>> state = ThermodynamicState(temperature=temperature)
    >>> # Create restraints.
    >>> restraints = ReceptorLigandRestraint(state, system, coordinates, receptor_atoms, ligand_atoms)
    >>> # Get standard state correction.
    >>> print restraints.getStandardStateCorrection()
    0.365059033167

    NOTES

    To create a subclass that uses a different restraint energy function, follow these steps:

    * Redefine class variable 'energy_function' with the energy function of choice.
    * Redefine class variable 'bond_parameter_names' to list the parameters in 'energy_function'
      that must be computed for each system.
    * Redefine the _determineBondParameters() member function to compute these parameters and
      return them in a list in the same order as 'bond_parameter_names'.

    """
    
    energy_function = 'restraint_lambda * (K/2) * r^2' # harmonic restraint
    bond_parameter_names = ['K'] # list of bond parameters that appear in energy function above

    def __init__(self, state, system, coordinates, receptor_atoms, ligand_atoms, verbose=False):
        """
        Initialize a receptor-ligand restraint class.

        ARGUMENTS
        
        state (thermodynamics.ThermodynamicState) - the thermodynamic state specifying tempearture, pressure, etc. to which restraints are to be added
        coordinates (simtk.unit.Quantity of natoms x 3 with units compatible with nanometers) - reference coordinates to use for imposing restraints
        receptor_atoms (list of int) - a complete list of receptor atoms
        ligand_atoms (list of int) - a complete list of ligand atoms

        OPTIONAL ARGUMENTS

        verbose (bool) - if True, will print verbose output

        """

        self.state = state
        self.system = system
        self.coordinates = units.Quantity(numpy.array(coordinates / coordinates.unit), coordinates.unit)
        self.receptor_atoms = list(receptor_atoms)
        self.ligand_atoms = list(ligand_atoms)

        self.verbose = verbose

        self.temperature = state.temperature
        self.kT = kB * self.temperature # thermal energy
        self.beta = 1.0 / self.kT # inverse temperature

        # Determine atoms closet to centroids on ligand and receptor.
        self.restrained_receptor_atom = self._closestAtomToCentroid(self.coordinates, self.receptor_atoms)
        self.restrained_ligand_atom = self._closestAtomToCentroid(self.coordinates, self.ligand_atoms) 
        
        if self.verbose: 
            print "restrained receptor atom: %d" % self.restrained_receptor_atom
            print "restrained ligand atom: %d" % self.restrained_ligand_atom

        # Determine parameters
        self.bond_parameters = self._determineBondParameters()

        # Determine standard state correction.
        self.standard_state_correction = self._computeStandardStateCorrection()

        return

    def _determineBondParameters(self):
        """
        Determine bond parameters for CustomBondForce between protein and ligand.

        RETURNS

        parameters (list) - list of parameters for CustomBondForce

        NOTE

        The spring constant is selected to give 1 kT at one standard deviation of
        receptor atoms about the receptor restrained atom.
        
        """

        unit = self.coordinates.unit
        
        # Get dimensionless receptor coordinates.
        x = self.coordinates[self.receptor_atoms,:] / unit
        
        # Get dimensionless restrained atom coordinate.
        xref = self.coordinates[self.restrained_receptor_atom,:] / unit # (3,) array
        xref = numpy.reshape(xref, (1,3)) # (1,3) array
        
        # Compute distances from restrained atom.
        natoms = x.shape[0]
        distances = numpy.sqrt(((x - numpy.tile(xref, (natoms, 1)))**2).sum(1)) # distances[i] is the distance from the centroid to particle i

        # Compute std dev of distances from restrained atom.
        sigma = distances.std() * unit 

        # Compute corresponding spring constant.
        K = self.kT / sigma**2

        # Assemble parameters.
        bond_parameters = [K]

        return bond_parameters

    def _createRestraintForce(self, particle1, particle2, mm=None):
        """
        Create a new copy of the receptor-ligand restraint force.

        RETURNS

        force (simtk.openmm.CustomBondForce) - a restraint force object

        """

        if mm is None: mm = openmm
        
        force = openmm.CustomBondForce(self.energy_function)
        force.addGlobalParameter('restraint_lambda', 1.0)
        for parameter in self.bond_parameter_names:
            force.addPerBondParameter(parameter)
        force.addBond(particle1, particle2, self.bond_parameters)
        
        return force

    def _computeStandardStateCorrection(self):
        """
        Compute the standard state correction for the arbitrary restraint energy function.
        
        RETURN VALUES

        DeltaG (float) - computed standard-state correction in dimensionless units (kT);

        NOTE

        Equivalent to the free energy of releasing restraints and confining into a box of standard state size.
                
        """

        verbose = False
        r_min = 0 * units.nanometers
        r_max = 100 * units.nanometers # TODO: Use maximum distance between atoms?

        # Create a System object containing two particles connected by the reference force
        system = openmm.System()
        system.addParticle(1.0 * units.amu)
        system.addParticle(1.0 * units.amu)
        force = self._createRestraintForce(0, 1)
        system.addForce(force)

        # Create a Reference context to evaluate energies on the CPU.
        integrator = openmm.VerletIntegrator(1.0 * units.femtoseconds)
        platform = openmm.Platform.getPlatformByName('Reference')
        context = openmm.Context(system, integrator, platform)
        
        # Set default positions.
        positions = units.Quantity(numpy.zeros([2,3]), units.nanometers)
        context.setPositions(positions)

        # Create a function to compute integrand as a function of interparticle separation.
        beta = self.beta
        def integrand(r):
            """
            ARGUMENTS
            
            r (float) - interparticle separation in nanometers
            
            RETURNS

            dI (float) - contribution to integrand (in nm^2)
            """
            positions[1,0] = r * units.nanometers
            context.setPositions(positions)
            state = context.getState(getEnergy=True)
            potential = state.getPotentialEnergy()
            dI = 4.0 * math.pi * r**2 * math.exp(-beta * potential)
            return dI

        (shell_volume, shell_volume_error) = scipy.integrate.quad(lambda r : integrand(r), r_min / units.nanometers, r_max / units.nanometers) * units.nanometers**3 # integrate shell volume
        if verbose: print "shell_volume = %f nm^3" % (shell_volume / units.nanometers**3)
        
        # Compute standard-state volume for a single molecule in a box of size (1 L) / (avogadros number)
        liter = 1000.0 * units.centimeters**3 # one liter        
        box_volume = liter / (units.AVOGADRO_CONSTANT_NA*units.mole) # standard state volume
        if verbose: print "box_volume = %f nm^3" % (box_volume / units.nanometers**3)
        
        # Compute standard state correction for releasing shell restraints into standard-state box (in units of kT).
        DeltaG = - math.log(box_volume / shell_volume)
        if verbose: print "Standard state correction: %.3f kT" % DeltaG
        
        # Return standard state correction (in kT).
        return DeltaG

    def getRestraintForce(self, mm=None):
        """
        Returns a new Force object that imposes the receptor-ligand restraint.
        
        OPTIONAL ARGUMENTS

        mm (simtk.openmm interface) - OpenMM implementation to use

        """

        return self._createRestraintForce(self.restrained_receptor_atom, self.restrained_ligand_atom, mm=mm)

    def getRestrainedSystemCopy(self):
        """
        Returns a copy of the restrained system.
        
        RETURNS

        system (simtk.openmm.System) - a copy of the restrained system
        
        """
        system = copy.deepcopy(self.system)
        force = self.getRestraintForce()
        system.addForce(force)

        return system
    
    def getStandardStateCorrection(self):
        """
        Return the standard state correction.

        RETURNS

        correction (float) - the standard-state correction, in kT

        """
        return self.standard_state_correction

    def _closestAtomToCentroid(self, coordinates, indices=None, masses=None):
        """
        Identify the closest atom to the centroid of the given coordinate set.

        ARGUMENTS
        
        coordinates (units.Quantity of natoms x 3 with units compatible with nanometers) - coordinates of object to identify atom closes to centroid

        OPTIONAL ARGUMENTS

        masses (units.Quantity of natoms with units compatible with amu) - masses of particles used to weight distance calculation, if not None (default: None)

        """

        if indices is not None:
            coordinates = coordinates[indices,:]

        # Get dimensionless coordinates.
        x_unit = coordinates.unit
        x = coordinates / x_unit
        
        # Determine number of atoms.
        natoms = x.shape[0]

        # Compute (natoms,1) array of normalized weights.
        w = numpy.ones([natoms,1])
        if masses:            
            w = masses / masses.unit # (natoms,) array
            w = numpy.reshape(w, (natoms,1)) # (natoms,1) array            
        w /= w.sum()

        # Compute centroid (still in dimensionless units).
        centroid = (numpy.tile(w, (1,3)) * x).sum(0) # (3,) array
        
        # Compute distances from centroid.
        distances = numpy.sqrt(((x - numpy.tile(centroid, (natoms, 1)))**2).sum(1)) # distances[i] is the distance from the centroid to particle i
        
        # Determine closest atom.
        closest_atom = int(numpy.argmin(distances))
        
        if indices is not None:
            closest_atom = indices[closest_atom]

        return closest_atom

#=============================================================================================
# Harmonic protein-ligand restraint.
#=============================================================================================

class FlatBottomReceptorLigandRestraint(ReceptorLigandRestraint):
    """
    An alternative choice to receptor-ligand restraints that uses a flat potential inside most of the protein volume
    with harmonic restraining walls outside of this.

    EXAMPLE
        
    >>> # Create a test system.
    >>> import testsystems
    >>> [system, coordinates] = testsystems.LysozymeImplicit()
    >>> # Identify receptor and ligand atoms.
    >>> receptor_atoms = range(0,2603)
    >>> ligand_atoms = range(2603,2621)
    >>> # Construct a reference thermodynamic state.
    >>> from thermodynamics import ThermodynamicState
    >>> temperature = 298.0 * units.kelvin
    >>> state = ThermodynamicState(temperature=temperature)
    >>> # Create restraints.
    >>> restraints = FlatBottomReceptorLigandRestraint(state, system, coordinates, receptor_atoms, ligand_atoms)
    >>> # Get standard state correction.
    >>> print restraints.getStandardStateCorrection()
    4.74636082004
    
    """

    energy_function = 'restraint_lambda * step(r-r0) * (K/2)*(r-r0)^2' # flat-bottom restraint
    bond_parameter_names = ['K', 'r0'] # list of bond parameters that appear in energy function above

    def _determineBondParameters(self):
        """
        Determine bond parameters for CustomBondForce between protein and ligand.

        RETURNS

        parameters (list) - list of parameters for CustomBondForce

        NOTE

        r0, the distance at which the harmonic restraint is imposed, is set at half the maximum distance between receptor atoms plus 5 A
        K, the spring constant, is set to 5.92 kcal/mol/A**2
        
        """

        unit = self.coordinates.unit
        
        # Get dimensionless receptor coordinates.
        x = self.coordinates[self.receptor_atoms,:] / unit        

        # Compute maximum interatomic distance in protein.
        # (Working in non-unit-bearing floats for speed.)
        max_distance = 0.0
        natoms = x.shape[0]
        for i in range(natoms):
            # Extract current atom position.
            xref = numpy.reshape(x[i,:], (1,3)) # (1,3) array            
            # Compute distances from all atoms to current atom.
            distances = numpy.sqrt(((x - numpy.tile(xref, (natoms, 1)))**2).sum(1)) # distances[i] is the distance from the centroid to particle i
            # Keep maximum distance.
            max_distance = max(max_distance, distances.max())
        # Convert back to unit-bearing quantity.
        max_distance *= unit     
           
        # Calculate r0, which is half of the longest distance in protein plus 5 A.
        r0 = max_distance/2.0 + 5.0 * units.angstroms

        # Set spring constant/
        K = (2.0 * 0.0083144621 * 5.0 * 298.0 * 100) * units.kilojoules_per_mole/units.nanometers**2

        # Assemble parameter vector.
        bond_parameters = [K, r0]

        return bond_parameters

#=============================================================================================
# DOCTEST HARNESS
#=============================================================================================

if __name__ == '__main__':
    import doctest
    doctest.testmod()
