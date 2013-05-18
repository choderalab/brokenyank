YANK design notes
=================

Manifest
--------
* `yank.py` - YANK main module and driver script
* `restraints/` - tools for imposing restraints
* `pyopenmm/` - Pythonic extensions to OpenMM
* `alchemy/` - tools for generating alchemically-modified systems
* `sampling/` - tools for sampling statistical ensembles
* `analysis/` - tools for analyzing alchemical free energy calculations

Files and classes
-----------------

`yank.py` - driver for alchemical binding free energy calculations

* `Yank`
  - Driver class for setting up and running alchemical free energy calculations from other programs

* `ModifiedHamiltonianExchange`
  - Modified HamiltonianExchange class that adds Monte Carlo moves to Langevin dynamics to sample replicas
  - (Temporary.  Will be replaced by general MCMC facility.)

* `__main__`
  - Command-line driver for Yank

--------------------------------------------------------------------------------

`restraints.py` - tool for selecting and imposing receptor-ligand restraints

* `ReceptorLigandRestraints`
  - Base class for automated selection and imposition of receptor-ligand restraints for standard-state binding free energy calculations

* `HarmonicReceptorLigandRestraints`
  - Imposition of a single harmonic restraint between receptor and ligand.

* `FlatBottomReceptorLigandRestraints`
  - Imposition of a flat-bottomed harmonic restraint between receptor and ligand.

--------------------------------------------------------------------------------

`thermodynamics.py` - useful classes for statistical mechanics

* `ThermodynamicState`
  - General thermodynamic state (e.g. Hamiltonian, temperature, pressure, pH) description.

--------------------------------------------------------------------------------

`analysis.py` - tools for analysis of free energy calculations

* `Analysis`
  - Abstract base class for analysis of simulation data

* `ReplicaExchangeAnalysis`
  Analysis of replica-exchange simulation data

--------------------------------------------------------------------------------

`alchemy.py` - alchemical modification of System objects

* `AlchemicalFactory`
  - Abstract base class for an enumerative factory for creating System objects for use in replica exchange alchemical free energy calculations

* `AbsoluteAlchemicalFactory`
  - Factory for creating alchemically-modified systems in which one component (e.g. a ligand or solute) is annihilated/decoupled from the environment

--------------------------------------------------------------------------------

`repex.py` - replica-exchange simulation algorithms

* `ReplicaExchange`
  - Replica-exchange simulation.

* `ParallelTempering`
  - Parallel tempering (replica-exchange among temperatures)

* `HamiltonianExchange`
  Hamiltonian exchange (replica-exchange among System objects)

--------------------------------------------------------------------------------

`expanded.py` - expanded-ensemble simulation algorithms

* `ExpandedEnsemble`
  - General expanded-ensemble simulation.

* `SimulatedTempering`
  - Simulated tempering simulation (expanded ensembles among temperatures)

* `SimulatedScaling` 
  - Simulated scaling simulation (expanded ensembles among Hamiltonians)

--------------------------------------------------------------------------------

Planned and future expansions
-----------------------------

`combinatorial_factories.py` - combinatorial System factories - INCOMPLETE

* `SubstituentFactory` (*)
  - Factory for producing small molecules with different substituents (R-groups)

* `MutationFactory` (*)
  - Factory for producing receptors with amino acid mutations

--------------------------------------------------------------------------------

`mcmc.py` - simulation tools based on Markov chain Monte Carlo (MCMC) - INCOMPLETE

* `MCMCSampler`
  - Markov chain Monte Carlo simulation driver

* `MarkovChainMonteCarloMove`
  - Abstract base class for MCMC moves

* `LangevinDynamicsMove`
  - Langevin dynamics move (without Metropolization); only approximate

* `HybridMonteCarloMove` (*)
  - Hybrid Monte Carlo move.

* `GeneralizedHybridMonteCarloMove` (*)
  - Generalized hybrid Monte Carlo (GHMC) move

* `MonteCarloVolumeMove` (*)
  - Monte Carlo system volume move for constant-pressure simulation

* `SidechainRotamerMove` (*)
  - Sidechain rotamer MCMC move using Straub method.

* `ProtonationStateMove` (*)
  - Protonation state MCMC move (for constant pH simulations)

* `TautomericStateMove` (*)
  - Ligand tautomeric state MCMC move 

* `CounterionSamplingMove` (*)
  - Counterion insertion/deletion MCMC move

--------------------------------------------------------------------------------

simulation.py - general simulation tools - INCOMPLETE

* `SamplerState`
  - Abstract base class for a 'state' that describes the probability density of a tuple of (configuration, volume, protonation state, number of ions, etc.)

--------------------------------------------------------------------------------

enumerative_factories.py - enumerative System factories - INCOMPLETE

* `AlchemicalFactory`
  - Abstract base class for an enumerative factory for creating System objects for use in replica exchange alchemical free energy calculations

* `AbsoluteAlchemicalFactory`
  - Factory for creating alchemically-modified systems in which one component (e.g. a ligand or solute) is annihilated/decoupled from the environment

* `RelativeAlchemicalFactory`
  - Factory for creating alchemically-modified systems in which one species is alchemically transmuated into another species using a "dual-topology" scheme

* `OneStepPerturbationFactory`
  - Factory for creating a family of chemical species that share one or more alchemical "soft-core" intermediates

--------------------------------------------------------------------------------

(*) indicates implementation will be postponed until this functionality is needed

