# UKRmol - UK R-matrix Method for Molecular Scattering Calculations

## Overview

UKRmol is an implementation of the UK R-matrix method for electron-molecule and positron-molecule scattering calculations. This project provides a comprehensive toolset for calculating electron scattering cross sections, photoionization cross sections, resonance positions and widths, and target state properties of molecules.

## Key Features

- **Electron Scattering Calculations**: Elastic and inelastic scattering cross sections
- **Photoionization**: Molecular photoionization cross section calculations
- **Resonance Analysis**: Identification of resonance positions and widths
- **Multiple Calculation Models**: SE, SEP, CAS, CHF, MAS, and more
- **Parallel Computing**: Support for parallel execution over geometries and symmetries

## Directory Structure

```
UKRmol/
├── basis.sets/          # Quantum chemistry basis sets
│   ├── sto-3g/         # Minimal basis
│   ├── 6-31g/          # Pople basis sets
│   ├── cc-pvtz/        # Correlation consistent basis
│   └── continuum/      # Continuum basis sets
├── input.templates/     # Calculation template files
│   ├── molpro.inp      # Molpro input template
│   ├── psi4.inp        # Psi4 input template
│   ├── scatci_integrals.inp  # Integral transformation
│   ├── congen.inp      # Configuration generation
│   └── scattering.*.inp # Scattering calculations
├── projects/           # Molecule-specific calculation projects
│   ├── Ne/            # Neon
│   ├── He/            # Helium
│   ├── Ar/            # Argon
│   ├── Kr/            # Krypton
│   ├── Xe/            # Xenon
│   └── H2O/           # Water molecule
├── resources/          # Libraries and utilities
│   ├── lib/           # Perl modules
│   └── scripts/       # Automation scripts
└── photoxrepo/        # Additional tools and BEB model implementation
```

## Prerequisites

### Required Software

- **Quantum Chemistry Packages** (at least one of):
  - Molpro
  - Psi4
  - Gaussian
  - ORCA
  - Molcas
- **R-matrix Codes**:
  - SCATCI (scattering calculations)
  - CONGEN (configuration generation)
  - RSOLVE (outer region)
- **Programming Languages**:
  - Perl 5.x
  - Python 3.x (for plotting)
- **Parallel Computing** (optional):
  - MPI implementation (OpenMPI, MPICH, etc.)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/UKRmol.git
cd UKRmol
```

2. Set environment variables:
```bash
export UKRMOL_HOME=/path/to/UKRmol
export PATH=$UKRMOL_HOME/resources/scripts:$PATH
```

3. Install required Perl modules:
```bash
cpan install Config::General
cpan install Math::Complex
```

## Usage

### Basic Workflow

1. **Create a project directory**:
```bash
cd projects/
mkdir my_molecule
cd my_molecule
```

2. **Create configuration files**:
   - `geometry.pl`: Define molecular structure
   - `model.pl`: Calculation model parameters
   - `config.pl`: Execution settings
   - `main.pl`: Main execution script

3. **Run the calculation**:
```bash
perl main.pl
```

### Configuration Example (H2O molecule)

`geometry.pl`:
```perl
$geometry = {
    atoms => [
        { element => 'O',  x =>  0.000,  y => 0.000,  z => 0.000 },
        { element => 'H',  x =>  0.757,  y => 0.586,  z => 0.000 },
        { element => 'H',  x => -0.757,  y => 0.586,  z => 0.000 }
    ],
    charge => 0,
    multiplicity => 1
};
```

`model.pl`:
```perl
$model = {
    type => 'SEP',              # Static Exchange + Polarization
    basis => 'cc-pvtz',         # Basis set
    active_space => [8, 4],     # CAS(8,4)
    states => 5,                # Number of states to calculate
    energy_range => [0, 20],    # Energy range (eV)
};
```

## Calculation Models

- **SE (Static Exchange)**: Static exchange approximation
- **SEP (Static Exchange + Polarization)**: Including polarization effects
- **CAS (Complete Active Space)**: Complete active space method
- **CHF (Coupled Hartree-Fock)**: Coupled Hartree-Fock method
- **MAS (Multiple Active Spaces)**: Multiple active spaces approach

## Output Files

After calculation, the following files are generated:

- `cross_sections.dat`: Scattering cross section data
- `eigenphases.dat`: Eigenphase shifts
- `resonances.dat`: Resonance parameters
- `photoionization.dat`: Photoionization cross sections

## Plotting

Visualize results:
```bash
python plot.py --input cross_sections.dat --output plot.png
```

## Troubleshooting

### Common Issues

1. **Out of memory errors**:
   - Increase memory allocation in `config.pl`
   - Use smaller basis sets

2. **Convergence problems**:
   - Improve initial guess
   - Relax convergence criteria

3. **Parallel execution issues**:
   - Verify MPI is correctly installed
   - Check inter-node communication

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

[Specify license type]

## Citation

If you use this software, please cite:
```
[Add appropriate citation information]
```

## Contact

For questions or issues, please contact [email address].