# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is the PHOTOX-knowhow repository, a comprehensive toolkit for computational chemistry calculations focusing on:
- UKRmol+ R-matrix calculations for electron scattering and photoionization
- Quantum chemistry interfaces (Gaussian, Molpro, ORCA, etc.)
- Binary Encounter Bethe (BEB) model for ionization cross sections

## Key Directories

- `/scripts/` - Modern Python implementation for UKRmol+ setup
- `/LAUNCH/` - Cluster job submission scripts for various programs
- `/INPUTS/` - Example input files for quantum chemistry programs
- `/projects/` - Actual calculation projects (atoms and molecules)
- `/MISC/BEB/` - BEB ionization model implementation

## Development Commands

### Python Scripts Setup (in `/scripts/`)
```bash
# Install dependencies
cd scripts
uv sync  # or: uv pip install -e .

# Format code
black .

# No tests currently exist
```

### Running UKRmol+ Calculations
```bash
# Generate configuration files using Python
python -c "from ukrmol_setup import create_atom_setup; create_atom_setup('Ne', 'cc-pVDZ').generate_all_files('output_dir')"

# Run calculation with Perl workflow
perl main.pl config.pl model.pl dirs.pl geometry.pl
```

## Architecture

### Python-Perl Integration
The modern Python scripts generate Perl configuration files for the established UKRmol+ workflow:

1. **Python Setup Layer** (`/scripts/`):
   - `ukrmol_setup.py` - High-level integrated setup class
   - `model.py` - Calculation model definitions (CAS, SE, CHF)
   - `config.py` - Runtime configuration
   - `geometry.py` - Molecular geometry using ASE
   - `dirs.py` - Directory structure management

2. **Perl Execution Layer**:
   - `main.pl` - Main workflow orchestrator
   - Generated `.pl` config files from Python scripts
   - Integration with quantum chemistry programs

### Calculation Workflow
1. Define calculation parameters using Python classes
2. Generate Perl configuration files
3. Execute main.pl with configs
4. Results organized in structured output directories

### Key Classes and Patterns

**UKRmolSetup** (ukrmol_setup.py:10):
- Central class integrating all configuration components
- Methods: `generate_all_files()`, `from_config_files()`
- Factory functions: `create_atom_setup()`, `create_h2_setup()`, `create_co_setup()`

**Config Classes**:
- `Model` (model.py:7) - Defines calculation types and parameters
- `Config` (config.py:7) - Runtime options and code selection
- `Dirs` (dirs.py:7) - Directory paths and structure
- `Geometry` (geometry.py:12) - Molecular structures via ASE

### Environment Requirements
- Python ≥ 3.11 with ASE and Black
- Perl with required modules
- Access to quantum chemistry programs (Molpro, Psi4, etc.)
- PHOTOX cluster environment (scripts assume specific paths and queuing system)

## Important Notes

- The system is designed for the PHOTOX computational cluster
- Launch scripts in `/LAUNCH/` handle environment setup and job submission
- Results include scattering cross sections, photoionization data, and resonance parameters
- The BEB model (`/MISC/BEB/beb.py`) provides electron impact ionization calculations