# Noble Gas BEB Ionization Cross Sections

This directory contains implementations for calculating Binary Encounter Bethe (BEB) ionization cross sections for noble gases (He, Ne, Ar, Kr, Xe).

## Directory Structure

```
noble_gas_beb/
├── common/              # Shared utilities
│   ├── beb_calculator.py   # Core BEB calculation class
│   └── plot_utils.py       # Plotting utilities
├── He/                  # Helium calculations
├── Ne/                  # Neon calculations
├── Ar/                  # Argon calculations
├── Kr/                  # Krypton calculations
├── Xe/                  # Xenon calculations
├── results/             # Output directory for all results
└── run_all_noble_gases.py  # Main workflow runner
```

## Usage

### Run all calculations:
```bash
python run_all_noble_gases.py
```

### Run individual element:
```bash
python He/he_beb_calc.py
python Ne/ne_beb_calc.py
# etc.
```

## Requirements

- Python 3.7+
- numpy
- matplotlib

## Gaussian Calculations

Each element directory contains:
- `*_beb_input.gjf` - Gaussian input file with IOp(6/80=1) for kinetic energies
- `*_beb_calc.py` - BEB calculation script

If Gaussian output files are not available, the scripts use literature values for orbital binding and kinetic energies.

## Output

Results are saved in the `results/` directory:
- `results/<element>/<element>_beb_results.dat` - Numerical cross section data
- `results/<element>/<element>_beb_cross_section.png/pdf/svg` - Cross section plots
- `results/<element>/<element>_orbital_contributions.png` - Orbital contribution plots
- `results/noble_gas_comparison.png` - Comparison plot of all elements
- `results/noble_gas_summary.png` - Summary grid plot

## References

Experimental data sources:
- He: Shah et al., J. Phys. B 20, 3501 (1987)
- Ne: Krishnakumar & Srivastava, J. Phys. B 21, 1055 (1988)
- Ar, Kr, Xe: Rejoub et al., Phys. Rev. A 65, 042713 (2002)