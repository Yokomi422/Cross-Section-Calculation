#!/usr/bin/env python3
"""
Calculate BEB ionization cross sections for He atom
"""
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import BEBCalculator, BEBPlotter
from pathlib import Path
import numpy as np
import re

def parse_he_gaussian_improved(calc, filename):
    """
    Parse Gaussian output with IOp(6/80=1) for orbital energies and kinetic energies
    """
    with open(filename, 'r') as f:
        content = f.read()
    
    # 新しい形式の軌道エネルギーと運動エネルギーを探す
    pattern = r'Orbital energies and kinetic energies \(alpha\):\s*\n\s*1\s+2\s*\n\s*1\s+\(A1G\)--O\s+([-\d.]+)\s+([\d.]+)'
    match = re.search(pattern, content)
    
    if match:
        E_orb = float(match.group(1))  # Hartree (negative value)
        Ekin_orb = float(match.group(2))  # Hartree
        
        # BEBでは正の束縛エネルギーが必要
        calc.E_orb = [-E_orb]  # 正の束縛エネルギー
        calc.Ekin_orb = [Ekin_orb]  # 運動エネルギー
        calc.occ = [2]  # ヘリウムの1s軌道の占有数
        
        print(f"\n# Improved parsing results for He:")
        print(f"# Orbital energy: {E_orb:.6f} Hartree = {E_orb * 27.2114:.2f} eV")
        print(f"# Binding energy: {-E_orb:.6f} Hartree = {-E_orb * 27.2114:.2f} eV")
        print(f"# Orbital kinetic energy: {Ekin_orb:.6f} Hartree = {Ekin_orb * 27.2114:.2f} eV")
        
        return True
    
    return False

def main():
    # Set up paths
    current_dir = Path(__file__).parent
    output_dir = current_dir.parent / 'results' / 'He'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize calculator
    calc = BEBCalculator('He')
    
    # Parse Gaussian output
    gaussian_file = current_dir / 'he_beb.log'
    
    if not gaussian_file.exists():
        print(f"Error: {gaussian_file} not found!")
        return
    
    # Try different parsing approaches
    try:
        # First try the standard parser
        calc.parse_gaussian_output(gaussian_file)
        print("# Standard parsing successful")
    except ValueError:
        # If that fails, try our improved parser
        print("# Standard parsing failed. Trying improved parser...")
        if not parse_he_gaussian_improved(calc, gaussian_file):
            # If that also fails, use manual values
            print("# Improved parsing failed. Using manual values...")
            calc.E_orb = [0.661440]  # From the output file
            calc.Ekin_orb = [1.432451]  # From the output file
            calc.occ = [2]
    
    # Print summary
    calc.print_summary()
    
    # Calculate cross sections
    E_incident, cross_sections, orbital_contributions = calc.calculate_cross_sections()
    
    # Save results
    output_file = output_dir / 'he_beb_results.dat'
    calc.save_results(E_incident, cross_sections, output_file)
    print(f"\n# Results saved to: {output_file}")
    
    # Create plots if BEBPlotter is available
    try:
        plotter = BEBPlotter('He')
        plotter.plot_cross_section(E_incident, cross_sections, output_dir)
        plotter.plot_orbital_contributions(E_incident, orbital_contributions,
                                         orbital_names=['1s'], output_dir=output_dir)
        print(f"# Plots saved to: {output_dir}")
    except Exception as e:
        print(f"# Warning: Could not create plots: {e}")
    
    # 実験値との比較
    from common.beb_calculator import get_experimental_data
    exp_data = get_experimental_data('He')
    if exp_data:
        print(f"\n# Experimental data reference: {exp_data['reference']}")
        print("# Comparison at selected energies:")
        for E_exp, sigma_exp in zip(exp_data['energy'][:5], exp_data['cross_section'][:5]):
            idx = np.argmin(np.abs(E_incident - E_exp))
            if abs(E_incident[idx] - E_exp) < 5:
                print(f"#   E = {E_exp:4.0f} eV: Calc = {cross_sections[idx]:.3f}, Exp = {sigma_exp:.3f} Å²")
    
    return E_incident, cross_sections

if __name__ == "__main__":
    main()
