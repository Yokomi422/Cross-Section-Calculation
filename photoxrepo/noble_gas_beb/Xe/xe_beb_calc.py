#!/usr/bin/env python3
"""
Calculate BEB ionization cross sections for Xe atom
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common import BEBCalculator, BEBPlotter
from pathlib import Path
import numpy as np
import re

def parse_xe_gaussian_improved(calc, filename):
    """
    Parse Gaussian output with IOp(6/80=1) for orbital energies and kinetic energies
    """
    with open(filename, 'r') as f:
        content = f.read()
    
    # 新しい形式の軌道エネルギーと運動エネルギーを探す
    pattern = r'Orbital energies and kinetic energies \(alpha\):\s*\n\s*1\s+2\s*\n((?:\s*\d+\s+\([^)]+\)--[OV]\s+[-\d.]+\s+[\d.]+\n)+)'
    match = re.search(pattern, content)
    
    if match:
        orbital_data = match.group(1)
        E_orb_list = []
        Ekin_orb_list = []
        occ_list = []
        
        # 各軌道の情報を抽出
        orbital_pattern = r'\s*(\d+)\s+\(([^)]+)\)--([OV])\s+([-\d.]+)\s+([\d.]+)'
        for orbital_match in re.finditer(orbital_pattern, orbital_data):
            orbital_type = orbital_match.group(3)
            if orbital_type == 'O':  # 占有軌道のみ
                E_orb = float(orbital_match.group(4))  # Hartree (negative value)
                Ekin_orb = float(orbital_match.group(5))  # Hartree
                symmetry = orbital_match.group(2)
                
                # 対称性から占有数を推定
                if 'S' in symmetry and 'G' in symmetry:
                    occ = 2  # s軌道
                elif 'P' in symmetry and 'U' in symmetry:
                    occ = 2  # p軌道の一つ（3つのp軌道で合計6電子）
                elif 'D' in symmetry:
                    occ = 2  # d軌道の一つ（5つのd軌道で合計10電子）
                else:
                    occ = 2  # デフォルト
                
                # BEBでは正の束縛エネルギーが必要
                E_orb_list.append(-E_orb)  # 正の束縛エネルギー
                Ekin_orb_list.append(Ekin_orb)  # 運動エネルギー
                occ_list.append(occ)
        
        calc.E_orb = E_orb_list
        calc.Ekin_orb = Ekin_orb_list
        calc.occ = occ_list
        
        print(f"\n# Improved parsing results for Xe:")
        for i, (E, K, occ) in enumerate(zip(calc.E_orb, calc.Ekin_orb, calc.occ)):
            print(f"# Orbital {i+1}: Binding E = {E:.6f} Hartree = {E * 27.2114:.2f} eV, Kinetic E = {K:.6f} Hartree = {K * 27.2114:.2f} eV, Occ = {occ}")
        
        return True
    
    return False

def main():
    # Set up paths
    current_dir = Path(__file__).parent
    output_dir = current_dir.parent / 'results' / 'Xe'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize calculator
    calc = BEBCalculator('Xe')
    
    # Parse Gaussian output (if available)
    gaussian_file = current_dir / 'xe_beb.log'
    if not gaussian_file.exists():
        # Use literature values for Xe
        print("# Warning: No Gaussian output found. Using literature values for Xe.")
        # Xe orbitals: 1s, 2s, 2p, 3s, 3p, 3d, 4s, 4p, 4d, 5s, 5p
        # Binding energies (eV) - approximate values
        calc.E_orb = [
            34561/27.2114,  # 1s
            5453/27.2114,   # 2s
            5104/27.2114, 5104/27.2114, 5104/27.2114,  # 2p
            1149/27.2114,   # 3s
            1002/27.2114, 1002/27.2114, 1002/27.2114,  # 3p
            689/27.2114, 689/27.2114, 689/27.2114, 689/27.2114, 689/27.2114,  # 3d
            213/27.2114,    # 4s
            146/27.2114, 146/27.2114, 146/27.2114,     # 4p
            67/27.2114, 67/27.2114, 67/27.2114, 67/27.2114, 67/27.2114,  # 4d
            23.4/27.2114,   # 5s
            12.1/27.2114, 12.1/27.2114, 12.1/27.2114   # 5p
        ]
        # Kinetic energies - estimated
        calc.Ekin_orb = [
            49700/27.2114,  # 1s
            7840/27.2114,   # 2s
            6800/27.2114, 6800/27.2114, 6800/27.2114,  # 2p
            1650/27.2114,   # 3s
            1340/27.2114, 1340/27.2114, 1340/27.2114,  # 3p
            820/27.2114, 820/27.2114, 820/27.2114, 820/27.2114, 820/27.2114,  # 3d
            300/27.2114,    # 4s
            195/27.2114, 195/27.2114, 195/27.2114,     # 4p
            80/27.2114, 80/27.2114, 80/27.2114, 80/27.2114, 80/27.2114,  # 4d
            27/27.2114,     # 5s
            13/27.2114, 13/27.2114, 13/27.2114         # 5p
        ]
        calc.occ = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]  # 1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6
    else:
        # Try different parsing approaches
        try:
            # First try the standard parser
            calc.parse_gaussian_output(gaussian_file)
            print("# Standard parsing successful")
        except ValueError:
            # If that fails, try our improved parser
            print("# Standard parsing failed. Trying improved parser...")
            if not parse_xe_gaussian_improved(calc, gaussian_file):
                # If that also fails, use manual values
                print("# Improved parsing failed. Using literature values...")
                calc.E_orb = [
                    34561/27.2114,  # 1s
                    5453/27.2114,   # 2s
                    5104/27.2114, 5104/27.2114, 5104/27.2114,  # 2p
                    1149/27.2114,   # 3s
                    1002/27.2114, 1002/27.2114, 1002/27.2114,  # 3p
                    689/27.2114, 689/27.2114, 689/27.2114, 689/27.2114, 689/27.2114,  # 3d
                    213/27.2114,    # 4s
                    146/27.2114, 146/27.2114, 146/27.2114,     # 4p
                    67/27.2114, 67/27.2114, 67/27.2114, 67/27.2114, 67/27.2114,  # 4d
                    23.4/27.2114,   # 5s
                    12.1/27.2114, 12.1/27.2114, 12.1/27.2114   # 5p
                ]
                calc.Ekin_orb = [
                    49700/27.2114,  # 1s
                    7840/27.2114,   # 2s
                    6800/27.2114, 6800/27.2114, 6800/27.2114,  # 2p
                    1650/27.2114,   # 3s
                    1340/27.2114, 1340/27.2114, 1340/27.2114,  # 3p
                    820/27.2114, 820/27.2114, 820/27.2114, 820/27.2114, 820/27.2114,  # 3d
                    300/27.2114,    # 4s
                    195/27.2114, 195/27.2114, 195/27.2114,     # 4p
                    80/27.2114, 80/27.2114, 80/27.2114, 80/27.2114, 80/27.2114,  # 4d
                    27/27.2114,     # 5s
                    13/27.2114, 13/27.2114, 13/27.2114         # 5p
                ]
                calc.occ = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    
    # Print summary
    calc.print_summary()
    
    # Calculate cross sections
    E_incident, cross_sections, orbital_contributions = calc.calculate_cross_sections()
    
    # Save results
    output_file = output_dir / 'xe_beb_results.dat'
    calc.save_results(E_incident, cross_sections, output_file)
    print(f"\n# Results saved to: {output_file}")
    
    # Create plots if BEBPlotter is available
    try:
        plotter = BEBPlotter('Xe')
        plotter.plot_cross_section(E_incident, cross_sections, output_dir)
        
        orbital_names = ['1s', '2s', '2p', '2p', '2p', '3s', '3p', '3p', '3p', 
                         '3d', '3d', '3d', '3d', '3d', '4s', '4p', '4p', '4p',
                         '4d', '4d', '4d', '4d', '4d', '5s', '5p', '5p', '5p']
        plotter.plot_orbital_contributions(E_incident, orbital_contributions, 
                                         orbital_names=orbital_names, 
                                         output_dir=output_dir)
        print(f"# Plots saved to: {output_dir}")
    except Exception as e:
        print(f"# Warning: Could not create plots: {e}")
    
    # 実験値との比較
    from common.beb_calculator import get_experimental_data
    exp_data = get_experimental_data('Xe')
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