#!/usr/bin/env python3
"""
Common BEB (Binary Encounter Bethe) calculator for noble gases
"""

import numpy as np
import re
import sys
from pathlib import Path

AU2EV = 27.2114
ANG2BOHR = 1.889726132873


class BEBCalculator:
    """Calculate BEB ionization cross sections for atoms"""
    
    def __init__(self, element_name, orbital_config=None):
        self.element = element_name
        self.orbital_config = orbital_config or self._get_default_config()
        self.E_orb = []
        self.Ekin_orb = []
        
    def _get_default_config(self):
        """Get default orbital configuration for noble gases"""
        configs = {
            'He': {'electrons': 2, 'orbitals': ['1s']},
            'Ne': {'electrons': 10, 'orbitals': ['1s', '2s', '2p']},
            'Ar': {'electrons': 18, 'orbitals': ['1s', '2s', '2p', '3s', '3p']},
            'Kr': {'electrons': 36, 'orbitals': ['1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p']},
            'Xe': {'electrons': 54, 'orbitals': ['1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p', '4d', '5s', '5p']},
        }
        return configs.get(self.element, {})
    
    def parse_gaussian_output(self, filename):
        """Parse Gaussian output file to extract orbital energies and kinetic energies"""
        self.E_orb = []
        self.Ekin_orb = []
        
        start_line = ' Orbital energies and kinetic energies (alpha):\n'
        dec = r' +-?\d+\.\d+'
        reg = re.compile(r'^ +[0-9]+ +O'+dec+dec+r'$')
        
        with open(filename, 'r') as f:
            read = False
            for line in f:
                if line == start_line:
                    read = True
                res = reg.search(line)
                if res and read:
                    # Convert to positive binding energy (Hartree)
                    self.E_orb.append(-float(line.split()[2]))
                    self.Ekin_orb.append(float(line.split()[3]))
        
        if not self.E_orb:
            raise ValueError(f"Could not find orbital energies in {filename}. "
                           "Make sure you used IOp(6/80=1) in your calculation")
        
        return self.E_orb, self.Ekin_orb
    
    def beb_cross_section(self, T, B, U, N=2):
        """
        Calculate BEB cross section for a single orbital
        T: kinetic energy of incident electron (au)
        B: binding energy of orbital (au)
        U: kinetic energy of orbital (au)
        N: occupation number (2 for closed shell)
        """
        a0 = 1  # Bohr radius in au
        R = 0.5  # Rydberg energy in au
        
        t = T / B
        u = U / B
        
        S = 4 * np.pi * a0**2 * N * (R/B)**2
        
        x1 = S / (t + u + 1)
        x2 = np.log(t) / 2 * (1 - 1/t**2)
        x3 = 1 - 1/t - np.log(t)/(1+t)
        
        sigma = x1 * (x2 + x3)
        
        return sigma
    
    def calculate_total_cross_section(self, T_eV):
        """Calculate total ionization cross section for all orbitals"""
        T = T_eV / AU2EV
        
        total_sigma = 0.0
        orbital_sigmas = []
        
        for B, U in zip(self.E_orb, self.Ekin_orb):
            if T > B:
                sigma = self.beb_cross_section(T, B, U)
                sigma_ang2 = sigma / (ANG2BOHR**2)
                orbital_sigmas.append(sigma_ang2)
                total_sigma += sigma_ang2
            else:
                orbital_sigmas.append(0.0)
        
        return total_sigma, orbital_sigmas
    
    def calculate_cross_sections(self, E_min=None, E_max=1000.0, n_points=200):
        """Calculate cross sections over energy range"""
        if E_min is None:
            E_min = min(self.E_orb) * AU2EV + 0.1
        
        E_incident = np.logspace(np.log10(E_min), np.log10(E_max), n_points)
        cross_sections = []
        orbital_contributions = []
        
        for E in E_incident:
            total_sigma, orbital_sigmas = self.calculate_total_cross_section(E)
            cross_sections.append(total_sigma)
            orbital_contributions.append(orbital_sigmas)
        
        return E_incident, np.array(cross_sections), orbital_contributions
    
    def save_results(self, E_incident, cross_sections, output_file):
        """Save results to file"""
        with open(output_file, 'w') as f:
            f.write(f"# BEB ionization cross sections for {self.element}\n")
            f.write("# Energy[eV]  Total_cross_section[Ang^2]\n")
            for E, sigma in zip(E_incident, cross_sections):
                f.write(f"{E:12.6f}  {sigma:12.6e}\n")
    
    def print_summary(self):
        """Print summary of orbital energies"""
        print(f"\n# {self.element} Orbital Summary:")
        print(f"# Found {len(self.E_orb)} occupied orbitals")
        print("# Orbital binding energies (eV):", 
              [f"{e*AU2EV:.2f}" for e in self.E_orb])
        print("# Orbital kinetic energies (eV):", 
              [f"{e*AU2EV:.2f}" for e in self.Ekin_orb])
        print(f"# Ionization threshold: {min(self.E_orb)*AU2EV:.2f} eV")


def get_experimental_data(element):
    """Get experimental ionization cross section data for noble gases"""
    data = {
        'He': {
            'energy': np.array([25, 30, 40, 50, 60, 80, 100, 150, 200, 300, 400, 600, 800, 1000]),
            'cross_section': np.array([0.354, 0.351, 0.332, 0.307, 0.283, 0.244, 0.214, 0.162, 0.131, 0.094, 0.074, 0.051, 0.039, 0.032]),
            'reference': 'Shah et al., J. Phys. B 20, 3501 (1987)'
        },
        'Ne': {
            'energy': np.array([22, 25, 30, 40, 50, 60, 80, 100, 150, 200, 300, 400, 600, 800, 1000]),
            'cross_section': np.array([0.30, 0.50, 0.72, 0.96, 1.05, 1.08, 1.04, 0.97, 0.78, 0.64, 0.47, 0.37, 0.26, 0.20, 0.16]),
            'reference': 'Krishnakumar & Srivastava, J. Phys. B 21, 1055 (1988)'
        },
        'Ar': {
            'energy': np.array([16, 18, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 300, 400, 500, 600, 800, 1000]),
            'cross_section': np.array([0.20, 0.52, 0.88, 1.60, 2.10, 2.42, 2.62, 2.75, 2.82, 2.85, 2.82, 2.75, 2.66, 2.56, 2.35, 2.16, 2.00, 1.86, 1.73, 1.49, 1.31, 1.05, 0.88, 0.76, 0.60, 0.50]),
            'reference': 'Rejoub et al., Phys. Rev. A 65, 042713 (2002)'
        },
        'Kr': {
            'energy': np.array([25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 250, 300, 400, 500, 600, 800, 1000]),
            'cross_section': np.array([1.0, 2.3, 3.2, 3.9, 4.6, 4.9, 5.0, 4.95, 4.85, 4.70, 4.35, 3.85, 3.15, 2.65, 2.30, 1.80, 1.50, 1.28, 1.00, 0.82]),
            'reference': 'Rejoub et al., Phys. Rev. A 65, 042713 (2002)'
        },
        'Xe': {
            'energy': np.array([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 250, 300, 400, 500, 600, 800, 1000]),
            'cross_section': np.array([0.8, 2.0, 3.2, 4.2, 5.0, 5.5, 6.2, 6.5, 6.6, 6.55, 6.45, 6.3, 5.9, 5.3, 4.4, 3.7, 3.2, 2.55, 2.10, 1.80, 1.40, 1.15]),
            'reference': 'Rejoub et al., Phys. Rev. A 65, 042713 (2002)'
        }
    }
    return data.get(element, None)


def get_ionization_thresholds():
    """Get ionization thresholds for noble gases (eV)"""
    return {
        'He': 24.59,
        'Ne': 21.56,
        'Ar': 15.76,
        'Kr': 14.00,
        'Xe': 12.13
    }