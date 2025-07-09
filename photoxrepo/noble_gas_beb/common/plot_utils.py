#!/usr/bin/env python3
"""
Common plotting utilities for BEB cross sections
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from .beb_calculator import get_experimental_data, get_ionization_thresholds


class BEBPlotter:
    """Plotting utilities for BEB cross sections"""
    
    def __init__(self, element):
        self.element = element
        self.exp_data = get_experimental_data(element)
        self.threshold = get_ionization_thresholds().get(element)
        
    def plot_cross_section(self, E_calc, sigma_calc, output_dir=None, show_exp=True):
        """Create publication-quality plot of cross sections"""
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Convert from Å² to m²
        sigma_calc_m2 = sigma_calc * 1e-20
        
        # Plot calculated results
        ax.loglog(E_calc, sigma_calc_m2, 'b-', linewidth=2, label=f'{self.element} BEB calculation')
        
        # Plot experimental data if available
        if show_exp and self.exp_data:
            exp_sigma_m2 = self.exp_data['cross_section'] * 1e-20
            ax.loglog(self.exp_data['energy'], exp_sigma_m2, 
                     'ro', markersize=8, label=f'Experiment ({self.exp_data["reference"].split(",")[0]})')
        
        # Add ionization threshold
        if self.threshold:
            ax.axvline(x=self.threshold, color='gray', linestyle='--', alpha=0.5, 
                      label=f'Ionization threshold ({self.threshold:.2f} eV)')
        
        # Mark maximum cross section
        max_sigma = np.max(sigma_calc_m2)
        max_E = E_calc[np.argmax(sigma_calc_m2)]
        ax.axvline(x=max_E, color='green', linestyle=':', alpha=0.5)
        ax.text(max_E*1.1, max_sigma*0.8, f'Max: {max_sigma:.2e} m² at {max_E:.1f} eV',
                fontsize=10, color='green')
        
        # Formatting
        ax.set_xlabel('Incident Electron Energy (eV)', fontsize=14)
        ax.set_ylabel('Total Ionization Cross Section (m²)', fontsize=14)
        ax.set_title(f'Electron Impact Ionization Cross Section of {self.element}', fontsize=16)
        ax.grid(True, which="both", ls="-", alpha=0.2)
        ax.legend(loc='upper right', fontsize=12)
        
        # Set appropriate axis limits
        ax.set_xlim(10, 2000)
        y_max = max(10e-20, max_sigma * 2)
        ax.set_ylim(1e-22, y_max)
        
        # Add calculation details
        textstr = f'BEB Model\n{self.element} atom\nB3LYP/aug-cc-pVTZ'
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.15, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)
        
        # Save figure
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(exist_ok=True)
            
            for fmt in ['png', 'pdf', 'svg']:
                output_file = output_dir / f'{self.element.lower()}_beb_cross_section.{fmt}'
                fig.savefig(output_file, dpi=300, bbox_inches='tight')
                print(f"# Plot saved to: {output_file}")
        
        return fig, ax
    
    def plot_orbital_contributions(self, E_calc, orbital_contributions, orbital_names=None, output_dir=None):
        """Plot orbital contributions to ionization cross section"""
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Total cross section (convert to m²)
        total_sigma = np.sum(orbital_contributions, axis=1) * 1e-20
        ax.semilogy(E_calc, total_sigma, 'k-', linewidth=2, label='Total')
        
        # Individual orbital contributions
        colors = plt.cm.tab10(np.linspace(0, 1, len(orbital_contributions[0])))
        
        for i, color in enumerate(colors[:len(orbital_contributions[0])]):
            orbital_sigma = np.array([contrib[i] for contrib in orbital_contributions]) * 1e-20
            label = orbital_names[i] if orbital_names and i < len(orbital_names) else f'Orbital {i+1}'
            ax.semilogy(E_calc, orbital_sigma, '--', color=color, label=label)
        
        ax.set_xlabel('Incident Electron Energy (eV)', fontsize=14)
        ax.set_ylabel('Cross Section (m²)', fontsize=14)
        ax.set_title(f'Orbital Contributions to {self.element} Ionization Cross Section', fontsize=16)
        ax.grid(True, which="both", ls="-", alpha=0.2)
        ax.legend(loc='upper right', fontsize=10)
        ax.set_xlim(10, 1000)
        
        if output_dir:
            output_dir = Path(output_dir)
            output_file = output_dir / f'{self.element.lower()}_orbital_contributions.png'
            fig.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"# Orbital contribution plot saved to: {output_file}")
        
        return fig, ax
    
    def create_comparison_plot(self, elements_data, output_file=None):
        """Create comparison plot for multiple elements"""
        fig, ax = plt.subplots(figsize=(12, 8))
        
        colors = ['blue', 'red', 'green', 'orange', 'purple']
        
        for i, (element, data) in enumerate(elements_data.items()):
            E_calc = data['energy']
            sigma_calc = data['cross_section'] * 1e-20  # Convert to m²
            color = colors[i % len(colors)]
            
            # Plot calculation
            ax.loglog(E_calc, sigma_calc, '-', color=color, linewidth=2, 
                     label=f'{element} BEB')
            
            # Plot experimental data if available
            exp_data = get_experimental_data(element)
            if exp_data:
                exp_sigma_m2 = exp_data['cross_section'] * 1e-20  # Convert to m²
                ax.loglog(exp_data['energy'], exp_sigma_m2, 
                         'o', color=color, markersize=6, alpha=0.7,
                         label=f'{element} Exp.')
        
        ax.set_xlabel('Incident Electron Energy (eV)', fontsize=14)
        ax.set_ylabel('Total Ionization Cross Section (m²)', fontsize=14)
        ax.set_title('Noble Gas Ionization Cross Sections: BEB vs Experiment', fontsize=16)
        ax.grid(True, which="both", ls="-", alpha=0.2)
        ax.legend(loc='upper right', fontsize=10, ncol=2)
        ax.set_xlim(10, 2000)
        ax.set_ylim(1e-22, 1e-19)
        
        if output_file:
            fig.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"# Comparison plot saved to: {output_file}")
        
        return fig, ax