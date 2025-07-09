#!/usr/bin/env python3
"""
Run BEB calculations for all noble gases and create comparison plots
"""

import sys
import subprocess
from pathlib import Path
import matplotlib.pyplot as plt

# Add common directory to path
sys.path.append(str(Path(__file__).parent))
from common import BEBPlotter


def run_calculation(element):
    """Run BEB calculation for a single element"""
    calc_script = Path(__file__).parent / element / f"{element.lower()}_beb_calc.py"
    
    if not calc_script.exists():
        print(f"Warning: Calculation script not found for {element}")
        return False
    
    print(f"\n{'='*60}")
    print(f"Running BEB calculation for {element}")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run([sys.executable, str(calc_script)], 
                              capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print("Errors:", result.stderr)
        return result.returncode == 0
    except Exception as e:
        print(f"Error running {element} calculation: {e}")
        return False


def load_results(element):
    """Load calculation results for an element"""
    results_file = Path(__file__).parent / 'results' / element / f'{element.lower()}_beb_results.dat'
    
    if not results_file.exists():
        print(f"Warning: Results file not found for {element}")
        return None, None
    
    energy = []
    cross_section = []
    
    with open(results_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                energy.append(float(parts[0]))
                cross_section.append(float(parts[1]))
    
    return energy, cross_section


def create_comparison_plot():
    """Create comparison plot for all noble gases"""
    elements = ['He', 'Ne', 'Ar', 'Kr', 'Xe']
    elements_data = {}
    
    # Load all results
    for element in elements:
        energy, cross_section = load_results(element)
        if energy and cross_section:
            elements_data[element] = {
                'energy': energy,
                'cross_section': cross_section
            }
    
    if not elements_data:
        print("No results found to plot")
        return
    
    # Create comparison plot
    plotter = BEBPlotter('Ar')  # Just for access to the method
    output_file = Path(__file__).parent / 'results' / 'noble_gas_comparison.png'
    output_file.parent.mkdir(exist_ok=True)
    
    plotter.create_comparison_plot(elements_data, output_file)
    
    # Also create individual summary plot
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    for i, element in enumerate(elements):
        if element in elements_data:
            ax = axes[i]
            data = elements_data[element]
            
            # Create individual plotter for experimental data
            element_plotter = BEBPlotter(element)
            
            # Plot calculation
            ax.loglog(data['energy'], data['cross_section'], 'b-', linewidth=2, label='BEB')
            
            # Plot experimental data if available
            if element_plotter.exp_data:
                ax.loglog(element_plotter.exp_data['energy'], 
                         element_plotter.exp_data['cross_section'], 
                         'ro', markersize=6, label='Experiment')
            
            # Add threshold
            if element_plotter.threshold:
                ax.axvline(x=element_plotter.threshold, color='gray', 
                          linestyle='--', alpha=0.5)
            
            ax.set_xlabel('Energy (eV)')
            ax.set_ylabel('Cross Section (Å²)')
            ax.set_title(f'{element}')
            ax.grid(True, which="both", ls="-", alpha=0.2)
            ax.legend(loc='upper right', fontsize=8)
            ax.set_xlim(10, 2000)
    
    # Hide the last subplot if we have 5 elements
    if len(elements) < 6:
        axes[-1].set_visible(False)
    
    plt.suptitle('Noble Gas Ionization Cross Sections: BEB Model', fontsize=16)
    plt.tight_layout()
    
    summary_file = Path(__file__).parent / 'results' / 'noble_gas_summary.png'
    plt.savefig(summary_file, dpi=300, bbox_inches='tight')
    print(f"\nSummary plot saved to: {summary_file}")


def main():
    """Main workflow"""
    print("Noble Gas BEB Ionization Cross Section Calculations")
    print("="*60)
    
    elements = ['He', 'Ne', 'Ar', 'Kr', 'Xe']
    
    # Run calculations for each element
    success_count = 0
    for element in elements:
        if run_calculation(element):
            success_count += 1
    
    print(f"\n{'='*60}")
    print(f"Completed {success_count}/{len(elements)} calculations successfully")
    
    # Create comparison plots
    if success_count > 0:
        print("\nCreating comparison plots...")
        create_comparison_plot()
    
    print("\n" + "="*60)
    print("Workflow completed!")
    print("Results are saved in the 'results' directory")
    print("Individual plots are in results/<element>/ directories")
    print("Comparison plots are in the results/ directory")
    print("="*60)


if __name__ == "__main__":
    main()