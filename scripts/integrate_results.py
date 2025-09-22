#!/usr/bin/env python3
"""
Main Results Integration and Visualization Script for Metastasis in Osteosarcomas Analysis

This script orchestrates the generation of all visualization plots and figures
to summarize the analysis results for the osteosarcoma metastasis study.
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path

# Add the scripts directory to Python path
script_dir = Path(__file__).parent
sys.path.append(str(script_dir))

# Import visualization modules
try:
    from volcano_plot import create_volcano_plot
    from pca_plot import create_pca_plot
    from heatmap import create_heatmap, create_clustered_heatmap
    from pathway_diagrams import create_pathway_enrichment_plot, create_metastasis_pathway_diagram
except ImportError as e:
    print(f"Warning: Could not import visualization modules: {e}")
    print("Running as standalone script using subprocess calls...")

def check_dependencies():
    """Check if required Python packages are installed."""
    required_packages = [
        'pandas', 'matplotlib', 'seaborn', 'numpy', 
        'scipy', 'sklearn'
    ]
    
    missing_packages = []
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing_packages.append(package)
    
    if missing_packages:
        print(f"Missing required packages: {', '.join(missing_packages)}")
        print("Please install them using: pip install -r requirements.txt")
        return False
    return True

def run_script(script_name, args_list):
    """Run a visualization script with given arguments."""
    script_path = script_dir / script_name
    cmd = ['python', str(script_path)] + args_list
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"✓ {script_name} completed successfully")
        if result.stdout:
            print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ {script_name} failed: {e}")
        if e.stderr:
            print(f"Error: {e.stderr}")
        return False
    except FileNotFoundError:
        print(f"✗ Script {script_name} not found")
        return False

def create_all_visualizations(data_dir="data", output_dir="figures", use_imports=True):
    """
    Create all visualization plots for the analysis.
    
    Parameters:
    - data_dir: Directory containing input data files
    - output_dir: Directory to save output figures
    - use_imports: Whether to use direct imports or subprocess calls
    """
    
    # Define input files
    deseq_file = os.path.join(data_dir, "deseq2_results.csv")
    pca_file = os.path.join(data_dir, "pca_results.csv")
    expression_file = os.path.join(data_dir, "expression_matrix.csv")
    
    # Check if input files exist
    input_files = {
        "DESeq2 results": deseq_file,
        "PCA results": pca_file,
        "Expression matrix": expression_file
    }
    
    missing_files = []
    for name, filepath in input_files.items():
        if not os.path.exists(filepath):
            missing_files.append(f"{name}: {filepath}")
    
    if missing_files:
        print("Warning: Missing input files:")
        for missing in missing_files:
            print(f"  - {missing}")
        print("Using sample data or creating plots without missing files...")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print("Starting visualization generation...")
    print("="*50)
    
    success_count = 0
    total_plots = 0
    
    # 1. Volcano Plot
    total_plots += 1
    print(f"\n{total_plots}. Generating Volcano Plot...")
    if use_imports and os.path.exists(deseq_file):
        try:
            create_volcano_plot(deseq_file, output_dir)
            success_count += 1
            print("✓ Volcano plot created successfully")
        except Exception as e:
            print(f"✗ Volcano plot failed: {e}")
    elif os.path.exists(deseq_file):
        if run_script("volcano_plot.py", ["-i", deseq_file, "-o", output_dir]):
            success_count += 1
    else:
        print("✗ DESeq2 results file not found, skipping volcano plot")
    
    # 2. PCA Plot
    total_plots += 1
    print(f"\n{total_plots}. Generating PCA Plot...")
    pca_input = pca_file if os.path.exists(pca_file) else expression_file
    if use_imports and os.path.exists(pca_input):
        try:
            create_pca_plot(pca_input, output_dir)
            success_count += 1
            print("✓ PCA plot created successfully")
        except Exception as e:
            print(f"✗ PCA plot failed: {e}")
    elif os.path.exists(pca_input):
        if run_script("pca_plot.py", ["-i", pca_input, "-o", output_dir]):
            success_count += 1
    else:
        print("✗ PCA or expression data not found, skipping PCA plot")
    
    # 3. Expression Heatmap
    total_plots += 1
    print(f"\n{total_plots}. Generating Expression Heatmap...")
    if use_imports and os.path.exists(expression_file):
        try:
            deseq_arg = deseq_file if os.path.exists(deseq_file) else None
            create_heatmap(expression_file, deseq_arg, output_dir)
            success_count += 1
            print("✓ Heatmap created successfully")
        except Exception as e:
            print(f"✗ Heatmap failed: {e}")
    elif os.path.exists(expression_file):
        args = ["-e", expression_file, "-o", output_dir]
        if os.path.exists(deseq_file):
            args.extend(["-d", deseq_file])
        if run_script("heatmap.py", args):
            success_count += 1
    else:
        print("✗ Expression matrix not found, skipping heatmap")
    
    # 4. Clustered Heatmap
    total_plots += 1
    print(f"\n{total_plots}. Generating Clustered Heatmap...")
    if use_imports and os.path.exists(expression_file):
        try:
            deseq_arg = deseq_file if os.path.exists(deseq_file) else None
            create_clustered_heatmap(expression_file, deseq_arg, output_dir)
            success_count += 1
            print("✓ Clustered heatmap created successfully")
        except Exception as e:
            print(f"✗ Clustered heatmap failed: {e}")
    elif os.path.exists(expression_file):
        args = ["-e", expression_file, "-o", output_dir, "-c"]
        if os.path.exists(deseq_file):
            args.extend(["-d", deseq_file])
        if run_script("heatmap.py", args):
            success_count += 1
    else:
        print("✗ Expression matrix not found, skipping clustered heatmap")
    
    # 5. Pathway Enrichment Plot
    total_plots += 1
    print(f"\n{total_plots}. Generating Pathway Enrichment Plot...")
    if use_imports and os.path.exists(deseq_file):
        try:
            create_pathway_enrichment_plot(deseq_file, output_dir)
            success_count += 1
            print("✓ Pathway enrichment plot created successfully")
        except Exception as e:
            print(f"✗ Pathway enrichment plot failed: {e}")
    elif os.path.exists(deseq_file):
        if run_script("pathway_diagrams.py", ["-d", deseq_file, "-o", output_dir]):
            success_count += 1
    else:
        print("✗ DESeq2 results not found, skipping pathway enrichment plot")
    
    # 6. Metastasis Pathway Diagram
    total_plots += 1
    print(f"\n{total_plots}. Generating Metastasis Pathway Diagram...")
    if use_imports:
        try:
            create_metastasis_pathway_diagram(output_dir)
            success_count += 1
            print("✓ Metastasis pathway diagram created successfully")
        except Exception as e:
            print(f"✗ Metastasis pathway diagram failed: {e}")
    else:
        if run_script("pathway_diagrams.py", ["--diagram-only", "-o", output_dir]):
            success_count += 1
    
    # Summary
    print("\n" + "="*50)
    print("VISUALIZATION SUMMARY")
    print("="*50)
    print(f"Total plots attempted: {total_plots}")
    print(f"Successfully created: {success_count}")
    print(f"Failed: {total_plots - success_count}")
    print(f"\nOutput directory: {os.path.abspath(output_dir)}")
    
    # List generated files
    if os.path.exists(output_dir):
        output_files = [f for f in os.listdir(output_dir) if f.endswith(('.png', '.pdf'))]
        if output_files:
            print(f"\nGenerated files:")
            for file in sorted(output_files):
                print(f"  - {file}")
        else:
            print("\nNo output files found.")
    
    return success_count, total_plots

def main():
    parser = argparse.ArgumentParser(
        description='Generate all visualization plots for Metastasis in Osteosarcomas analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                                    # Use default directories
  %(prog)s --data ./data --output ./figures  # Specify custom directories
  %(prog)s --check-deps                      # Check dependencies only
        """
    )
    
    parser.add_argument('--data', '-d', default='data', 
                       help='Directory containing input data files (default: data)')
    parser.add_argument('--output', '-o', default='figures', 
                       help='Output directory for figures (default: figures)')
    parser.add_argument('--check-deps', action='store_true', 
                       help='Check dependencies and exit')
    parser.add_argument('--subprocess', action='store_true', 
                       help='Force use of subprocess calls instead of direct imports')
    
    args = parser.parse_args()
    
    if args.check_deps:
        if check_dependencies():
            print("✓ All required dependencies are installed")
            return 0
        else:
            return 1
    
    # Check dependencies
    if not check_dependencies():
        print("Please install missing dependencies before proceeding.")
        return 1
    
    # Set up paths
    data_dir = args.data
    output_dir = args.output
    
    print("Metastasis in Osteosarcomas - Results Integration and Visualization")
    print("="*70)
    print(f"Data directory: {os.path.abspath(data_dir)}")
    print(f"Output directory: {os.path.abspath(output_dir)}")
    
    # Create visualizations
    use_imports = not args.subprocess
    success_count, total_plots = create_all_visualizations(data_dir, output_dir, use_imports)
    
    # Return appropriate exit code
    return 0 if success_count == total_plots else 1

if __name__ == "__main__":
    exit(main())