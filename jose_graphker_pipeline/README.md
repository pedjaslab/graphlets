
# Jose GraphKer Pipeline

This repository provides tools and scripts for working with graphlet kernels based on structural data.

## Folder Structure

- `pdb_to_sgk_input.ipynb`: A Jupyter Notebook that takes a PDB file as input and creates the required input files for the graphlet kernel. Modify this notebook to suit your specific data needs.
- `graphlet_kernels_framework_v1.0/`: Contains the scripts and resources for generating graphlets.

## Usage

### Generating Graphlets for Undirected Edges
1. Navigate to the graphlet kernel framework directory:
   ```bash
   cd graphlet_kernels_framework_v1.0/graphlet_kernels_v1.0/
   ```
2. Run the standard graphlet kernel script:
   ```bash
   ./run_std.sh
   ```
   - Adjust the paths in `run_std.sh` to match your data locations.
   - Refer to `run_example.sh` and the `README` file in the `graphlet_kernels_v1.0` directory for additional usage instructions.

## Notes
- Ensure all dependencies and environment configurations are set up correctly before running the scripts.
- Feel free to customize the provided scripts and notebooks for your specific use case.
