
# Jose's Graphlet Kernel Pipeline

This repository provides tools and scripts for working with graphlet kernels from PDBs.

## Folder Structure

- `pdb_to_sgk_input.ipynb`: A Jupyter Notebook that takes a PDB file as input and creates the required input files for the graphlet kernel. Modify to suit your specific data needs.
- `graphlet_kernels_framework_v1.0/`: Contains Jose's scripts and resources for generating graphlets.

## Usage once input files are generated

### Generating Graphlets for Undirected Edges
1. Navigate to the graphlet kernel framework directory:
   ```bash
   cd graphlet_kernels_framework_v1.0/graphlet_kernels_v1.0/
   ```
2. Adjust the paths in `run_std.sh` to match your data locations.
3. Run the standard graphlet kernel script:
   ```bash
   ./run_std.sh
   ```
   - Refer to `run_example.sh` and the `README.txt` file in `graphlet_kernels_framework_v1.0/graphlet_kernels_v1.0` for additional usage instructions, e.g. label substitution, edge indel, or edit distance kernels.

## Notes
- Ensure all dependencies and environment configurations are set up correctly before running the input file generation:
```python
import numpy as np
import scipy.sparse as sp
from scipy.io import savemat
import os
from Bio import PDB
```
