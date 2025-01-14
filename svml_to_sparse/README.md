
# svml_to_sparse

This project contains code for converting `.svml` graphlet output per protein into one big `.npz` SciPy sparse matrix format.

## Folder Structure

- **graphlet_svml_to_sparse.py**: A Python script that performs the conversion in one go.  
- **graphlet_svml_to_sparse.ipynb**: A Jupyter notebook for incremental steps of the conversion process.
- **example_1arz_A.svml**: An example `.svml` file that would be outputted from the graphlet counting for a single protein.

## Required Libraries

The following libraries are required to run the scripts. I used Python 3.11.8:

```python
import os
import glob
import pandas as pd
import numpy as np
import sys
from scipy.sparse import coo_matrix, csr_matrix, save_npz
```

## Usage

```bash
python graphlet_svml_to_sparse.py <path_to_svml_dir> <path_to_results_dir>
```

Where:
- `<path_to_svml_dir>` is the directory containing the `.svml` files. They should be named 
- `<path_to_results_dir>` is the directory where the resulting `.npz` matrix file will be saved, along with a labels file and optional pdb ids file (useful for per-protein stratification).
