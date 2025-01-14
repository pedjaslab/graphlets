
# svml_to_sparse

This project contains code for converting `.svml` graphlet output per protein into one big `.npz` SciPy sparse matrix format.

## Folder Structure

- **graphlet_svml_to_sparse.py**: A Python script that performs the conversion in one go.  
  **Usage**:  
  ```bash
  python graphlet_svml_to_sparse.py <path_to_svml_dir> <path_to_results_dir>
  ```
- **graphlet_svml_to_sparse.ipynb**: A Jupyter notebook for incremental steps of the conversion process.
- **example_1arz_A.svml**: An example `.svml` file that would be outputted from the graphlet counting for a single protein.

## Required Libraries

The following libraries are required to run the scripts:

```python
import os
import glob
import pandas as pd
import numpy as np
import sys
from scipy.sparse import coo_matrix, csr_matrix, save_npz
```

## Python Version

This code was developed using Python 3.11.8.

## Usage

To convert your `.svml` files into a `.npz` sparse matrix, use the following command:

```bash
python graphlet_svml_to_sparse.py <path_to_svml_dir> <path_to_results_dir>
```

Where:
- `<path_to_svml_dir>` is the directory containing the `.svml` files.
- `<path_to_results_dir>` is the directory where the resulting `.npz` matrix file will be saved.

## Example

You can use the provided `example_1arz_A.svml` as a sample input to see how the conversion works.

## Notes

- The code will read all `.svml` files in the provided directory and generate a sparse matrix stored in `.npz` format.
