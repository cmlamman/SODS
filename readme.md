Functions for efficiently simulating 3D galaxy light profiles and exploring orientation-dependent selection effects.

*This is written specifically for generating a DESI LRG-like sample on NERSC from Abacus halo catalogs. Functions can be generalized with user edits.*


## Main Functions

1. `getting_abacus_halo-info.py`\
    Getting halo information (including shapes) from Abacus\
    Requirements:\ 
    CompaSOHaloCatalog from [abacusnbody](https://abacusutils.readthedocs.io/en/latest/compaso.html) `pip install abacusutils`

2. `map-catalog-to-sky.py`\
    Mapps comoving coordinates and 3D shapes to RA, DEC, Z, and projected shapes (E1, E2)\
    Requirements: 
    - [nbodykit](https://nbodykit.readthedocs.io/en/latest/getting-started/install.html):
    ```
    conda create --name nbodykit-env python=3 # or use python=2 for python 2.7*
    source activate nbodykit-env
    conda install -c bccp nbodykit
    ```

3. `simulating-selection.py`\
    Takes any catalog with 3D shapes and redshifts, generates columns for projected light profiles\
    *Note: this necessarily loops over every galaxy and will take a long time for large catalogs.*