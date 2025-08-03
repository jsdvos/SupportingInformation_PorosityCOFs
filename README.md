This file relates to the research done in the following paper:

Rasha G. AbdulHalim, Bikash Garai, Juul S. De Vos, Sander Borgmans, Lydia Gkoura, Sabu Varghese, Farah Benyettou, Mark A. Olson, Sven M. J. Rogge, Ali Trabolsi, _Maximizing Porosity and Water Sorption in Covalent Organic Frameworks via β-Ketoenamine Linkages_ (2025).

This file is part of the midterm storage (publically available under the [CC BY-SA license](https://creativecommons.org/licenses/by-sa/4.0/)) of the input files relevant in this work. In the remainder of this file, we outline in detail the workflow used in this work including references to the relevant input and output files.

# Software

The following software packages are used to perform all relevant calculations.

- Gaussian (Gaussian 16, Revision C.01)
- HORTON (version 2.0.0)
- Molden (version 6.8)
- QuickFF (version 2.2.4)
- TAMkin (version 1.2.6)
- Yaff (version 1.6.0)
- Zeo++ (version 0.3)
- GPXRDpy


# Workflow

Some tar-files are splitted in multiple smaller files. In these cases, the original file can be restored by running:

`cat file.tar.gz.* | tar xvzf -`
`cat input.tar.xz.a* | tar xJf -`

## STEP 1 - Cluster force field development

### Step 1a - QuickFF cluster force fields

The QuickFF force fields for the clusters visualized in Fig. S7 are derived using the protocol described in [this article](https://doi.org/10.1039/D3TA00470H). For more details on the workflow, we refer to the data archive of the article.

For each cluster, the _ab initio_ Hessian matrix is calculated using Gaussian. These can be found in the `SBU_freq.fchk` files. From this Hessian, the force field parameters, which can be found in the `pars_yaff.txt`, `pars_ei.txt`, and `pars_mm3.txt` files, are derived with QuickFF. Finally, the system is optimized using the Yaff software and stored in the `system_opt.chk` file.

### Step 1b - Additional force field terms

As explained in the SI, the default torsional terms do not properly describe the out-of-plane behavior of the imine and the amine dihedral angles. Therefore, an _ab initio_ rotation scan was performed to fit a more accurate dihedral term, replacing the original torsion term. The [pyiron](https://pyiron.org/) workflow is outlined in the `rotational_barrier_rasha.ipynb` Python Notebook. The `plot.py` script can help to reproduce Fig. S9.

## STEP 2 - Structure generation

### Step 2a - Initial structure

Starting from a single layer structure of COF-TP-100, eleven new structures are generated using the `create_init.py` script. Each structure contains a unit cell with 12 layers and a total of 96 three-connected SBUs (either TFB or TP). The suffix in the names of the structures indicates the number of resulting TFB units. An overview of the structure names and corresponding TP content is given below:

- ABCDEF_none: 100.0% TP
- ABCDEF_full10: 89.6% TP
- ABCDEF_full19: 80.2% TP
- ABCDEF_full29: 69.8% TP
- ABCDEF_full38: 60.4% TP
- ABCDEF_full48: 50.0% TP
- ABCDEF_full58: 39.6% TP
- ABCDEF_full67: 30.2% TP
- ABCDEF_full77: 19.8% TP
- ABCDEF_full86: 10.4% TP
- ABCDEF_all: 0.0% TP

### Step 2b - Optimization

Once an initial structure is generated and stored in an initial `.chk` file, it is relaxed using the earlier derived force field parameters (collected in the `pars.txt` file) during a Yaff optimization. The optimization can be executed by running one of the following commands:

`qsub opt.sh`
`python opt.py`

The optimized structure is stored in the `struct_opt.chk` file.

## STEP 3 - Molecular dynamics simulations

### Step 3a - Running the MD simulation

For each structure, an MD simulation is run on four different temperatures (77K, 293K, 400K, and 500K):

**input**
`init.chk`, `pars.txt`

**command line**
`qsub md.sh` or `python md.py`

**output**
`traj.h5`, `md.log`

### Step 3b - Analyzing the MD simulation

For each MD simulation, 500 snapshots are extracted from the last 75 ps of the simulation using the `extract_frames.py` script. Three sets of characterizations are performed on each snapshot:

- Geometrical characterization using Zeo++ of the accessible surface area and pore volume.
- Calculation of the PXRD pattern using GPXRDpy.
- Determination of the layer offset and interlayer distance.

Whereas all resulting log and data files are stored in the `calculations` folder, the relevant data per MD run is collected in some files in the `data` folder. Subsequently, Fig. S10, S11, and S12 can be reproduced by running `plot_pxrd.py`, `plot_geo.py`, and `plot_layer.py`, respectively.

## STEP 4 - Static energy scan

As explained in the SI, a static energy scan of the interlayer distance is performed for COF-TP-0 (ABCDEF_all) and COF-TP-100 (ABCDEF_none). From the twelve layered structure, two layers are extracted and subsequently optimized using the `prepare_system.py` script. By running the `plot.py` script, Fig. S13 can be reproduced.
