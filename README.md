# Molecular Vibration Explorer

**Molecular Vibration Explorer** is an interactive tool for exploring molecular vibrational spectra and tensorial light-vibration coupling strength for applications in the field of surface-enhanced spectroscopy. The Gold database gathers the results from density functional theory calculations on 2’800 commercially available thiol compounds linked to a gold atom, with the main motivation to screen the best molecules for THz and mid-infrared to visible upconversion. The Thiol database contains calculations on 1’900 commercially available thiol compounds. The different tools available to analyze the database were previously introduced in the following references: 
* [(1) Z. Koczor-Benda, A. L. Boehmke, A. Xomalis, R. Arul, C. Readman, J. J. Baumberg & E. Rosta "Molecular Screening for Terahertz Detection with Machine-Learning-Based Methods", Phys. Rev. X, 11(4):041035 (2021)](https://link.aps.org/doi/10.1103/PhysRevX.11.041035)
* [(2) P. Roelli, D. Martin-Cano, T.J. Kippenberg & C. Galland, "Molecular Platform for Frequency Upconversion at the Single-Photon Level", Phys. Rev. X, 10(3):031057 (2020)](https://link.aps.org/doi/10.1103/PhysRevX.10.031057)


## Description of jupyter notebooks

* **Database scan**

The *Gold* and *Thiol* databases can be explored separately. The database notebook allows for ranking molecules within the database and selecting them for further analysis. The user may choose the *target property* (A, R, or P) and set a frequency range of interest. The direction of the field polarization vectors can be specified independently, among three possible orthogonal directions. 

The distribution of the target property across the entire database, for mode frequencies within the user-specified range, is plotted as a histogram, and a table is generated below that shows all molecules sorted according to decreasing value of the target. Further properties of the molecules and their individual normal modes can be explored by navigating the links displayed in the table.


* **Molecule analysis**

In  the  tab  labeled  “Set molecule orientation” on the left, the molecular orientation can be modified at will by specifying the rotational angles with respect to the Cartesian coordinate system defined along the possible polarisation vectors. The spectra for the selected orientation can be directly compared with those for the orientation average by checking “Show full orientation average” in the “Set plotting parameters” tab on the left. For Raman scattering and SFG, the temperature and laser wavelength can be set in  the “Set experimental parameters” tab, while the frequency range, scaling factor, style of the spectrum (stick, broadened, stick+broadened) and corresponding broadening parameters can all be tuned in the “Set plotting  parameters” tab.  Finally, field polarization directions can be specified in the “Set field polarization” tab. 

In addition to optical properties, physical properties of the molecule relevant to surface-enhanced spectroscopy are listed below the 3D drawing. The length of the projection of the molecule along the z-axis is referred to as “molecule  height”. Indeed, it corresponds to the approximate thickness of a molecular monolayer assuming binding through the thiol group onto a flat gold surface spanning the x-y plane.  Finally, the anisotropy of the polarizability tensor is displayed, as it was found to be a good indication of the sensitivity of the Raman and SFG signals to molecular orientation. 

Another feature that is especially useful for navigating the database is the table of similar molecules that is displayed at the bottom of the notebook. Molecules susceptible to form self-assembled monolayers (SAMs) and similar to the chosen molecule are listed to help optimisation of surface-enhanced applications. 

* **Vibrational mode inspection**

The user can explore all normal modes of the selected molecule through 3D visualization.  The dependence of the IR, Raman, and conversion/SFG intensity of a specific normal mode on molecular orientation is represented by a projection onto the plane perpendicular to  the  axis  of  rotation. The  direction  of  the  field  polarization vectors of the IR and Raman (excitation and scattered) beams can also be set in this notebook.


## Methods

DFT calculations are performed at the B3LYP+D3/def2-SVP level using the Gaussian program package. Analytic formulas for the orientation averages of intensities were derived using Mathematica.
The Molecular Vibration Explorer web application is created using Voila, and is based on interactive Jupyter notebooks.
Molecular similarity scores are generated using RDKit. NGLview is used for the 3D visualization of molecules.


## How to cite

If this database was useful for your work, please cite it in the following way:
[Z. Koczor-Benda, P. Roelli, C. Galland & E. Rosta, “Molecular Vibration Explorer: an online database and toolbox for surface-enhanced frequency conversion and infrared and Raman spectroscopy”, J. Phys. Chem. A 2022, 126, 28, 4657–4663](https://doi.org/10.1021/acs.jpca.2c03700)


## How to go further / Advanced options

Both the DFT results used in this app and the notebooks scripts can be found following the link to the dedicated GitHub repository: 
[MVE GitHub repository link](https://github.com/zskb/molecular-vibration-explorer)

There, we also included a ***dft2app*** folder that contains: 
    * Python scripts used to convert Gaussian output files to the input files required by the molecular-vibration-explorer application 
    * Local notebooks that enable to quickly visualize DFT calculations’ results not included in the molecular-vibration-explorer database 
    * A few molecules’ DFT outputs to test the scripts

To use the application on one’s own DFT calculations, two scripts need to be run before visualizing the results via the local notebooks: 
*create_database_files.py*: creates files for the database analysis tool, only needed if one wants to compare with other molecules in the database or extend it. This generates the 2D pictures for the molecules
*create_molecular_data_files.py*: creates files used by the molecule and normal mode tools, such as dat files and 3D geometries


## Acknowledgements

This work received funding from the European Union’s Horizon 2020 research and innovation program under Grant Agreement No. 829067 (FET Open THOR). The authors want to thank Giovanni Pizzi and Dou Du for their technical support during the development of the website. 
