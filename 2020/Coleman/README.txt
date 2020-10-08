% README file

To run the model:
Open the run file in Matlab, select values for the variable parameter inputs
Run the file (this will call the model file automatically)


Below details and explains the section of the model file

Global definitions contains:
- Diffusion coefficient equations
- Partition coefficient equations
- Compound data - information related to the chemicals being modelled
- Other parameters - information related to the model set-up, and parameter input for the coefficient equations above
- Correction factors/Constants
- Partition coefficient ratio input for boundaries - ratios required for input of partition coefficients

Aromatic amine data contains:
- All the required input parameters for the 6 aromatic amines detailed in Manwaring et al. 2015
- Stratum corneum diffusion coefficients for the 6 aromatic amines calculated using the Chen et al. model

Unknown parameters for analysis contains:
- The parameters that have been commented to be adjusted/fitted
- This links these parameters to values input in the run file

Geometry contains:
- Details the domains and their lengths (1D)

Definitions contains:
- Integrations of the results for analysis

Transport of dilute species contains:
- Defining the location of the physics and concentrations wrt. the domains created in geometry
- Details  the related diffusion and partition, and metabolism where needed, wrt. the domains created in geometry

Mesh contains:
- Details of the mesh sizes in each domain

Study contains:
- The creation and definition of the study

Create plots contains
- The plots the define what data to output related to which domain and which format

Transient contains: 
- The study study length and timestep

Export contains: 
- Exports the plots defined above