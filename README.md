Description of folders:
'data' contains all model simulation output used in the figures and supplemental material for Le Roy et al. 2024 "Impact of climate variability and change on the ozone response to NOx emission reductions"
'code' contains all Jupyter notebooks used to generate the figures and supplemental materials for Palmo et al. 2025
'ozclim' contains the project package which contains useful functions

Description of files:
'setup.py' contains instructions for installing the ozclim package
'environment.yml' contains the project-specific Python virtual environment ("OZCLIM-env")

File naming convention in 'data':
GCHP.<realization>_<scenario>_c48.<variable_type>.nc4

<realization> is the shortname for each of the 5 initial condition realizations: w10, w13, w14, w26, w28
<scenario>* is one of the three climate scenarios: ref (present), refshort (future-high), pol3.7short (future-low)
<variable_type>** refers to one of either: Emissions, SpeciesConc, or MDA8_O3

*the designation 'SNOx' ('LNOx' in the manuscript) means that anthropogenic NOx emissions were reduced by 10%: refSNOx (present-LNOx), refSNOxshort (future-high-LNOx), pol3.7SNOxshort (future-low-LNOx)

**Emissions contains soil NO emissions (EmisNO_Soil), total NO emissions (EmisNO_Total), biogenic isoprene (EmisISOP_Biogenic), and total isoprene emissions (EmisISOP_Total);
SpeciesConc contains mixing ratios of: H2O (SpeciesConc_H2O), HO2 (SpeciesConc_HO2), OH (SpeciesConc_OH), PAN (SpeciesConc_PAN);
MDA8_O3 contains maximum daily 8-hour average ozone

The 'tools' directory within 'data' contains useful files for the ozclim.tools package (cube-sphere landmask).
