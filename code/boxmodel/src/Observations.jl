module Observations
using RateLaws
using Statistics, Random, Distributions, DataFrames, CSV
export near_field, far_field, dVOC, dNOxCO, dCH2O, dHNO2, dASA, dkOH, dkO3, dalpha, djNO2, dCO

### load the observations
df = DataFrame(CSV.File("/home/jpalmo/fs09/Projects/O3Fire/model_eval/output/model_obs_2024.7.31.csv", dateformat="yyyy-mm-dd"));
### split by regime
fire_df = df[df.regime .== "fire", :];
mixed_df = df[df.regime .== "heavy_mixed", :];

### derive ASA from mass concentration assuming a constant 
# specific surface area of 4 m2 g−1 (following Lindsay et al 2022)
fire_df.ASA = fire_df.OA_AMS .* 4;
### Filter out rows with missing values in the "age" column
fire_df = dropmissing(fire_df, :age);
fire_df = dropmissing(fire_df, :NOx);
### scale concentrations by CO
fire_df.NOxCO = fire_df.NOx ./ fire_df.CO;
fire_df.CH2OCO = fire_df.CH2O ./ fire_df.CO;
fire_df.HNO2CO = fire_df.HNO2 ./ fire_df.CO;

### filter non-physical values
fire_df = fire_df[fire_df.NOx .> 0, :];

### split by age
near_field = fire_df[fire_df.age .< 20, :];
far_field = fire_df[fire_df.age .> 20, :];

### drop missing
near_field = dropmissing(near_field, :NOxCO);
near_field = disallowmissing!(near_field, :NOxCO);

### Fit Distributions to observations for bootstrapping
# Taken both from my dataset and Gkatzelis et al. 2024

### NOx/CO
#from Gkatzelis
dNOxCO = truncated(Normal{Float64}(5.37e-3, 4.92e-3), lower=0)

### VOC
#from Gkatzelis
dVOC = truncated(Normal{Float64}(0.13424, 0.01823), lower=0);

### CH2O
#from Gkatzelis
dCH2O = truncated(Normal{Float64}(17.92e-3, 4.31e-3), lower=0)

### HONO
#from Gkatzelis
dHNO2 = truncated(Normal{Float64}(1.89e-3, 1.61e-3), lower=0)

### jNO2
# from near-field photolysis observations
near_field = dropmissing(near_field, :jNO2);
near_field = disallowmissing!(near_field, :jNO2);
djNO2 = truncated(fit(Normal{Float64}, near_field.jNO2), lower=0)

### CO
# from near-field CO observations
near_field = dropmissing(near_field, :CO);
near_field = disallowmissing!(near_field, :CO);
dCO = truncated(fit(Exponential{Float64}, near_field.CO), lower=0)

### ASA
# ASA_df = DataFrame(CSV.File("ASA.csv",));
# ASA_df = dropmissing(ASA_df, :sSMPS_stdPT_MOORE);
# ASA_df = disallowmissing!(ASA_df, :sSMPS_stdPT_MOORE);
# this is made up
dASA = truncated(Normal{Float64}(317.3, 148.9), lower=0);


end
