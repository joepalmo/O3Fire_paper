push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using RateLaws
using BoxModel
using Observations
using ModelingToolkit, DifferentialEquations, Catalyst, Unitful, MonteCarloMeasurements, LaTeXStrings, Distributions, Bootstrap, Random, CSV, DataFrames, Statistics

### Bootstrap values from the distributions defined in the Observations module
# these represent the "at emission" scenario
# most are drawn from emission ratios reported by Gkatzelis et al. 2024,
# except for ASA, which is arbitrarily defined as a normal distribution with a mean that is 
bootstrap_N = 500

VOC_bootstrap = rand(dVOC, bootstrap_N,)
NOxCO_bootstrap = rand(dNOxCO, bootstrap_N,)
CH2O_bootstrap = rand(dCH2O, bootstrap_N,)
HNO2_bootstrap = rand(dHNO2, bootstrap_N,);
ASA_bootstrap = rand(dASA, bootstrap_N,).*4;

kHO2NO = GCARR_ac(3.30e-12, 270.0e0)*2.46e10;
kNORO2 = GC_RO2NO_B2_aca(2.90e-12, 350.0e0, 3.0e0)*2.46e10;

tspan = (0, 24*60*60); # 24 hours, this doesn't matter since we are running to steady-state

# number of bootstrapped samples per sensitivity test
# runtime will largely depend on this number
N = 100

println("Starting sensitivity testing (N = $N)")

#-------------------------------------#
#######################################
####### RUN SENSITIVITY TESTING #######
#######################################
#-------------------------------------#

#-------------------------------------#
############# NEAR FIELD ##############
#-------------------------------------#

# define parameter and concentration priors to bootstrap from:
#### observations
CO_bootstrap = sample(collect(skipmissing(near_field.CO)), N, replace=true)./1000;
jNO2_bootstrap = sample(collect(skipmissing(near_field.jNO2)), N, replace=true)./photolysisMCM(1.165E-02, 0.244, 0.267, zenith=0);
ASA_bootstrap = sample(collect(skipmissing(near_field.ASA)), N, replace=true);
#### Gkatzelis
VOCCO_bootstrap = rand(dVOC, N,)
NOxCO_bootstrap =  sample(collect(skipmissing(near_field.NOx./near_field.CO)), N, replace=true);
CH2OCO_bootstrap = sample(collect(skipmissing(near_field.CH2O./near_field.CO)), N, replace=true);
HNO2CO_bootstrap = sample(collect(skipmissing(near_field.HNO2./near_field.CO)), N, replace=true);

# define default parameter priors
mean_co = (mean(skipmissing(near_field.CO))./1000)

default_parameters = Dict("CO" => fill(mean_co, N), "VOCCO" => fill(mean(dVOC), N), "NOxCO" => fill(mean(skipmissing(near_field.NOx./near_field.CO)), N), "CH2OCO" => fill(mean(skipmissing(near_field.CH2O./near_field.CO)), N), "HNO2CO" => fill(mean(skipmissing(near_field.HNO2./near_field.CO)), N), "I" => fill(1, N), "ASA" => fill(0, N), "alpha" => fill(0.05, N), "kOH" => fill(GCARR_ac(2.54e-11, 410.0e0), N), "kO3" => fill(1.3e-17, N),);
#perturbed parameters
perturbed_parameters = Dict("default" => [], "CO" => CO_bootstrap, "VOCCO" => VOCCO_bootstrap, "NOxCO" => NOxCO_bootstrap, "CH2OCO" => CH2OCO_bootstrap, "HNO2CO" => HNO2CO_bootstrap, "I" => jNO2_bootstrap, "ASA" => ASA_bootstrap,);

# analyze sensitivity of the model to the parameters one-by-one
sensitivity_results = Dict()
simprob_results = Dict()

# Define the ODE problem
oprob = ODEProblem(rxn_sys, [], tspan, [], combinatoric_ratelaws=false);
sprob = SteadyStateProblem(oprob);

for x in perturbed_parameters
    # Use a loop to update the values in place
    updated_parameters = copy(default_parameters);
    for (key, new_value) in Dict(x)
        if haskey(updated_parameters, key)
            updated_parameters[key] = new_value
        end
    end

    # Define function for building EnsembleProblem
    simprobs = []
    function prob_func(prob, i, repeat)
        u_tmp = [:NO => updated_parameters["NOxCO"][i].*updated_parameters["CO"][i], :RH => updated_parameters["VOCCO"][i].*updated_parameters["CO"][i], :CO => updated_parameters["CO"][i], :CH2O => updated_parameters["CH2OCO"][i].*updated_parameters["CO"][i], :HNO2 => updated_parameters["HNO2CO"][i].*updated_parameters["CO"][i],];
        p_tmp = [:I => updated_parameters["I"][i], :ASA => updated_parameters["ASA"][i], :α => updated_parameters["alpha"][i], :kOH => updated_parameters["kOH"][i], :kO3 => updated_parameters["kO3"][i]];
        oprob = ODEProblem(rxn_sys, u_tmp, tspan, p_tmp, combinatoric_ratelaws=false);
        push!(simprobs, SteadyStateProblem(oprob))
        SteadyStateProblem(oprob);
    end;
    ensemble_prob = EnsembleProblem(sprob, prob_func = prob_func);
    sim = solve(ensemble_prob, DynamicSS(Rosenbrock23()), EnsembleThreads(), trajectories = N);
    sensitivity_results[x[1]] = sim
    simprob_results[x[1]] = simprobs
end;

sensitivity_results_near_field = deepcopy(sensitivity_results);
simprob_results_near_field = deepcopy(simprob_results);

println("Finished near field sensitivity testing")

#-------------------------------------#
############# FAR FIELD ###############
#-------------------------------------#

# define parameter and concentration priors to bootstrap from:
#### observations
CO_bootstrap = sample(collect(skipmissing(far_field.CO)), N, replace=true)./1000;
jNO2_bootstrap = sample(collect(skipmissing(far_field.jNO2)), N, replace=true)./photolysisMCM(1.165E-02, 0.244, 0.267, zenith=0);
ASA_bootstrap = sample(collect(skipmissing(far_field.ASA)), N, replace=true);
#### Gkatzelis
VOCCO_bootstrap = rand(dVOC, N,)
NOxCO_bootstrap =  sample(collect(skipmissing(far_field.NOx./far_field.CO)), N, replace=true);
CH2OCO_bootstrap = sample(collect(skipmissing(far_field.CH2O./far_field.CO)), N, replace=true);
HNO2CO_bootstrap = sample(collect(skipmissing(far_field.HNO2./far_field.CO)), N, replace=true);

# define default parameter priors
mean_co = (mean(skipmissing(far_field.CO))./1000)

default_parameters = Dict("CO" => fill(mean_co, N), "VOCCO" => fill(mean(dVOC), N), "NOxCO" => fill(mean(skipmissing(far_field.NOx./far_field.CO)), N), "CH2OCO" => fill(mean(skipmissing(far_field.CH2O./far_field.CO)), N), "HNO2CO" => fill(mean(skipmissing(far_field.HNO2./far_field.CO)), N), "I" => fill(1, N), "ASA" => fill(0, N), "alpha" => fill(0.05, N), "kOH" => fill(GCARR_ac(2.54e-11, 410.0e0), N), "kO3" => fill(1.3e-17, N),);
#perturbed parameters
perturbed_parameters = Dict("default" => [], "CO" => CO_bootstrap, "VOCCO" => VOCCO_bootstrap, "NOxCO" => NOxCO_bootstrap, "CH2OCO" => CH2OCO_bootstrap, "HNO2CO" => HNO2CO_bootstrap, "I" => jNO2_bootstrap, "ASA" => ASA_bootstrap,);

# analyze sensitivity of the model to the parameters one-by-one
sensitivity_results = Dict()
simprob_results = Dict()

# Define the ODE problem
oprob = ODEProblem(rxn_sys, [], tspan, [], combinatoric_ratelaws=false);
sprob = SteadyStateProblem(oprob);

for x in perturbed_parameters
    # Use a loop to update the values in place
    updated_parameters = copy(default_parameters);
    for (key, new_value) in Dict(x)
        if haskey(updated_parameters, key)
            updated_parameters[key] = new_value
        end
    end

    # Define function for building EnsembleProblem
    simprobs = []
    function prob_func(prob, i, repeat)
        u_tmp = [:NO => updated_parameters["NOxCO"][i].*updated_parameters["CO"][i], :RH => updated_parameters["VOCCO"][i].*updated_parameters["CO"][i], :CO => updated_parameters["CO"][i], :CH2O => updated_parameters["CH2OCO"][i].*updated_parameters["CO"][i], :HNO2 => updated_parameters["HNO2CO"][i].*updated_parameters["CO"][i],];
        p_tmp = [:I => updated_parameters["I"][i], :ASA => updated_parameters["ASA"][i], :α => updated_parameters["alpha"][i], :kOH => updated_parameters["kOH"][i], :kO3 => updated_parameters["kO3"][i]];
        oprob = ODEProblem(rxn_sys, u_tmp, tspan, p_tmp, combinatoric_ratelaws=false);
        push!(simprobs, SteadyStateProblem(oprob))
        SteadyStateProblem(oprob);
    end;
    ensemble_prob = EnsembleProblem(sprob, prob_func = prob_func);
    sim = solve(ensemble_prob, DynamicSS(Rosenbrock23()), EnsembleThreads(), trajectories = N);
    sensitivity_results[x[1]] = sim
    simprob_results[x[1]] = simprobs
end;

sensitivity_results_far_field = deepcopy(sensitivity_results);
simprob_results_far_field = deepcopy(simprob_results);

println("Finished far field sensitivity testing")

#-------------------------------------#
############ AT EMISSION ##############
#-------------------------------------#

# define parameter and concentration priors to bootstrap from:
#### observations
CO_bootstrap = fill(4000, N)
# CO_bootstrap = sample(collect(skipmissing(far_field.CO)), N, replace=true)./1000;
jNO2_bootstrap = sample(collect(skipmissing(near_field.jNO2)), N, replace=true)./photolysisMCM(1.165E-02, 0.244, 0.267, zenith=0);
ASA_bootstrap = sample(collect(skipmissing(near_field.ASA)), N, replace=true);
#### Gkatzelis
VOCCO_bootstrap = rand(dVOC, N,);
NOxCO_bootstrap =  rand(dNOxCO, N,);
CH2OCO_bootstrap = rand(dCH2O, N,);
HNO2CO_bootstrap = rand(dHNO2, N,);

# define default parameter priors
mean_co = 4000

default_parameters = Dict("CO" => fill(mean_co, N), "VOCCO" => fill(mean(dVOC), N), "NOxCO" => fill(mean(dNOxCO), N), "CH2OCO" => fill(mean(dCH2O), N), "HNO2CO" => fill(mean(dHNO2), N), "I" => fill(1, N), "ASA" => fill(0, N),);
#perturbed parameters
perturbed_parameters = Dict("default" => [], "CO" => CO_bootstrap, "VOCCO" => VOCCO_bootstrap, "NOxCO" => NOxCO_bootstrap, "CH2OCO" => CH2OCO_bootstrap, "HNO2CO" => HNO2CO_bootstrap, "I" => jNO2_bootstrap, "ASA" => ASA_bootstrap,);

# analyze sensitivity of the model to the parameters one-by-one
sensitivity_results = Dict()
simprob_results = Dict()

# Define the ODE problem
oprob = ODEProblem(rxn_sys, [], tspan, [], combinatoric_ratelaws=false);
sprob = SteadyStateProblem(oprob);

for x in perturbed_parameters
    # Use a loop to update the values in place
    updated_parameters = copy(default_parameters);
    for (key, new_value) in Dict(x)
        if haskey(updated_parameters, key)
            updated_parameters[key] = new_value
        end
    end

    # Define function for building EnsembleProblem
    simprobs = []
    function prob_func(prob, i, repeat)
        u_tmp = [:NO => updated_parameters["NOxCO"][i].*updated_parameters["CO"][i], :RH => updated_parameters["VOCCO"][i].*updated_parameters["CO"][i], :CO => updated_parameters["CO"][i], :CH2O => updated_parameters["CH2OCO"][i].*updated_parameters["CO"][i], :HNO2 => updated_parameters["HNO2CO"][i].*updated_parameters["CO"][i],];
        p_tmp = [:I => updated_parameters["I"][i], :ASA => updated_parameters["ASA"][i],];
        oprob = ODEProblem(rxn_sys, u_tmp, tspan, p_tmp, combinatoric_ratelaws=false);
        push!(simprobs, SteadyStateProblem(oprob))
        SteadyStateProblem(oprob);
    end;
    ensemble_prob = EnsembleProblem(sprob, prob_func = prob_func);
    sim = solve(ensemble_prob, DynamicSS(Rosenbrock23()), EnsembleThreads(), trajectories = N);
    sensitivity_results[x[1]] = sim
    simprob_results[x[1]] = simprobs
end;

sensitivity_results_toe = deepcopy(sensitivity_results);
simprob_results_toe = deepcopy(simprob_results);

println("Finished at emission sensitivity testing")

#-------------------------------------#
########### COMBINE RESULTS ###########
#-------------------------------------#

# Combine the sensitivity results into a DataFrame
cases = ["near_field", "far_field", "toe"] 
factors = keys(Dict(sensitivity_results_near_field))

dataframe = DataFrame()
for key in factors
    for (j, case) in enumerate(cases)
        result = Dict(eval(Symbol("sensitivity_results_$case")))[key]
        values = zeros(Float64, N)
        for (i,s) in enumerate(result)
            tmp_prob = Dict(eval(Symbol("simprob_results_$case")))[key][i]
            values[i] = kHO2NO.*s[:HO2].*(ModelingToolkit.getp(tmp_prob, :NO)(tmp_prob)) + kNORO2.*s[:RO2].*(ModelingToolkit.getp(tmp_prob, :NO)(tmp_prob)).*(1 - (ModelingToolkit.getp(tmp_prob, :α)(tmp_prob)))*3600
        end
        col_name = Symbol("$key - $case")
        dataframe[!, col_name] = values;
    end
end

# Melt the DataFrame for plotting
melted_data = stack(dataframe)
rename!(melted_data, Dict(:variable => :Key, :value => :Sensitivity));
# Split the Key column into separate Factor and Case columns
split_keys = split.(melted_data.Key, " - ")
melted_data[!, :Factor] = first.(split_keys)
melted_data[!, :Case] = last.(split_keys)


CSV.write("sensitivity_results.csv", melted_data)

println("Results saved to sensitivity_results.csv")