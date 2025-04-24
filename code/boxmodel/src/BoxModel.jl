module BoxModel
using RateLaws
using ModelingToolkit, DifferentialEquations, Catalyst, Unitful
export photolysisMCM, rxn_sys, sources, sinks, rxn_rate, dxdt

# MCM parameterization
#j(k) = l(k)*cosx**( mm(k))*exp(-nn(k)*secx) 
function photolysisMCM(l, mm, nn; zenith=0)
    return (l*cosd(zenith)^(mm)*exp(-1.0*nn*secd(zenith)))
end

################################################################
######### DEFINE MECHANISM AND BOX MODEL WITH CATALYST #########
################################################################

# time variable
@parameters t
D = Differential(t)

# number density -> ppb
@parameters c = 2.46e10
@parameters per_sec = 1.0

# uptake parameters
@parameters mean_molecular_speed = sqrt((8*1.38e-23*288.15)/(pi*5.48e-26))
@parameters γHO2 = 0.2
@parameters ASA = 200

##### reaction rate parameters
# gas phase reactions

# photolysis reactions
@parameters jO3_O1D = photolysisMCM(6.073E-05, 1.743, 0.474)
@parameters jO3_O3P = photolysisMCM(6.073E-05, 1.743, 0.474)
@parameters jH2O2 = photolysisMCM(1.041E-05, 0.723, 0.279)
@parameters jNO2 = photolysisMCM(1.165E-02, 0.244, 0.267)
@parameters jNO3_NO = photolysisMCM(2.485E-02, 0.168, 0.108)
@parameters jNO3_NO2 = photolysisMCM(1.747E-01, 0.155, 0.125)
@parameters jHONO = photolysisMCM(2.644E-03, 0.261, 0.288)
@parameters jHNO3 = photolysisMCM(9.312E-07, 1.23, 0.307)
@parameters jHCHO_H = photolysisMCM(4.642E-05, 0.762, 0.353)
@parameters jHCHO_H2 = photolysisMCM(6.853E-05, 0.477, 0.323)
@parameters jCH3OOH = photolysisMCM(7.649E-06, 0.682, 0.279)

##### species

#standard
@parameters CH2O        = 1  [isconstantspecies=true]  #CH2O; Formaldehyde
@parameters CH4         = 1800.0 [isconstantspecies=true]  #CH4; Methane
@parameters CO          = 1000 [isconstantspecies=true]  #CO; Carbon monoxide
@parameters CO2         = 355000.0  [isconstantspecies=true] #CO2; Carbon dioxide
@parameters H2O2        = 4e-6 [isconstantspecies=true]  #H2O2; Hydrogen peroxide
@parameters HCOOH       = 0.0  [isconstantspecies=true] #HCOOH; Formic acid
@parameters HNO2        = 0.0  [isconstantspecies=true] #HONO; Nitrous acid
@parameters HNO3        = 0.0  [isconstantspecies=true] #HNO3; Nitric acid
@parameters HNO4        = 2e-4  [isconstantspecies=true] #HNO4; Pernitric acid
@species HO2(t)         = 4e-6   #HO2; Hydroperoxyl radical
@parameters N2O         = 300.0  [isconstantspecies=true] #N2O; Nitrous oxide
@parameters N2O5        = 0.0  [isconstantspecies=true] #N2O5; Dinitrogen pentoxide
@parameters NO          = 0.0004  [isconstantspecies=true] #NO; Nitric oxide
@species NO2(t)         = 0.0004   #NO2; Nitrogen dioxide
@species NO3(t)         = 0.0   #NO3; Nitrate radical
@species O(t)           = 0.0   #O(3P); Ground state atomic oxygen
@species O1D(t)         = 1e-6   #O(1D); Excited atomic oxygen
@species O3(t)          = 40.0   #O3; Ozone
@species OH(t)          = 3e-4   #OH; Hydroxyl radical

#methane oxidation
@species MO2(t)         = 4e-6   #CH3O2; Methylperoxy radical
@parameters MO          = 4e-6  [isconstantspecies=true] #CH3O; Methoxy radical
@parameters MOH         = 0.0  [isconstantspecies=true] #CH3OH; Methanol
@parameters MP          = 4e-6  [isconstantspecies=true] #CH3OOH; Methylhydroperoxide

#atomic
@parameters H           = 2e-6  [isconstantspecies=true]  #H; Atomic hydrogen
@parameters N           = 4e-6  [isconstantspecies=true]  #N; Atomic nitrogen

#quasi-constant
@parameters H2O         = 18390000.0  [isconstantspecies=true] #H2O; Water vapor
@parameters H2          = 500.0 [isconstantspecies=true]  #H2; Molecular hydrogen
@parameters N2          = 780800000.0  [isconstantspecies=true] #N2; Molecular nitrogen
@parameters O2          = 209500000.0  [isconstantspecies=true] #O2; Molecular oxygen

### custom VOC scheme to simplify smoke chemistry
@parameters RH          = 8e-3  [isconstantspecies=true] #total VOC population
@species RO2(t)         = 0.0000   #lumped peroxy radical
@parameters RCHO        = 0.0000  [isconstantspecies=true] #
@parameters RCO3        = 0.0000  [isconstantspecies=true] #
@parameters RONO2       = 0.0000  [isconstantspecies=true] #
@parameters ROOH        = 0.0000  [isconstantspecies=true] #

### VOC speciation determines these parameters
@parameters α            = 0.05 # RO2 + NO branching ratio
@parameters kOH          = GCARR_ac(2.54e-11, 410.0e0) # oxidation rate
@parameters kO3          = 1.3e-17 # ozonolysis rate

# light intensity depends on time
@parameters I = 1.0

rxs = [
        ##### GAS-PHASE
        Reaction(GCARR_ac(6E-34, 2.4)*NUMDEN*c, [O, O2], [O3], [1, 1], [1])  #  O + O2 = O3
        Reaction(GCARR_ac(8.00e-12, -2060.0e0)*c, [O, O3], [O2], [1, 1], [2.000])  #  O + O3 = 2.000O2
        Reaction(GCARR_ac(2.15e-11, 110.0e0)*c, [O1D, N2], [O, N2], [1, 1], [1, 1])  #  O1D + N2 = O + N2   
        Reaction(GCARR_ac(3.30e-11, 55.0e0)*c, [O1D, O2], [O, O2], [1, 1], [1, 1])  #  O1D + O2 = O + O2
        Reaction(GCARR_ac(1.63e-10, 60.0e0)*c, [O1D, H2O], [OH], [1, 1], [2.000])  #  O1D + H2O = 2.000OH
        Reaction(1.20e-10*c, [O1D, H2], [H, OH], [1, 1], [1, 1])  #  O1D + H2 = H + OH   
        Reaction(GCARR_ac(2.80e-12, -1800.0e0)*c, [OH, H2], [H2O, H], [1, 1], [1, 1])  #  OH + H2 = H2O + H 
        Reaction(GCARR_ac(1.80e-11, 180.0e0)*c, [O, OH], [O2, H], [1, 1], [1, 1])  #  O + OH = O2 + H   
        Reaction(GCARR_ac(3.00e-11, 200.0e0)*c, [HO2, O], [OH, O2], [1, 1], [1, 1])  #  HO2 + O = OH + O2
        Reaction(GCARR_ac(1.70e-12, -940.0e0)*c, [O3, OH], [HO2, O2], [1, 1], [1, 1])  #  O3 + OH = HO2 + O2
        Reaction(GCARR_ac(1.00e-14, -490.0e0)*c, [O3, HO2], [OH, O2], [1, 1], [1, 2])  #  O3 + HO2 = OH + O2 + O2
        Reaction(GC_HO2HO2_acac(3.00e-13, 460.0e0, 2.1e-33, 920.0e0)*c, [HO2], [H2O2, O2], [2], [1, 1])  #  HO2 + HO2 = H2O2 + O2
        Reaction(1.80e-12*c, [OH, H2O2], [H2O, HO2], [1, 1], [1, 1])  #  OH + H2O2 = H2O + HO2
        Reaction(GCARR_ac(4.80e-11, 250.0e0)*c, [OH, HO2], [H2O, O2], [1, 1], [1, 1])  #  OH + HO2 = H2O + O2
        Reaction(1.80e-12*c, [OH], [H2O, O], [2], [1, 1])  #  OH + OH = H2O + O
        Reaction(GCJPLPR_aba(6.90e-31, 1.0e+00, 2.6e-11, 0.6e0)*c, [OH], [H2O2], [2], [1])  #  OH + OH = H2O2
        Reaction(GCARR_ac(4.63e-11, 20.0e0)*c, [O1D, N2O], [N2, O2], [1, 1], [1, 1])  #  O1D + N2O = N2 + O2   
        Reaction(GCARR_ac(7.25e-11, 20.0e0)*c, [O1D, N2O], [NO], [1, 1], [2.000])  #  O1D + N2O = 2.000NO
        Reaction(GCARR_ac(3.30e-12, 270.0e0)*c, [HO2, NO], [OH, NO2], [1, 1], [1, 1])  #  HO2 + NO = OH + NO2
        Reaction(GCARR_ac(3.00e-12, -1500.0e0)*c, [O3, NO], [NO2, O2], [1, 1], [1, 1])  #  O3 + NO = NO2 + O2
        Reaction(GCARR_ac(5.10e-12, 210.0e0)*c, [NO2, O], [NO, O2], [1, 1], [1, 1])  #  NO2 + O = NO + O2
        Reaction(GCARR_ac(1.20e-13, -2450.0e0)*c, [O3, NO2], [O2, NO3], [1, 1], [1, 1])  #  O3 + NO2 = O2 + NO3
        Reaction(3.50e-12*c, [HO2, NO3], [OH, NO2, O2], [1, 1], [1, 1, 1])  #  HO2 + NO3 = OH + NO2 + O2

        Reaction(GCJPLPR_abab(2.40e-30, 3.0e+00, 1.6e-12, -0.1e0, 0.6e0)*c, [NO2, NO3], [N2O5], [1, 1], [1])  #  NO2 + NO3 = N2O5
        Reaction(GCARR_ac(4.35e-14, -1335.0e0)*c, [NO2, NO3], [NO, NO2, O2], [1, 1], [1, 1, 1])  #  NO2 + NO3 = NO + NO2 + O2
        Reaction(GCJPLPR_abcabc(4.14e-04, 3.0e0, -10840.0e0, 2.76e14, -0.1e0, -10840.0e0, 0.6e0)*per_sec, [N2O5], [NO2, NO3], [1], [1, 1])  #  N2O5 = NO2 + NO3
        Reaction(GCJPLPR_aba(1.80e-30, 3.0e+00, 2.8e-11, 0.6e0)*c, [NO2, OH], [HNO3], [1, 1], [1])  #  NO2 + OH = HNO3
        Reaction(GC_OHHNO3_acacac(2.41e-14, 460.0e0, 2.69e-17, 2199.0e0, 6.51e-34, 1335.0e0)*c, [HNO3, OH], [H2O, NO3], [1, 1], [1, 1])  #  HNO3 + OH = H2O + NO3

        ## added HONO chem myself
        Reaction(GCJPLPR_abab(7.00e-31, 2.6e+00, 3.60e-11, 0.1e0, 0.6e0)*c, [NO, OH], [HNO2], [1, 1], [1])  #  NO + OH = HNO2
        Reaction(GCARR_ac(1.80e-11, -390.0e0)*c, [HNO2, OH], [H2O, NO2], [1, 1], [1, 1])  #  HNO2 + OH = H2O + NO2

        Reaction(GCARR_ac(1.70e-11, 125.0e0)*c, [NO, NO3], [NO2], [1, 1], [2.000])  #  NO + NO3 = 2.000NO2
        Reaction(GCJPLPR_abab(1.90e-31, 3.4e+00, 4.0e-12, 0.3e0, 0.6e0)*c, [HO2, NO2], [HNO4], [1, 1], [1])  #  HO2 + NO2 = HNO4
        Reaction(GCJPLPR_abcabc(9.05e-05, 3.4e0, -10900.0e0, 1.90e15, 0.3e0, -10900.0e0, 0.6e0)*per_sec, [HNO4], [HO2, NO2], [1], [1, 1])  #  HNO4 = HO2 + NO2
        Reaction(GCARR_ac(1.30e-12, 380.0e0)*c, [HNO4, OH], [H2O, NO2, O2], [1, 1], [1, 1, 1])  #  HNO4 + OH = H2O + NO2 + O2

        Reaction(GCARR_ac(2.45e-12, -1775.0e0)*c, [OH, CH4], [MO2, H2O], [1, 1], [1, 1])  #  OH + CH4 = MO2 + H2O
        Reaction(1.31e-10*c, [O1D, CH4], [MO2, OH], [1, 1], [1, 1])  #  O1D + CH4 = MO2 + OH   
        Reaction(0.09e-10*c, [O1D, CH4], [CH2O, H2], [1, 1], [1, 1])  #  O1D + CH4 = CH2O + H2   
        Reaction(0.35e-10*c, [O1D, CH4], [CH2O, H, HO2], [1, 1], [1, 1, 1])  #  O1D + CH4 = CH2O + H + HO2 

        Reaction(GC_RO2NO_B1_ac(2.80e-12, 300.0e0)*c, [MO2, NO], [CH2O, HO2, NO2], [1, 1], [1, 1, 1])  #  MO2 + NO = CH2O + HO2 + NO2
        Reaction(GCARR_abc(4.10e-13, 0.0e0, 750.0e0)*c, [MO2, HO2], [MP, O2], [1, 1], [1, 1])  #  MO2 + HO2 = MP + O2
        Reaction(GC_TBRANCH_1_acac(9.50e-14, 390.0e0, 2.62e1, -1130.0e0)*c, [MO2], [MOH, CH2O, O2], [2], [1, 1, 1])  #  MO2 + MO2 = MOH + CH2O + O2
        Reaction(GC_TBRANCH_1_acac(9.50e-14, 390.0e0, 4.0e-2, 1130.0e0)*c, [MO2], [CH2O, HO2], [2], [2.000, 2.000])  #  MO2 + MO2 = 2.000CH2O + 2.000HO2
        Reaction(1.60e-10*c, [MO2, OH], [MOH, CH2O, HO2], [1, 1], [0.13, 0.87, 1.74])  #  MO2 + OH = 0.13MOH + 0.87CH2O + 1.74HO2
        Reaction(GCARR_ac(2.66e-12, 200.0e0)*c, [MP, OH], [MO2, H2O], [1, 1], [1, 1])  #  MP + OH = MO2 + H2O   
        Reaction(GCARR_ac(3.80e-12, 200.0e0)*c, [MP, OH], [MO2, CH2O, OH, H2O], [1, 1], [0.700, 0.3, 0.3, 1])  #  MP + OH = MO2 + CH2O + OH + H2O

        Reaction(GCARR_ac(5.50e-12, 125.0e0)*c, [CH2O, OH], [CO, HO2, H2O], [1, 1], [1, 1, 1])  #  CH2O + OH = CO + HO2 + H2O
        Reaction(5.80e-16*c, [NO3, CH2O], [HNO3, HO2, CO], [1, 1], [1, 1, 1])  #  NO3 + CH2O = HNO3 + HO2 + CO

        Reaction(GC_OHCO_a(1.50e-13)*c, [OH, CO], [HO2, CO2], [1, 1], [1, 1])  #  OH + CO = HO2 + CO2
        Reaction(GCARR_ac(2.90e-12, -345.0e0)*c, [MOH, OH], [HO2, CH2O], [1, 1], [1, 1])  #  MOH + OH = HO2 + CH2O

        ##### VOC
        Reaction(kOH*c, [RH, OH], [RO2], [1, 1], [1])  #  RH + OH = RO2
        Reaction(GCARR_ac(2.80e-12, -3280.0e0)*c, [RH, NO3], [RO2, HNO3], [1, 1], [1, 1])  #  RH + NO3 = RO2 + HNO3
        Reaction(GC_RO2NO_B2_aca(2.90e-12, 350.0e0, 3.0e0)*c, [RO2, NO], [NO2, RCHO, HO2, RONO2], [1, 1], [(1-α), (1-α), (1-α), α])  #  RO2 + NO = NO2 + RCHO + HO2(1-α) + RONO2(α)
        Reaction(GCARR_ac(7.40e-13, 700.0e0)*c, [RO2, HO2], [ROOH], [1, 1], [1])  #  RO2 + HO2 = ROOH
        Reaction(GCARR_ac(6.00e-12, 410.0e0)*c, [RCHO, OH], [RCO3, H2O], [1, 1], [1, 1])  #  RCHO + OH = RCO3 + H2O
        Reaction(6.50e-15*c, [RCHO, NO3], [HNO3, RCO3], [1, 1], [1, 1])  #  RCHO + NO3 = HNO3 + RCO3   
        Reaction(kO3*c, [RH, O3], [RO2, OH, CO2, CO, HO2, CH2O, H2O2], [1,1], [0.2, 0.28, 0.407, 0.407, 0.16, 0.827, 0.013])  #  RH + O3 = RO2 + .... 

        ##### PHOTOLYSIS
        Reaction(jO3_O1D*min.(1,I), [O3], [O1D, O2], [1], [1, 1])  #  O3 -> O2 + O(1D)
        Reaction(jO3_O3P*min.(1,I), [O3], [O, O2], [1], [1, 1])  #  O3 -> O2 + O(3P)
        Reaction(jH2O2*min.(1,I), [H2O2], [OH], [1], [2])  #  H2O2 -> 2OH
        Reaction(jNO2*min.(1,I), [NO2], [NO, O], [1], [1, 1])  #  NO2 -> NO + O(3P)
        Reaction(jNO3_NO*min.(1,I), [NO3], [NO, O2], [1], [1, 1])  #  NO3 -> NO + O2
        Reaction(jNO3_NO2*min.(1,I), [NO3], [NO2, O], [1], [1, 1])  #  NO3 -> NO2 + O(3P)
        Reaction(jHONO*min.(1,I), [HNO2], [NO, OH], [1], [1, 1])  #  HNO2 -> NO + OH
        Reaction(jHNO3*min.(1,I), [HNO3], [NO2, OH], [1], [1, 1])  #  HNO3 -> NO2 + OH
        Reaction(jHCHO_H*min.(1,I), [CH2O], [CO, HO2], [1], [1, 2])  #  HCHO -> CO + 2HO2
        Reaction(jHCHO_H2*min.(1,I), [CH2O], [CO, H2], [1], [1, 1])  #  HCHO -> CO + H2
        Reaction(jCH3OOH*min.(1,I), [MP], [MO, OH], [1], [1, 1])  #  CH3OOH -> CH3O + OH
        Reaction(jHCHO_H*min.(1,I), [RCHO], [RO2, CO, HO2], [1], [1, 1, 1])  #  RCHO -> RO2 + CO + HO2

        ##### DEPOSITION
        # Reaction((1/56/3600)*per_sec, [O3], nothing, [1], nothing)

        ##### aerosol uptake
        Reaction(0.25*mean_molecular_speed*γHO2*(ASA*(1e-6)^2*(0.01)^(-3))*per_sec, [HO2], nothing, [1], nothing)  #  HO2 uptake
];

const rxn_sys = ReactionSystem(rxs, t; name=:smokechem);

# Mark the reaction system as complete
rxn_sys = complete(rxn_sys);

(convert(ODESystem, rxn_sys; combinatoric_ratelaws=false), rxn_sys);


################################################
############# QUERY REACTIONSYSTEM #############
################################################ 

function sources(compound)
    # Find reactions involving the target species as a product
    prod_rxs = filter(reaction -> any(substrate -> isequal(substrate, compound), reaction.products), reactions(rxn_sys));
    
    return prod_rxs
end

function sinks(compound)
    # Find reactions involving the target species as a substrate
    loss_rxs = filter(reaction -> any(substrate -> isequal(substrate, compound), reaction.substrates), reactions(rxn_sys))
    
    return loss_rxs
end
function rxn_rate(rxn, rxn_sys, u, p)
    # symbolic rate law for a given Reaction
    rl = oderatelaw(rxn, combinatoric_ratelaw=false)

    # Julia function that evaluates the ratelaw's value
    f =  build_function(rl, states(rxn_sys), parameters(rxn_sys), independent_variable(rxn_sys), expression=Val{false})

    # evaluating for a specific u and t

    t = 1.0
    rl_value = f(u, p, t)
    
    
    return rl_value
end
function dxdt(rxns, rxn_sys, u, p)
    total_rate = sum([rxn_rate(rxn, rxn_sys, u, p) for rxn in rxns]);

    return (total_rate)
end

end