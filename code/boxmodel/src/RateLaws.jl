module RateLaws
export NUMDEN, GCARR_ab, GCARR_ac, GCARR_abc, GC_HO2HO2_acac, GC_TBRANCH_1_acac, GC_GLYXNO3_ac, GC_GLYCOH_A_a, GC_GLYCOH_B_a, GC_OHCO_a, GC_RO2NO_B1_ac, GCJPLPR_aba, GCJPLPR_abab, GCJPLPR_abcabc, GC_RO2NO_A1_ac, GC_RO2NO_A2_aca, GC_RO2NO_B2_aca, GC_OHHNO3_acacac, GCJPLEQ_acabab


##########################################################################
# these rate laws are taken from GEOS-Chem and translated from FORTRAN
##########################################################################


temp = 288.15
INV_TEMP = 1.0 / temp
temp_over_K300 = temp / 300.0
k300_over_temp = 300.0 / temp
SR_TEMP = sqrt(temp)
NUMDEN = 2.7e19 #1
NSPEC = 1
NREACT = 1
H2O = 1


##########################################################################
######                   ARRHENIUS FUNCTIONS                         #####
##########################################################################

function GCARR_ab(a0, b0)
  # Arrhenius function, skipping computation of EXP( c0 / T ),
  # which evaluates to 1 when c0 = 0.0  This avoids excess CPU
  # cycles. (bmy, 12 / 18 / 20)

  k = a0 * k300_over_temp^b0
  return k 
end

function GCARR_ac(a0, c0)
  # Arrhenius function, skipping computation of ( 300 / T )^b0,
  # which evaluates to 1 when b0 = 0.0  This avoids excess CPU
  # cycles (bmy, 12 / 18 / 20)
  #


  k = a0 * exp(c0 / temp)
  return k
end

function GCARR_abc(a0, b0, c0)
  # Arrhenius function, using all 3 terms.
  # Use this when a0, b0, c0 are all nonzero.
  #

  k = a0 * exp(c0 / temp) * k300_over_temp^b0
  return k
end



function GC_HO2HO2_acac(a0, c0, a1, c1)
    # Used to compute the rate for these reactions:
    #    HO2 + HO2 = H2O2 + O2
    #
    # For this reaction, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1
    # because b0 = b1 = 0.  Therefore we can skip computing these
    # terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0, a1, c1
    # REAL(dp)             :: k0, k1, k
    #
    k0 = a0 * exp(c0 / temp)
    k1 = a1 * exp(c1 / temp)
    k = (k0 + k1 * NUMDEN) * (1.0 + 1.4E-21 * exp(2200.0 / temp))
    return k
end

function GC_TBRANCH_1_acac(a0, c0, a1, c1)
    # temperature Dependent Branching Ratio, used for reactions:
    #    MO2 + MO2 = CH2O  + MOH + O2
    #    MO2 + MO2 = 2CH2O + 2HO2
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1
    # because b0 = b1 = 0.  Therefore we can skip computing these
    # terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0, a1, c1
    # REAL(dp)             :: k0, k1, k
    #
    k0 = a0 * exp(c0 / temp)
    k1 = a1 * exp(c1 / temp)
    k = k0 / (1.0 + k1)
    return k
end


function GC_GLYXNO3_ac(a0, c0)
    # Reaction rate for:
    #    GLYX + NO3 = HNO3 + HO2 + 2CO
    #    i.e. the HO2 + 2*CO branch
    #
    # For this reaction, this Arrhenius term evaluates to 1:
    #    (300/T)^b0
    # because b0 = 0.  Therefore we can skip computing this
    # term.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0
    # REAL(dp)             :: O2, k
    #
    # ---  K = k1*([O2]+3.5D18)/(2*[O2]+3.5D18)
    O2 = NUMDEN * 0.2095
    k = a0 * exp(c0 / temp)
    k = k * (O2 + 3.5e+18) / (2.0 * O2 + 3.5e+18)
    return k
end


function GC_GLYCOH_A_a(a0)
    # Used to compute the rate for this reaction:
    #    GLYC + OH = 0.732CH2O + 0.361CO2  + 0.505CO    + 0.227OH
    #              + 0.773HO2  + 0.134GLYX + 0.134HCOOH
    # which is the "A" branch of GLYC + OH.
    #
    # For this reaction, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0 * exp(c0/T)
    # Because b0 = c0 = 0.  Therefore we can skip computing these
    # terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0
    # REAL(dp)             :: glyc_frac, k
    # REAL(dp), PARAMETER  :: 
    exp_arg = -1.0 / 73.0
    # #
    glyc_frac = 1.0 - 11.0729 * exp(exp_arg * temp)
    glyc_frac = max(glyc_frac, 0.0)
    k = a0 * glyc_frac
    return k
end

function GC_GLYCOH_B_a(a0)
    # Used to compute the rate for this reaction:
    #    GLYC + OH = HCOOH + OH + CO
    # which is the "B" branch of GLYC + OH.
    #
    # For this reaction, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0 * exp(c0/T)
    # Because b0 = c0 = 0.  Therefore we can skip computing these
    # terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0
    # REAL(dp)             :: glyc_frac, k
    # REAL(dp), PARAMETER  :: 
    exp_arg = -1.0 / 73.0
    #
    glyc_frac = 1.0 - 11.0729 * exp(exp_arg * temp)
    glyc_frac = max(glyc_frac, 0.0)
    k = a0 * (1.0 - glyc_frac)
    return k
end


function GC_OHCO_a(a0)
    # Reaction rate for:
    #    OH + CO = HO2 + CO2 (cf. JPL 15-10)
    #
    # For this reaction, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0 * exp(c0/T)
    # because b0 = c0 = 0.  Therefore we can skip computing these
    # terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0
    # #
    # REAL(dp)             :: klo1,   klo2,   khi1,  khi2
    # REAL(dp)             :: xyrat1, xyrat2, blog1, blog2,   fexp1
    # REAL(dp)             :: fexp2,  kco1,   kco2,  temp300, k
    #
    klo1 = 5.9E-33 * k300_over_temp
    khi1 = 1.1E-12 * k300_over_temp^(-1.3)
    xyrat1 = klo1 * NUMDEN / khi1
    blog1 = log10(xyrat1)
    fexp1 = 1.0 / (1.0 + blog1 * blog1)
    kco1 = klo1 * NUMDEN * 0.6^fexp1 / (1.0 + xyrat1)
    klo2 = 1.5E-13
    khi2 = 2.1E+09 * k300_over_temp^(-6.1)
    xyrat2 = klo2 * NUMDEN / khi2
    blog2 = log10(xyrat2)
    fexp2 = 1.0 / (1.0 + blog2 * blog2)
    kco2 = klo2 * 0.6^fexp2 / (1.0 + xyrat2)
    k = kco1 + kco2
    return k
end

function GC_RO2NO_B1_ac(a0, c0)
    # Reaction rate for the "B" branch of these RO2 + NO reactions:
    #    MO2 + NO = CH2O + NO2 + HO2
    # in which the "a1" parameter equals exactly 1.
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1 * exp(c1/T)
    # because b0 = c0 = c1 = 0.  Therefore we can skip computing
    # these terms.  This avoids excess CPU cycles. (bmy, 1/4/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0
    # REAL(dp), PARAMETER  :: 
    one_minus_fyrno3 = 1.0 - 3.0e-4
    # REAL(dp)             :: k
    #
    k = a0 * exp(c0 / temp) * one_minus_fyrno3
    return k
end


function GCJPLPR_aba(a1, b1, a2, fv)
    # Third body effect for pressure dependence of rate coefficients.
    # a1, b1 are the Arrhenius parameters for the lower-limit rate.
    # a2     is  the Arrhenius parameters for the upper-limit rate.
    # fv     is the falloff curve paramter, (see ATKINSON ET. AL (1992)
    #        J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.
    #
    # Used to compute the rate for these reactions:
    #    OH  + OH  {+M} = H2O2
    #    NO2 + OH  {+M} = HNO3       {+M}
    #    Cl  + O2  {+M} = ClOO       {+M}
    #    SO2 + OH  {+M} = SO4  + HO2
    #    Br  + NO2 {+M} = BrNO2      {+M}
    #    NO  + O   {+M} = NO2        {+M}
    #    I   + NO2 {+M} = IONO       {+M}
    #    I   + NO  {+M} = INO        {+M}
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    exp(c1/T)
    #    (300/T)^b2 * exp(c2/T)
    # because b2 = c1 = c2 = 0.  Therefore we can skip computing these
    # terms.  Also, fct1 = fct2 = 0, so we will skip computing these
    # terms as well.  This is more computationally efficient.
    # (bmy, 1/25/20)
    #
    # REAL(dp), INTENT(IN) :: a1,   b1,    a2,   fv
    # REAL(dp)             :: rlow, xyrat, blog, fexp, k
    #
    rlow = a1 * (k300_over_temp^1) * NUMDEN
    xyrat = rlow / a2   #rhigh = a2
    blog = log10(xyrat)
    fexp = 1.0 / (1.0 + (blog * blog))
    k = rlow * (fv^fexp) / (1.0 + xyrat)
    return k
end

function GCJPLPR_abab(a1, b1, a2, b2, fv)
    # Third body effect for pressure dependence of rate coefficients.
    # a1, b1 are the Arrhenius parameters for the lower-limit rate.
    # a2, b2 are the Arrhenius parameters for the upper-limit rate.
    # fv     is the falloff curve paramter, (see ATKINSON ET. AL (1992)
    #        J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.
    #
    # Used to compute the rate for these reactions:
    #    NO   + OH  {+M} = HNO2  {+M}
    #    HO2  + NO2 {+M} = HNO4
    #    NO2  + NO3 {+M} = N2O5
    #    ClO  + NO2 {+M} = ClNO3 {+M}
    #    MCO3 + NO2 {+M} = PAN
    #    RCO3 + NO2 {+M} = PPN
    #    PRPE + OH  {+M} = PO2
    #    MO2  + NO2 {+M} = MPN   {+M}
    #    BrO  + NO2 {+M} = BrNO3 {+M}
    #    NO2  + O   {+M} = NO3   {+M}
    #    H    + O2  {+M} = HO2   {+M}
    #    IO   + NO2 {+M} = IONO2 {+M}
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    exp(c1/T)
    #    exp(c2/T)
    # because c1 = c2 = 0.  Therefore we can skip computing these
    # terms.  Also, fct1 = fct2 = 0, so we will skip computing these
    # terms as well.  This is more computationally efficient.
    # (bmy, 1/25/20)
    #
    # REAL(dp), INTENT(IN) :: a1,   b1,    a2,    b2,   fv
    # REAL(dp)             :: rlow, rhigh, xyrat, blog, fexp, k
    #
    rlow = a1 * (k300_over_temp^b1) * NUMDEN
    rhigh = a2 * (k300_over_temp^b2)
    xyrat = rlow / rhigh
    blog = log10(xyrat)
    fexp = 1.0 / (1.0 + (blog * blog))
    k = rlow * (fv^fexp) / (1.0 + xyrat)
    return k
end

function GCJPLPR_abcabc(a1, b1, c1, a2, b2, c2, fv)
    # Third body effect for pressure dependence of rate coefficients.
    # a1, b1, c1 are the Arrhenius parameters for the lower-limit rate.
    # a2, b2, c2 are the Arrhenius parameters for the upper-limit rate.
    # fv         is the falloff curve paramter, (see ATKINSON ET. AL (1992)
    #           J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.
    #
    # Used to compute the rate for these reactions:
    #    HNO4 {+M} = HO2 + NO2
    #    N2O5 {+M} = NO2 + NO3
    #    MPN  {+M} = MO2 + NO2
    #
    # REAL(dp), INTENT(IN) :: a1,   b1,    c1,    a2,   b2,   c2,  fv
    # REAL(dp)             :: rlow, rhigh, xyrat, blog, fexp, k
    #
    rlow = a1 * (k300_over_temp^b1) * exp(c1 / temp) * NUMDEN
    rhigh = a2 * (k300_over_temp^b2) * exp(c2 / temp)
    xyrat = rlow / rhigh
    blog = log10(xyrat)
    fexp = 1.0 / (1.0 + (blog * blog))
    k = rlow * (fv^fexp) / (1.0 + xyrat)
    return k
end


function GC_RO2NO_A1_ac(a0, c0)
    # Reaction rate for the "A" branch of these RO2 + NO reactions:
    #    MO2  + NO = MENO3
    # in which the "a1" parameter equals exactly 1.
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1 * exp(c1/T)
    # because b0 = b1 = c1 = 0.  Therefore we can skip computing
    # these terms.  This avoids excess CPU cycles. (bmy, 1/4/20)
    #
    # Special treatment for methyl nitrate based on observations
    # as Carter and Atkinson formulation does not apply to C1.
    # Value based on upper limit of Flocke et al. 1998 as applied
    # in Fisher et al. 2018
    #
    # REAL(dp), INTENT(IN) :: a0, c0
    # REAL(dp)             :: k
    #
    k = a0 * exp(c0 / temp) * 3.0e-4
    return k
end

function GC_RO2NO_A2_aca(a0, c0, a1)
    # Reaction rate for the "A" branch of these RO2 + NO reactions,
    #    ETO2 + NO = ETNO3
    #    A3O2 + NO = NPRNO3
    #    R4O2 + NO = R4N2
    #    B3O2 + NO = IPRNO3
    # in which the "a1" parameter is greater than 1.0.
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1 * exp(c1/T)
    # because b0 = b1 = c1 = 0.  Therefore we can skip computing
    # these terms.  This avoids excess CPU cycles. (bmy, 1/4/20)
    #
    # REAL(dp), INTENT(IN) :: a0,  c0,   a1
    # REAL(dp)             :: k0,  k, yyyn, xxyn
    # REAL(dp)             :: aaa, rarb, zzyn, fyrno3
    #
    k0 = a0 * exp(c0 / temp)
    xxyn = 1.94e-22 * exp(0.97 * a1) * NUMDEN
    yyyn = 0.826 * ((300.0 / temp)^8.1)
    aaa = log10(xxyn / yyyn)
    zzyn = (1.0 / (1.0 + (aaa * aaa)))
    rarb = (xxyn / (1.0 + (xxyn / yyyn))) * (0.411^zzyn)
    fyrno3 = (rarb / (1.0 + rarb))
    k = k0 * fyrno3
    return k
end


function GC_RO2NO_B2_aca(a0, c0, a1)
    # Reaction rate for the "B" branch of these RO2 + NO reactions:
    #    ETO2 + NO = NO2 +     HO2 + ...
    #    A3O2 + NO = NO2 +     HO2 + ...
    #    R4O2 + NO = NO2 + 0.27HO2 + ...
    #    B3O2 + NO = NO2 +     HO2 + ...
    # in which the "a1" parameter is greater than 1.0.
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1 * exp(c1/T)
    # because b0 = c0 = c1 = 0.  Therefore we can skip computing
    # these terms.  This avoids excess CPU cycles. (bmy, 1/4/20)
    #
    # Use this function when a1 input argument is greater than 1.0.
    # This avoids IF statements, which saves CPU cycles (bmy, 1/4/20)
    #
    # REAL(dp), INTENT(IN) :: a0,  c0,   a1
    # REAL(dp)             :: k0,  k, yyyn, xxyn
    # REAL(dp)             :: aaa, rarb, zzyn, fyrno3
    #
    k0 = a0 * exp(c0 / temp)
    xxyn = 1.94e-22 * exp(0.97 * a1) * NUMDEN
    yyyn = 0.826 * (k300_over_temp^8.1)
    aaa = log10(xxyn / yyyn)
    zzyn = (1.0 / (1.0 + (aaa * aaa)))
    rarb = (xxyn / (1.0 + (xxyn / yyyn))) * (0.411^zzyn)
    fyrno3 = (rarb / (1.0 + rarb))
    k = k0 * (1.0 - fyrno3)
    return k
end


function GC_OHHNO3_acacac(a0, c0, a1, c1, a2, c2)
    # Used to compute the rate for these reactions:
    #    HNO3  + OH = H2O + NO3
    #    HONIT + OH = NO3 + HAC
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1
    #    (300/T)^b2
    # Because b0 = b1 = b2 = 0.  Therefore we can skip computing
    # these terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0, a1, c1, a2, c2
    # REAL(dp)             :: k0, k1, k2, k
    #
    # ---  OH + HNO3:   K = k0 + K3[M] / (1 + K3[M]/K2)  ------
    k0 = a0 * exp(c0 / temp)
    k1 = a1 * exp(c1 / temp)
    k2 = NUMDEN * (a2 * exp(c2 / temp))
    k = k0 + k2 / (1.0 + k2 / k1)
    return k
end


function GCJPLEQ_acabab(a0, c0, a1, b1, a2, b2, fv)
    # Calculates the equilibrium constant
    # Find the backwards reaction by K=kforward/kbackwards
    # Calculates the rate constant of the forward reaction
    #
    # Used to compute the rate for these reactions:
    #    PPN        = RCO3 + NO2
    #    PAN        = MCO3 + NO2
    #    ClOO  {+M} = Cl   + O2 {+M}
    #    Cl2O2 {+M} = 2ClO      {+M}
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    exp(c1/T)
    #    exp(c2/T)
    # because b0 = c1 = c2 = 0.  Therefore we can skip computing these terms.
    # Also, fct1 = fct2 = 0, so we will skip those terms as well.  This is
    # more computationally efficient. (bmy, 1/25/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0, a1, b1, a2, b2, fv
    # REAL(dp)             :: k0, k1, k
    #
    k0 = a0 * exp(c0 / temp)               # backwards rxn rate
    k1 = GCJPLPR_abab(a1, b1, a2, b2, fv)  # forwards rxn rate
    k = k1 / k0
    return k
end

end