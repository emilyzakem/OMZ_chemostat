%traits (parameters) for metabolisms: uptake kinetics and yields

%1. Uptake params:

%Oxygen:
po_coef=2329.100;%1/uM/day - see Zakem and Follows 2016 for diffusive uptake parameterization

%organic matter (uncertain/somewhat arbitrary)
VmaxS = 1; %mol org N/mol cell N/day
Ks=0.1; %uM

%DIN:
%high affinity AOA:
VmaxN_AOO = 50.8; %mol DIN/mol cell N/day at 20C for AOA (see Zakem et al 2018)
Kn_AOO=0.133; %uM for AOA
%general for all:
VmaxN = VmaxN_AOO;
Kn=Kn_AOO;

%2. yields y and products/excretions e:
%excretion of NH4 (e) for all of the hets is always (1/ynd - 1)

%denominator coefs for OrgM half-reactions (see Zakem et al 2019 ISME supplement)
dC = 20; %denominator for biomass eqnas
dD = 29.1; %denominator for OM eqns for Redfieldian OM

%Aerobic hetetrotrophy
yod=0.14; %OrgM yield, &!mol biomass N/mol organic N - open oc avg Robinson 2009
f=yod*dC/dD; %associated electron fraction
yoe=4*f/dC/(1-f); %O2 yield, &!mol biomass N/mol O2

%Nitrate reduction: NO3 to NO2
yn1d=yod*0.9; %OrgM yield, assume 10% as efficient (energetics: 95%. obs: lower)
f=yn1d*dC/dD;
yn1e=2*f/dC/(1-f); %NO3 yield
en1=1/yn1e; %excretion of NO2 (mol NO2 formed per mol biomass N synthesized)

%Denitrification: NO2 to N2
yn23d=yod*0.9; %OrgM yield
f=yn23d*dC/dD;
yn23e=3*f/dC/(1-f); %NO2 yield
en23=1/yn23e; %mol N (as N2) formed per mol biomass N synthesized (careful w N2 vs N)

%FULL denitrification: NO3 to N2
yn123d=yod*0.9; %OrgM yield
f=yn123d*dC/dD;
yn123e=5*f/dC/(1-f); %NO3 yield
en123=1/yn123e; %mol N (as N2) formed per mol biomass N synthesized (careful w N2 vs N)

%Chemoautotr ammonia oxidation (nitrification step 1: NH4 to NO2) --
%see Zakem et al 2018 Nat Comm for derivations
ynh4_bnh4=1/112; %mol biomass N per mol NH4 consumed
yo_bnh4=1/162; %mol biomass N per mol O2 consumed

%Chemoautotr nitrite oxidation (nitrification step 2: NO2 to NO3)
yno2_bno2=1/334; %mol biomass N per mol NH4 consumed
yo_bno2=1/162; %mol biomass N per mol O2 consumed

%Chemoautotr anammox (NH4 + NO2 -> NO3 + N2)
%using x=0.5 and f = 0.03 (see Zakem et al ISME 2019)
ynh4_a = 1/154; %mol biomass N per mol NH4 consumed
yno2_a = 1/216; %mol biomass N per mol NO2 consumed
ena = 163*2; %mol N (as N2) formed per mol biomass N synthesized (careful w N2 vs N)
eno3a = 42; %mol NO3 formed per mol biomass N synthesized

%Facultative heterotroph with penalties for mixotrophy:
facpen = 0.8; %penalty w/r/to aerobic heterotroph
yodfac=yod*facpen;%, & !mol biomass N/mol organic matter - aerobic bacteria yield
yoefac=4*yodfac/dC/(1-yodfac);%, & !mol B/mol O2 -where cells have 1mol N
ynfacd=yn123d*facpen; %ynd/113.22*82.65 %mol B/mol OM from van de leemput
ynface=5*ynfacd/dC/(1-ynfacd); %mol B/mol NO3
enfac=1/ynface; %mol NO2 formed/mol B (no NO2 if it's 123)

%anitox (not developed)
yno2_ani = 0;