%Created by: Emily Zakem, 11/22/19
%(note to self: copy from Mean1D_anxNO3_Jul2019_Ch3Fig4_ChemostatLoopAll.m)
%This version was published (as an ensemble with parameter variation) in:
%Zakem et al. "Stable aerobic and anaerobic coexistence in anoxic marine zones"
%ISME 2019

%This code solves for the steady state N cycle by numerically integrating
%forward in time in a virtual chemostat for multiple specified ratios of O2 to
%organic matter (OM) supply. It is currently set up so that O2 supply varies (looping through
%varying values of O2:OM). As is, OM supply is constant, so total biomass remains
%relatively constant, but can also explore 2D solutions with OM varying as
%well.

%Variables:

%nutrients:
%Sd %organic matter substrate (originally "dissolved")
%S %a second org matter type, set to 0 here (originally "particulate")
%O2, NH4, NO2, NO3, and N2 (the latter is not treated as a gas here, so not necessarily helpful)

%metabolic functional type possibilities as biomass "XQ":
%IMPORTANT: to turn on and off types, set initial biomass conditions below (l. 143)
%right now, 6 types are on, as in Fig. 3 in Zakem et al ISME 2019

%Heterotrophs:

%1. aerobic hets: XQo
%2. nitrate reducing hets (NO3 to NO2): XQn1
%3. denitrifying hets (NO2 to N2): XQn23
%4. full denitrification hets (NO3 to N2): XQn123 (this label was used for DNRA in other version, fyi)

%Redundant heterotrophs for organic substrate 2 ("P" was originally for particle)
%5. aerobic hets: XQoP
%6. nitrate reducing hets (NO3 to NO2): XQn1P
%7. denitrifying hets (NO2 to N2): XQn23P
%8. full denitrification hets (NO3 to N2): XQn123P

%Chemoautotrophs
%9. ammonia oxidizer: XQbnh4
%10. nitrite oxidizer: XQbno2
%11. annamox: XQa

%Additional types that are not used here and may be incompletey represented
%check over before using (anitox is definitely not represented throughout)
%12. "anitox": XQani 
%13. facultative type XQfac that carries out 1 and 2 above (aerobic het and
%nitrate-reducing het)

%%

%options for running:
runnewloop=1; %1 to run a new loop, 0 to load saved variables from directory dirload
savenewfiles=0; %1 to save new files in directory dirsave (does not run if runnewloop=0)
plotindrun = 0; %plot the results over time for each run

if runnewloop == 0
    dirload = 'Run_Ex';
    disp(['Loading and plotting from: ',dirload])
else
    if savenewfiles==1
        dirsave = 'Run_Ex2'; %this needs to not already exist, else it gives the following error
        if isdir(dirsave) %checks for existence to avoid copying over old runs
            error(['Directory ',dirsave,' (dirsave) already exists: Stopping.', ...
                ' To proceed, change dirsave or rename directory'])
        else
            mkdir(dirsave)
        end
    end
end

%Set incoming supply concentration(s) of OM in uM nitrogen (note this could be a vector)
Sdoall=1; %muM N
Soall=0;
%makePOMforHigherO2crit=0; %1 for o2crit to affect hets (as if inside particles), 0 for not.

%Set incoming supply concentration(s) of O2 in log10(uM)
%O2o_exp=[1.5 0.5 -2]; %quick 3 oxygen thresholds (decreasing supply) to test
O2o_exp=[2 1 0 -1 -2]; %quick 5
%O2o_exp=2:-0.1:-2; %higher res
%this one published in ISME:
%O2o_exp=cat(2,[2:-0.1:1.1],[1:-0.025:0],[-0.1:-0.1:-2]); %zooms in
% to intermediate state

O2oall=10.^O2o_exp; %actual O2 incoming concs in uM

vary1=O2oall;
vary2=Sdoall;

%model params:
D=0.05; %dilution rate (1/day)

ndays=1e4; %total # of days to run chemostat. 1e4 d is WAY more than enough time to equilibrate 
dt=0.001; %timestep (days)
tstepmax = ndays/dt; %# of timesteps
trecd=200; %(days) record the time resolution of individual runs this often
trec=tstepmax/trecd; %for recording time resolution

%load in the traits: uptake kinetic parameters and yields
chemostat_traits

%%
if runnewloop==1
    
    %set up empties for loops
    %Nutrients:
    O2all=nan*ones(length(vary1),length(vary2));
    Sall=nan*ones(length(vary1),length(vary2));
    Sdall=nan*ones(length(vary1),length(vary2));
    NO3all=nan*ones(length(vary1),length(vary2));
    NO2all=nan*ones(length(vary1),length(vary2));
    NH4all=nan*ones(length(vary1),length(vary2));
    N2all=nan*ones(length(vary1),length(vary2));
    %Biomasses:
    %Hets:
    XQoall=nan*ones(length(vary1),length(vary2));
    XQn1all=nan*ones(length(vary1),length(vary2));
    XQn23all=nan*ones(length(vary1),length(vary2));
    XQn123all=nan*ones(length(vary1),length(vary2));
    %PA:
    XQoPall=nan*ones(length(vary1),length(vary2));
    XQn1Pall=nan*ones(length(vary1),length(vary2));
    XQn23Pall=nan*ones(length(vary1),length(vary2));
    XQn123Pall=nan*ones(length(vary1),length(vary2));
    %Chemoautos:
    XQbnh4all=nan*ones(length(vary1),length(vary2));
    XQbno2all=nan*ones(length(vary1),length(vary2));
    XQaall=nan*ones(length(vary1),length(vary2));
    %Others:
    XQfacall=nan*ones(length(vary1),length(vary2));
    XQaniall=nan*ones(length(vary1),length(vary2));
    %to track facultative average respiration (i.e. use of oxygen vs nitrate as e- acceptor)
    facaerall=nan*ones(length(vary1),length(vary2));
    
%start loop
    for k=1:length(vary1)
        for m=1:length(vary2)
            
            %Chemostat input supply concentrations (all in uM):
            O2o=O2oall(k); %incoming O2 supply for loop k
            Sdo=Sdoall(m); %incoming org matter supply for loop m
            So=0; %incoming org matter #2 -- turning this one off here
            NO3o=30; %simulates deep anoxic zone
            NO2o=0;
            NH4o=0;
            
            %Initial conditions for biomass: Set to 0 to turn off a type
            %heterotroph:
            XQo=0.1;
            XQn1=0.1;
            XQn23=0.1;
            XQn123=0;%.1;
            %PA heterotroph:
            XQoP=0;%.1;
            XQn1P=0;%.1;
            XQn23P=0;%.1;
            XQn123P=0;%.1;
            %chemoauto:
            XQbnh4=0.1;
            XQbno2=0.1;
            XQa=0.1;
            %other:
            XQfac=0;%.1;
            XQani=0;%.1;
            
            %the model:
            chemostat_core
            
            %record solutions:
            %Nutrients:
            O2all(k,m)=O2;
            Sall(k,m)=S;
            Sdall(k,m)=Sd;
            NO3all(k,m)=NO3;
            NO2all(k,m)=NO2;
            NH4all(k,m)=NH4;
            N2all(k,m)=N2;
            %Biomasses
            XQoall(k,m)=XQo;
            XQn1all(k,m)=XQn1;
            XQn23all(k,m)=XQn23;
            XQn123all(k,m)=XQn123;
            XQoPall(k,m)=XQoP;
            XQn1Pall(k,m)=XQn1P;
            XQn23Pall(k,m)=XQn23P;
            XQn123Pall(k,m)=XQn123P;
            XQbnh4all(k,m)=XQbnh4;
            XQbno2all(k,m)=XQbno2;
            XQaall(k,m)=XQa;
            XQfacall(k,m)=XQfac;
            XQaniall(k,m)=XQani;
            %fac resp
            facaerall(k,m)=facaer;
            
        end %end m loop
    end %end k loop
    
    if savenewfiles==1
        chemostat_savefiles
    end 
    
else
    chemostat_loadfiles
end

%%
%calculate subsistence concentrations from parameter values (note that we
%don't need solutions for this)

O2s_bo=D/po_coef/yoe;
Ss_bo=D*Ks/(VmaxS*yod-D);

Ss_n1=D*Ks/(VmaxS*yn1d-D);
NO3s_n1=D*Kn/(VmaxN*yn1e-D);

Ss_dnra=D*Ks/(VmaxS*yn1d-D);
NO3s_dnra=D*Kn/(VmaxN*yn1e-D);

Ss_n23=D*Ks/(VmaxS*yn23d-D);
NO2s_n23=D*Kn/(VmaxN*yn23e-D);

NH4s_bnh4=D*Kn/(VmaxN*ynh4_bnh4-D);
O2s_bnh4=D/po_coef/yo_bnh4;

NO2s_bno2=D*Kn/(VmaxN*yno2_bno2-D);
O2s_bno2=D/po_coef/yo_bno2;

NH4s_a=D*Kn/(VmaxN*ynh4_a-D);
NO2s_a=D*Kn/(VmaxN*yno2_a-D);

%calculate phi, the threshold for anaerobic metabolism (crude: based on aerobic heterotroph bc
%we want to describe ultimately as being controlled by organic matter.
%anx pops up just before this threshold bc of higher O2* of AOA)

phi = O2s_bo + yod/yoe*(Sdoall - Ss_bo); %uM O2


%%
%Plot:
chemostat_plot

