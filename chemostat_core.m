%this is the core of the chemostat model:
%differential equations and their numerical integration

%Initial Conditions from main script set here:
O2=O2o;
Sd=Sdo;
S=So;
NO3=NO3o;
NO2=NO2o;
NH4=NH4o;
N2=0; %so no confusion with nitrous oxide N2O

%store the initial condition for all in the time resolution:
St=nan*ones(1,trecd);
Sdt=nan*ones(1,trecd);
NO3t=nan*ones(1,trecd);
NO2t=nan*ones(1,trecd);
NH4t=nan*ones(1,trecd);
N2t=nan*ones(1,trecd);
O2t=nan*ones(1,trecd);
XQot=nan*ones(1,trecd);
XQn1t=nan*ones(1,trecd);
XQn23t=nan*ones(1,trecd);
XQoPt=nan*ones(1,trecd);
XQn1Pt=nan*ones(1,trecd);
XQn23Pt=nan*ones(1,trecd);
XQbnh4t=nan*ones(1,trecd);
XQbno2t=nan*ones(1,trecd);
XQat=nan*ones(1,trecd);
XQfact=nan*ones(1,trecd);
XQfacrespt=nan*ones(1,trecd);
XQanit=nan*ones(1,trecd);
timet=nan*ones(1,trecd);

%record the first step of time resolution
St(1)=S;
Sdt(1)=Sd;
NO3t(1)=NO3;
NO2t(1)=NO2;
NH4t(1)=NH4;
N2t(1)=N2;
O2t(1)=O2;
XQot(1)=XQo;
XQn1t(1)=XQn1;
XQn23t(1)=XQn23;
XQoPt(1)=XQoP;
XQn1Pt(1)=XQn1P;
XQn23Pt(1)=XQn23P;
XQbnh4t(1)=XQbnh4;
XQbno2t(1)=XQbno2;
XQat(1)=XQa;
XQfact(1)=XQfac;
XQanit(1)=XQani;
timet(1)=0;

%need to record fac resp at all timesteps
XQfacresp=nan*ones(1,tstepmax);

%%
%start time loop
for t=1:tstepmax
    
    
    %uptake rates:
    pO2=po_coef*O2;
    pS=VmaxS*S./(Ks+S);
    pSd=VmaxS*Sd./(Ks+Sd);
    pNO3=VmaxN*NO3./(Kn+NO3);
    pNO2=VmaxN*NO2./(Kn+NO2);
    pNH4_AOO=VmaxN_AOO*NH4./(Kn_AOO+NH4);
    pNH4=VmaxN*NH4./(Kn+NH4);
    
    %growth rates:
    uo=min(pSd*yod,pO2*yoe);
    un1=min(pSd*yn1d,pNO3*yn1e);
    un23=min(pSd*yn23d,pNO2*yn23e);
    un123=min(pSd*yn123d,pNO3*yn123e);
    %PA-assoc
    uoP=min(pS*yod,pO2*yoe);
    un1P=min(pS*yn1d,pNO3*yn1e);
    un23P=min(pS*yn23d,pNO2*yn23e);
    un123P=min(pS*yn123d,pNO3*yn123e);
    ubnh4=min(pNH4_AOO*ynh4_bnh4,pO2*yo_bnh4);
    ubno2=min(pNO2*yno2_bno2,pO2*yo_bno2);
    ua=min(pNO2*yno2_a,pNH4*ynh4_a);
    uani=pNO2*yno2_ani; %only NO2 limited
    
    %fac:
    uofac=min(pSd*yodfac,pO2*yoefac);
    unfac = min(pSd*ynfacd,pNO3*ynface);
    
    %facultative growth rate, using either O2 or NO3: (could add NO2, too!)
    if uofac>unfac
        ufac=uofac;
        yfacd=yodfac;
        facO2use=ufac*XQfac/yoefac;
        facNO3use=0;
        facNO2prod=0;
        XQfacresp(t)=1; %1 for using oxygen
    elseif unfac>=uofac
        ufac=unfac;
        yfacd=ynfacd; %for NH4 excr
        facO2use=0;
        facNO3use=ufac*XQfac/ynface;
        facNO2prod=ufac*XQfac*enfac;
        XQfacresp(t)=0; %1 for using oxygen
    end
    
    %change in external nutrients:
    ddtO2 = D*(O2o - O2) ...
        - uo*XQo/yoe ...
        - uoP*XQoP/yoe ...
        - facO2use ...
        - ubnh4*XQbnh4/yo_bnh4 ...
        - ubno2*XQbno2/yo_bno2 ...
        ;
    O2 = O2 + ddtO2*dt;
    
    ddtSd = D*(Sdo - Sd) ...
        - uo*XQo/yod ...
        - ufac*XQfac/yfacd ...
        - un1*XQn1/yn1d ...
        - un23*XQn23/yn23d ...
        - un123*XQn123/yn123d ...
        ;
    Sd = Sd + ddtSd*dt;
    
    ddtS = D*(So - S) ...
        - uoP*XQoP/yod ...
        - un1P*XQn1P/yn1d ...
        - un23P*XQn23P/yn23d ...
        - un123P*XQn123P/yn123d ...
        ;
    S = S + ddtS*dt;
    
    ddtNO3 = D*(NO3o - NO3) ...
        - un1*XQn1/yn1e ...
        - un123*XQn123/yn123e ...
        - un1P*XQn1P/yn1e ...
        - facNO3use ...
        + ubno2*XQbno2/(1/yno2_bno2-1) ...
        + uani*XQani/yno2_ani ...
        + ua*XQa*eno3a ...
        ;
    NO3 = NO3 + ddtNO3*dt;
    
    ddtNO2 = D*(NO2o - NO2) ...
        - un23*XQn23/yn23e ...
        - un23P*XQn23P/yn23e ...
        - ubno2*XQbno2*(1/yno2_bno2-1) ... %bc it uses 1 mol nh4 for synth
        - ua*XQa/yno2_a ...
        - uani*XQani/yno2_ani ...
        + un1*XQn1*en1 ...
        + un1P*XQn1P*en1 ...
        + facNO2prod ...
        + ubnh4*XQbnh4*(1/ynh4_bnh4-1) ...
        ;
    NO2 = NO2 + ddtNO2*dt;
    
    ddtNH4 = D*(NH4o - NH4) ...
        - ubnh4*XQbnh4/ynh4_bnh4 ...
        - ua*XQa/ynh4_a ...
        - ubno2*XQbno2 ... %bc it uses NH4 for synthesis
        + uo*XQo*(1/yod-1) ...
        + un1*XQn1*(1/yn1d-1) ...
        + un23*XQn23*(1/yn23d-1) ...
        + un123*XQn123*(1/yn123d-1) ...
        + uoP*XQoP*(1/yod-1) ...
        + un1P*XQn1P*(1/yn1d-1) ...
        + un23P*XQn23P*(1/yn23d-1) ...
        + un123P*XQn123P*(1/yn123d-1) ...
        + ufac*XQfac*(1/yfacd-1) ...
        ;
    NH4 = NH4 + ddtNH4*dt;
    
    ddtN2 =  -D*N2 ...
        + un23*XQn23*en23 ...
        + un123*XQn123*en123 ...
        + un23P*XQn23P*en23 ...
        + un123P*XQn123P*en123 ...
        + ua*XQa*ena ...
        ;
    N2 = N2 + ddtN2*dt;
    
    %change in biomass:
    ddtXQo=uo*XQo - D*XQo;% - 0.05*XQo;
    XQo=XQo+ddtXQo*dt;
    
    ddtXQn1=un1*XQn1 - D*XQn1;% - .05*XQn1;
    XQn1=XQn1+ddtXQn1*dt;
    
    ddtXQn23=un23*XQn23 - D*XQn23;% - .05*XQn23;
    XQn23=XQn23+ddtXQn23*dt;
    
    ddtXQn123=un123*XQn123 - D*XQn123;% - .05*XQn23;
    XQn123=XQn123+ddtXQn123*dt;
    
    ddtXQoP=uoP*XQoP - D*XQoP;
    XQoP=XQoP+ddtXQoP*dt;
    
    ddtXQn1P=un1P*XQn1P - D*XQn1P;% - .05*XQn1;
    XQn1P=XQn1P+ddtXQn1P*dt;
    
    ddtXQn23P=un23P*XQn23P - D*XQn23P;% - .05*XQn23;
    XQn23P=XQn23P+ddtXQn23P*dt;
    
    ddtXQn123P=un123P*XQn123P - D*XQn123P;% - .05*XQn23;
    XQn123P=XQn123P+ddtXQn123P*dt;
    
    ddtXQbnh4=ubnh4*XQbnh4 - D*XQbnh4;
    XQbnh4=XQbnh4+ddtXQbnh4*dt;
    
    ddtXQbno2=ubno2*XQbno2 - D*XQbno2;
    XQbno2=XQbno2+ddtXQbno2*dt;
    
    ddtXQa=ua*XQa - D*XQa;
    XQa=XQa+ddtXQa*dt;
    
    ddtXQani=uani*XQani - D*XQani;
    XQani=XQani+ddtXQani*dt;
    
    ddtXQfac=ufac*XQfac - D*XQfac;
    XQfac=XQfac+ddtXQfac*dt;
    
    %record at trec for time-varying resolution
    if mod(t,trec)==0
        
        j=t/trec;
        timet(j)=t*dt;
        
        disp(['O2o=' num2str(O2o) ' ' num2str(t*dt)])
        %disp(j)
        
        St(j)=S;
        Sdt(j)=Sd;
        NO3t(j)=NO3;
        NO2t(j)=NO2;
        NH4t(j)=NH4;
        N2t(j)=N2;
        O2t(j)=O2;
        XQot(j)=XQo;
        XQn1t(j)=XQn1;
        XQn23t(j)=XQn23;
        XQoPt(j)=XQoP;
        XQn1Pt(j)=XQn1P;
        XQn23Pt(j)=XQn23P;
        XQbnh4t(j)=XQbnh4;
        XQbno2t(j)=XQbno2;
        XQat(j)=XQa;
        XQanit(j)=XQani;
        XQfact(j)=XQfac;
        
    end
    
end %end population model time t loop

%%
%display end state solution
disp('IC')
disp(['O2 in: ' num2str(O2o)])
disp(['OrgM in: ' num2str(So)])
disp('NUTRIENTS')
disp(['O2: ' num2str(O2)])
disp(['OrgM: ' num2str(S)])
disp(['DOrgM: ' num2str(Sd)])
disp(['NO3: ' num2str(NO3)])
disp(['NH4: ' num2str(NH4)])
disp(['NO2: ' num2str(NO2)])
disp(['N2: ' num2str(N2)])
disp('BIOMASSES')
disp(['Aerobic het: ' num2str(XQo)])
disp(['Facult. het: ' num2str(XQfac)])
disp(['Nitr red het: ' num2str(XQn1)])
disp(['Denitr het: ' num2str(XQn23)])
disp(['Full Denitr het: ' num2str(XQn123)])
disp(['AOO: ' num2str(XQbnh4)])
disp(['NOO: ' num2str(XQbno2)])
disp(['anammox: ' num2str(XQa)])

%Total fraction of aerobic faculative resp:
avgl=length(XQfacresp); %length of the record
facaer=nanmean(XQfacresp(avgl/2:end-1));
disp(['Fac: fraction of aerobic resp: ' num2str(facaer)])

%%
%THIS took time, showed a very smooth cumulative mean convering
%to the actual mean!
% %cumulative mean
% mfac=[];
% for t=1:100:length(XQfacrespt)
%     temp=XQfacrespt(1:t);
%     mfac(end+1)=sum(temp)/length(temp);
% end

%%

%plot individual runs over time
if plotindrun==1
    
    figure;
    g=2;
    subplot(4,1,1);plot(timet,[O2t; Sdt; St],'LineWidth',g);
    title('External Substrates ','FontSize',16);
    ylabel('uM ','FontSize',14);grid on;
    legend('O2 ','S ')
    ylim([0 .2])
    
    subplot(4,1,2);plot(timet,[NH4t; NO2t; NO3t; N2t],'LineWidth',g);
    title('External Substrates ','FontSize',16);
    ylabel('uM ','FontSize',14);grid on;
    legend('NH4 ', 'NO2 ', 'NO3 ', 'N2 ')
    ylim([0 5])
    
    %Biomass/abundances
    subplot(4,1,3);plot(timet,[XQot; XQfact; XQn1t; XQn23t; XQbnh4t; XQbno2t; XQat; XQanit],'LineWidth',g);
    %subplot(4,1,4);plot(time,[ninds(1,:)/vol(1); X(1,:)],'b',time,[ninds(2,:)/vol(2); X(2,:)],'g',time,[ninds(3,:)/vol(3); X(3,:)],'r','LineWidth',g);
    title('Biomass ','FontSize',16);
    %ylabel('cells/m^3 ','FontSize',14);grid on;
    ylabel('\muM N ','FontSize',14);grid on;
    %legend('Aerobes ','n123 ','n1 ','n2 ','n3 ','n12 ','n23 ','Anammox ',12);
    legend('Het: Aerobic','Het: Facult.','Het: NO3-red','Het: Denitr','AOO','NOO','Anammox','Anitox');
    
    %Rates
    subplot(4,1,4);plot(timet,[XQot*uo/yoe+XQfact*uo/yoe.*XQfacrespt; XQn1t*un1/yn1e+XQfact*un1/yn1e.*(-XQfacrespt+1); XQn23t*un23/yn23e; XQbnh4t*ubnh4/ynh4_bnh4; XQbno2t*ubno2/yno2_bno2; XQat*ua*ena; XQanit*uani/yno2_ani]*1e3,'LineWidth',g);
    %subplot(4,1,4);plot(time,[ninds(1,:)/vol(1); X(1,:)],'b',time,[ninds(2,:)/vol(2); X(2,:)],'g',time,[ninds(3,:)/vol(3); X(3,:)],'r','LineWidth',g);
    title('Rates ','FontSize',16);
    %ylabel('cells/m^3 ','FontSize',14);grid on;
    ylabel('nM/d ','FontSize',14);grid on;
    %legend('Aerobes ','n123 ','n1 ','n2 ','n3 ','n12 ','n23 ','Anammox ',12);
    legend('Het: Aerobic Resp rate (incl Fac) (nM O2/d?)','Het: NO3-red (incl fac) (nM N/d)','Het: Denitr rate','NH4 oxid rate','NO2 oxid rate','Anammox rate');
    
end

%%
%some diagnostics:

disp('frac Anammox of tot denitr rate ratio:')
disp((XQa*ua*ena)/(XQn123*un123/yn123e+XQn23*un23/yn23e+XQa*ua*ena))

disp('Tot N loss (uM N/d):')
disp((XQn123*un123/yn123e+XQn23*un23/yn23e+XQa*ua*ena))

