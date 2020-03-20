%1D:
if length(Sdoall)==1
    
    figure;
    set(gcf,'color','w')
    %%
    g=2;
    subplot(5,2,1);
    %plot(O2oall,[O2all'; Sdall'; Sall'; NO3all'],'-','LineWidth',g);
    plot(O2oall,O2all','b-','LineWidth',g);
    hold on
    plot(O2oall,repmat(O2s_bo,1,length(O2oall)),'b--','LineWidth',2);
    %plot(O2oall,repmat(Ss_bo,1,length(O2oall)),'b--','LineWidth',1);
    %plot(O2oall,repmat(Ss_n1,1,length(O2oall)),'r--','LineWidth',1);
    %plot(O2oall,repmat(Ss_n23,1,length(O2oall)),'r--','LineWidth',1);
    %plot(O2oall,repmat(O2s_bnh4,1,length(O2oall)),'m--','LineWidth',1);
    %light blue pop color:
    %        plot(O2oall,repmat(O2s_bno2,1,length(O2oall)),'--','color',[30,144,255]/255,'LineWidth',2);
    %nut color:
    plot(O2oall,repmat(O2s_bno2,1,length(O2oall)),'b:','LineWidth',2);
    plot([(phi) (phi)],[1e-5 1e3],'k-.','LineWidth',1)
    %plot([(phi_bnh4) (phi_bnh4)],[1e-5 1e3],'b-.','LineWidth',1)
    %plot([(phi_bno2) (phi_bno2)],[1e-5 1e3],'r-.','LineWidth',1)
    hold off
    title('a. Oxygen','FontSize',16);
    %ylabel('\muM O_2 ','FontSize',14);grid on;
    %legend('O2 ','Sdom ','Spom','NO3','\Phi=1')
    legend('O_2 ','O_2^* (B_H_e_t_O) ','O_2^* (B_A_O_O and B_N_O_O) ','\Phi=1')
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    set(gca,'xdir','reverse')
    set(gca,'fontsize',16)
    ylim([1e-4 1e2])
    %legend(['D=' num2str(D(1))],['The Blend '],['D=' num2str(D(3))]); %ylim([0 1000])
    set(gca,'xtick',[1e-2,0.1,1,10,100])
    set(gca,'xticklabel','')
    legend boxoff
    set(gca,'ytick',[1e-3,0.1,10])
    set(gca,'yTickLabel',str2mat('1 nM','100 nM','10 \mu M'))
    grid on
    
    subplot(5,2,3);
    
    %        hold on
    plot(O2oall,NH4all','k','LineWidth',g);%orange: [1,0.5,0]
    hold on
    plot(O2oall,NO2all','r','LineWidth',g);
    %pop colors:
    %plot(O2oall,repmat(NH4s_bnh4,1,length(O2oall)),'--','color',[30,144,255]/255,'LineWidth',2);
    %plot(O2oall,repmat(NO2s_bno2,1,length(O2oall)),'--','color',[135,206,250]/250,'LineWidth',2);
    %plot(O2oall,repmat(NH4s_a,1,length(O2oall)),'--','color',[255,192,203]/255,'LineWidth',2);
    %nut colors:
    plot(O2oall,repmat(NH4s_bnh4,1,length(O2oall)),'k--','LineWidth',2);
    plot(O2oall,repmat(NO2s_bno2,1,length(O2oall)),'r--','LineWidth',2);
    plot(O2oall,repmat(NH4s_a,1,length(O2oall)),'k:','LineWidth',2);
    %almost indistinguishable from NH4s_a:
    plot(O2oall,repmat(NO2s_a,1,length(O2oall)),'r:','LineWidth',1);
    hold off
    title('b. DIN ','FontSize',16);
    ylabel('\muM N ','FontSize',14);grid on;
    legend('NH_4 ','NO_2 ','NH_4^* (B_A_O_O)','NO_2^* (B_N_O_O)','NO_2^* and NH_4^* (B_a_n_x)')
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    set(gca,'xdir','reverse')
    ylim([1e-2 2])%1e1])
    set(gca,'ytick',[1e-2,1e-1,1,10])
    set(gca,'yTickLabel',str2mat('0.01','0.1','1','10'))
    set(gca,'fontsize',16)
    set(gca,'xtick',[1e-2,0.1,1,10,100])
    set(gca,'xticklabel','')
    legend boxoff
    grid on
    
    subplot(5,2,5);
    %plot(O2oall,([XQoall'+XQoPall'; XQn1all'+XQn1Pall'; XQn23all'+XQn23Pall'; ...
    %XQbnh4all'; XQbno2all'; XQaall']),'-','LineWidth',g);
    plot(O2oall,(XQoall'+XQoPall'),'b','LineWidth',g);
    hold on
    plot(O2oall,XQbnh4all','color',[30,144,255]/255,'LineWidth',g);
    plot(O2oall,XQbno2all','color',[135,206,250]/250,'LineWidth',g);
    plot(O2oall,(XQn1all'+XQn1Pall'),'r','LineWidth',g+1);
    plot(O2oall,(XQn23all'+XQn23Pall'),'color',[250,128,114]/250,'LineWidth',g);
    plot(O2oall,XQaall','color',[255,192,203]/255,'LineWidth',g);
    hold off
    title('c. Biomass ','FontSize',16);
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    set(gca,'xdir','reverse')
    %ylabel('cells/m^3 ','FontSize',14);grid on;
    ylabel('\muM N ','FontSize',14);grid on;
    %legend('Aerobes ','n123 ','n1 ','n2 ','n3 ','n12 ','n23 ','Anammox ',12);
    legend('B_H_e_t_O','B_A_O_O','B_N_O_O','B_H_e_t_N_O_3','B_H_e_t_N_O_2','B_a_n_x');
    %legend('Aerobic heterotrophy','Nitrate reduction','Denitrification','Ammonia oxidation','Nitrite oxidation','Anammox');
    ylim([1e-4 1]) %for log conc
    set(gca,'fontsize',16)
    set(gca,'xtick',[1e-2,0.1,1,10,100])
    set(gca,'xticklabel','')
    set(gca,'ytick',[1e-3,1e-2,0.1])
    %set(gca,'yTickLabel',str2mat('10^-^5','10^-^3','0.1'))
    legend boxoff
    grid on
    
    
    %RATES:
    subplot(5,2,7);
    %plot(O2oall,([D/yoe*(XQoall'+XQoPall'); D/yn1e*(XQn1all'+XQn1Pall'); ...
    %    D/yn23e*(XQn23all'+XQn23Pall'); D/ynh4_bnh4*XQbnh4all'; ...
    %    D/yno2_bno2*XQbno2all'; D*ena*XQaall']*1e3),'-','LineWidth',g);
    plot(O2oall,1e3*D/yoe*(XQoall'+XQoPall'),'b','LineWidth',g);
    hold on
    plot(O2oall,1e3*D/ynh4_bnh4*XQbnh4all','color',[30,144,255]/255,'LineWidth',g);
    plot(O2oall,1e3*D/yno2_bno2*XQbno2all','color',[135,206,250]/250,'LineWidth',g);
    plot(O2oall,1e3*D/yn1e*(XQn1all'+XQn1Pall'),'r','LineWidth',g+1);
    plot(O2oall,1e3*D/yn23e*(XQn23all'+XQn23Pall'),'color',[250,128,114]/250,'LineWidth',g);
    plot(O2oall,1e3*D*ena*XQaall','color',[255,192,203]/255,'LineWidth',g);
    hold off
    title('d. Respiration rates ','FontSize',16);
    %set(gca,'yscale','log')
    set(gca,'xscale','log')
    set(gca,'xdir','reverse')
    %ylabel('cells/m^3 ','FontSize',14);grid on;
    ylabel('nM d^-^1 ','FontSize',14);grid on;
    %legend('Aerobes ','n123 ','n1 ','n2 ','n3 ','n12 ','n23 ','Anammox ',12);
    legend('Aerobic heterotrophic O_2 cons.','NH_4^+ oxidation','NO_2^- oxidation','NO_3^- reduction','Denitrification','Anammox');
    legend boxoff
    ylim([0 400]) %for rates
    %xlim([-2 2])
    set(gca,'fontsize',16)
    set(gca,'xtick',[1e-2,0.1,1,10,100])
    set(gca,'xticklabel','')
    %set(gca,'ytick',[0,125,250,375,500])
    %set(gca,'yTickLabel',str2mat('0','','250','','500'))
    %set(gca,'ytick',[0,125,250,375,500])
    %set(gca,'yTickLabel',str2mat('0','','250','','500'))
    
    %legend boxoff
    grid on
    
    fanx=ena*XQaall./(ena*XQaall + 1/yn23e*(XQn23all+XQn23Pall));
    subplot(5,2,9);
    plot(O2oall,fanx','k-','linewidth',2)
    hold on;
    plot(O2oall,((D/yoe*(XQoall'+XQoPall') + D/yo_bno2*XQbno2all' + ...
        D/yo_bno2*XQbno2all')*1e3/400),'b-','LineWidth',2);
    plot(O2oall,((D/yn23e*(XQn23all'+XQn23Pall') + D*ena*XQaall')...
        *1e3/400),'r-','LineWidth',2);
    hold off
    title('e. Bulk properties ','FontSize',16);
    %set(gca,'yscale','log')
    set(gca,'xscale','log')
    set(gca,'xdir','reverse')
    %ylabel('cells/m^3 ','FontSize',14);grid on;
    ylabel('Fraction / nM d^-^1 ','FontSize',14);grid on;
    %legend('Aerobes ','n123 ','n1 ','n2 ','n3 ','n12 ','n23 ','Anammox ',12);
    legend('Fraction anammox','Total O_2 consumption','Total N loss');
    ylim([0 1]) %for rates
    %xlim([-2 2])
    xlabel('O_2 : OM supply (mol/mol)','fontsize',16)
    set(gca,'fontsize',16)
    set(gca,'ytick',[0,.5,1])
    set(gca,'yTickLabel',str2mat('0','0.5 / 200','1 / 400'))
    legend boxoff
    set(gca,'xtick',[1e-2,0.1,1,10,100])
    %set(gca,'xTickLabel',str2mat('1 nM','100 nM','10 \mu M'))
    grid on
    %%
    %save:
    
    %export_fig Figures/AllSix_Jul2018.png -m4
end