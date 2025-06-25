%Graphs for a FAIR2 emission scearios corresponding to RCP's
%CO2 is in ppm and GtC (=PgC) and methane is in ppb and MtCH4 (=TgCH4)
close all; clear;
cd(fileparts(mfilename('fullpath')));

%parameters
kappa=0.0077; % output loss for 1°C, damage parameter from Howard & Sterner (fit on 3°C damage), 
gamma=2*kappa; %total damages=exp(-gamma/2^*T²)
discount=0.025; %discount rate
g=0.02; %growth rate of the economy
GDP0=106000; %Billion$=G$ in 2023 (world bank)
year_n=81;%number of years on graph
year=2020:1:2020+year_n-1 ; %periods for graph
decennia=2020:10:2100; %dates of RCP emission data.
horizon=500; %number of years to calculate SCC and SCM
sc=2; % Scenario 2 is RCP2.6, change here for graphs with other background emissions

%load parameters
paramR = readmatrix('Complete_gas_cycle_params.csv', 'Range', 'B3:CD19');
paramR=paramR(:,[9,37,41])'; %CO2,CH4,N20
paramT=[0.903, 7.92, 355, 0.180, 0.297, 0.386];

%load future RCP paths from 2020:2100
txt_RCP=readcell('iamc_db_Emissions.xlsx','Range','E2:E8'); 
CO2=readmatrix('iamc_db_Emissions.xlsx','Range','G2:O8');
CO2=CO2/1000/3.66; % convert MtCO2 into GtC
CH4=readmatrix('iamc_db_Methane.xlsx','Range','G2:O8'); %already in MtCH4 
TEMP = readmatrix('iamc_db_Temp.xlsx','Range','H2:P8');
F_TOT = readmatrix('iamc_db_TotalF.xlsx','Range','H2:P8');
F_CO2 = readmatrix('iamc_db_CO2F.xlsx','Range','H2:P8');

%load historic emissions
E_historic=readmatrix('ghg-emissions.csv','Range','C2:FT5'); %PIK emissions from 1850 until 2023 downloaded from climatewatchdata.org 
E_historic=E_historic([2,1,3],:); % set CO2 on line 1, Methane on line 2 (NO2 on line 3)
E_historic(1,:)=E_historic(1,:)/1000/3.66; % convert from MtCO2 to PgC (or GtC)
E_historic(2,:)=E_historic(2,:)/21; %convert from MtCO2eq to TgCH4 (or MtCH4)
E_historic(3,:)=E_historic(3,:)/273; %convert from MtCO2eq to TgN20 (or MtN20)

%find initial values in 2023
F_ext=zeros(1,size(E_historic,2));
R0=zeros(5,3);
S0=zeros(3,1);
[T,F,S,R]=FAIR2(E_historic, F_ext, R0, S0, paramR, paramT);
R0=squeeze(R(1:5,end,1:2));
S0=S(:,end);

%Baseline emission and forcing paths
E_base_CO2=zeros(7,year_n); 
E_base_CH4=zeros(7,year_n); 
F_ext_CO2=zeros(7,year_n); 
F_ext_TOT=zeros(7,year_n); 
%interpolate between decennia
for i=1:7 
    E_base_CO2(i,:)=interp1(decennia,CO2(i,:),year, 'spline'); 
    E_base_CH4(i,:)=interp1(decennia,CH4(i,:),year, 'spline'); 
    F_ext_CO2(i,:)=interp1(decennia,F_CO2(i,:),year, 'spline'); 
    F_ext_TOT(i,:)=interp1(decennia,F_TOT(i,:),year, 'spline'); 
end
%extend to horizon
E_base_CO2(:,82:horizon)=E_base_CO2(:,year_n)*ones(1,horizon-year_n);
E_base_CH4(:,82:horizon)=E_base_CH4(:,year_n)*ones(1,horizon-year_n);
F_ext_CO2(:,82:horizon) = F_ext_CO2(:,year_n)*ones(1,horizon-year_n);
F_ext_TOT(:,82:horizon) = F_ext_TOT(:,year_n)*ones(1,horizon-year_n);

F_ext_NONCO2=F_ext_TOT - F_ext_CO2;

%LOOP TO MAKE FIGURES
Projectlength=[20,30,40,100];
Equiv_distr=zeros(length(Projectlength),5);
for i=1:length(Projectlength)

    %emissions path with the project
    %CHANGE HERE FOR A PROJECT WITH GRADUAL REMOVAL
    projectstart=2;%absorption in 2021
    Projectsize=1/3.66; % project of  1GtC02 
    Shock=zeros(1,horizon);
    Shock(:,projectstart)=-Projectsize;  
    Shock(:,projectstart+Projectlength(i))=+Projectsize;
    E_project_CO2=E_base_CO2+Shock;
    
    E_project_CH4=E_base_CH4;
    E_project_CH4(:,projectstart)=E_project_CH4(:,projectstart) + 1000; %emission of 1GtCH4
     
    %IRF methane
    [T_IRF_CH4, F_IRF_CH4, Damage_IRF_CH4,PV_Damage_IRF_CH4, GWP_CH4, SCM]=IRF_FAIR2(E_base_CH4(sc,:), E_project_CH4(sc,:), F_ext_TOT*0.9, R0(:,2), S0, paramR(2,:), paramT, gamma, discount,g);
    %IRF carbon
    [T_IRF_CO2, F_IRF_CO2, Damage_IRF_CO2,PV_Damage_IRF_CO2, GWP_CO2, SVO]=IRF_FAIR2(E_base_CO2(sc,:), E_project_CO2(sc,:), F_ext_NONCO2, R0(:,1), S0, paramR(1,:), paramT, gamma, discount,g);
    
    %scale the carbon project to make it welfare neutral
    Equivalence=-SCM/SVO; %tonnes of CO2 offsets per tCH4
    T_IRF_CO2 = T_IRF_CO2 * Equivalence;
    Damage_IRF_CO2 = Damage_IRF_CO2 * Equivalence;
    PV_Damage_IRF_CO2 = PV_Damage_IRF_CO2 * Equivalence;
    Damage_net = Damage_IRF_CO2 + Damage_IRF_CH4;
    PV_Damage_net = PV_Damage_IRF_CO2 + PV_Damage_IRF_CH4;
    
    
    %MONTE CARLO
    % pct = 0.1:0.1:0.9;% in case we want deciles for only 1 stochstic parameter.
    % z_CH4f = norminv(pct, 1, 0.170); %variance from Table 6, "Effective Radiative Forcing".
    % z_CH4tau= norminv(pct,1,0.1); %not available for methane but e-folding lifetime is 9.1+-0.9 years in WGI chapter 6 Table 6.2
    % z_CO2f = norminv(pct, 1, 0.122); % variance from Table 6, "Effective Radiative Forcing".
    % z_CO2r0 = norminv(pct, 1, 0.154); %st dev from Table 5
    % z_CO2ru = logninv(pct, 0, 0.44);
    % z_CO2rT = norminv(pct, 1, 0.615);
    
    n=5000; % number of Monte Carlo runs
    z = 1 + randn(n,6).*[0.170,0.1,0.122,0.154,0 ,0.615];
    z(:,5)=lognrnd(0,0.44,n,1);
    z(z<0)=0.001; %truncate in case there would be neg values
    %preallocate
    T_MC_CO2=zeros(n,horizon);
    T_MC_CH4=zeros(n,horizon);
    SCM_MC=zeros(n,1);
    SVO_MC=zeros(n,1);
    paramR_MC=paramR;
    for j=1:n
        %IRF methane
        paramR_MC(2,15:17)=paramR(2,15:17) .* z(j,1) ; % z_CH4f(j)
        paramR_MC(2,4)=paramR(2,4)*z(j,2);
        [T_MC_CH4(j,:), ~, ~,~, ~, SCM_MC(j)]=IRF_FAIR2(E_base_CH4(sc,:), E_project_CH4(sc,:), F_ext_TOT*0.9, R0(:,2), S0, paramR_MC(2,:), paramT, gamma, discount,g);
        %IRF carbon
        paramR_MC(1,15:17)=paramR(1,15:17) * z(j,3); %   z_CO2f(j)
        paramR_MC(1,9:11)=paramR(1,9:11) .* z(j,4:6);     
        [T_MC_CO2(j,:), ~, ~,~, ~, SVO_MC(j)]=IRF_FAIR2(E_base_CO2(sc,:), E_project_CO2(sc,:), F_ext_NONCO2, R0(:,1), S0, paramR_MC(1,:), paramT, gamma, discount,g);
    end
    
    T_MC_CO2=T_MC_CO2* Equivalence;
    T_MC_net=T_MC_CH4+T_MC_CO2;
    %probability distribution of equivalence
    Equiv_MC=-SCM_MC./SVO_MC;
    figure();
    histogram(Equiv_MC,50);
    %print(['EquivalenceHistogram',num2str(Projectlength(i)),'y'],'-fillpage','-dpdf' );  % use 5000 runs to create smooth histogram.
    writematrix(Equiv_MC, ['Equivalence_MonteCarlo_',num2str(Projectlength(i)),'.csv']);
    exportgraphics(gcf,['EquivalenceHistogram',num2str(Projectlength(i)),'y.png'], 'Resolution', 300);
    Equiv_distr(i,1)=Projectlength(i);
    Equiv_distr(i,2)=Equivalence;
    Equiv_distr(i,3)=median(Equiv_MC);
    Equiv_distr(i,4)=mean(Equiv_MC);
    Equiv_distr(i,5)=std(Equiv_MC);
    
    %Sort to extract percentiles
    T_MC_CO2=sort(T_MC_CO2,1);
    T_MC_CH4=sort(T_MC_CH4,1);
    T_MC_net=sort(T_MC_net,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FIGURES
    %figure temperature
    figure();    hold on;
        ar1=area(transpose(year),transpose(T_MC_CH4(n*.9,1:year_n)), 'FaceColor' , [0.9 0.9 1],'basevalue',-1);% [0.9 1 0.9]
        ar1.EdgeColor= 'none';
        area(transpose(year),transpose(T_MC_CH4(n*.8,1:year_n)), 'FaceColor' , [0.7 0.7 0.96],'basevalue',-1);% [0.7 0.96 0.7]
        area(transpose(year),transpose(T_MC_CH4(n*.7,1:year_n)), 'FaceColor' , [0.7 0.7 0.90],'basevalue',-1);% [0.7 0.90 0.7]
        area(transpose(year),transpose(T_MC_CH4(n*.6,1:year_n)), 'FaceColor' , [0.7 0.7 0.80],'basevalue',-1);% [0.7 0.80 0.7]
        area(transpose(year),transpose(T_MC_CH4(n*.4,1:year_n)), 'FaceColor' , [0.7 0.7 0.90],'basevalue',-1);% [0.7 0.90 0.7]
        area(transpose(year),transpose(T_MC_CH4(n*.3,1:year_n)), 'FaceColor' , [0.7 0.7 0.96],'basevalue',-1);% [0.7 0.96 0.7]
        area(transpose(year),transpose(T_MC_CH4(n*.2,1:year_n)), 'FaceColor' , [0.9 0.9 1],'basevalue',-1);% [0.9 1 0.9]
        ar2=area(transpose(year),transpose(T_MC_CH4(n*.1,1:year_n)), 'FaceColor' , [1 1 1],'basevalue',-1);% 
        ar2.EdgeColor= 'none';
    
        ar3=area(transpose(year),transpose(T_MC_net(n*.9,1:year_n)), 'FaceColor' , [1 0.9 0.9 ],'basevalue',-1,'FaceAlpha',0.3);
        ar3.EdgeColor= 'none';
        area(transpose(year),transpose(T_MC_net(n*.8,1:year_n)), 'FaceColor' , [0.96 0.7 0.7 ],'basevalue',-1,'FaceAlpha',0.3)
        area(transpose(year),transpose(T_MC_net(n*.7,1:year_n)), 'FaceColor' , [0.90 0.7 0.7 ],'basevalue',-1,'FaceAlpha',0.3)
        area(transpose(year),transpose(T_MC_net(n*.6,1:year_n)), 'FaceColor' , [0.80 0.7 0.7 ],'basevalue',-1)
        area(transpose(year),transpose(T_MC_net(n*.4,1:year_n)), 'FaceColor' , [0.90 0.8 0.8 ],'basevalue',-1)
        area(transpose(year),transpose(T_MC_net(n*.3,1:year_n)), 'FaceColor' , [0.96 0.8 0.8 ],'basevalue',-1)
        area(transpose(year),transpose(T_MC_net(n*.2,1:year_n)), 'FaceColor' , [1 0.9 0.9 ],'basevalue',-1)
        ar4=area(transpose(year),transpose(T_MC_net(n*.1,1:year_n)), 'FaceColor' , [1 1 1],'basevalue',-1);
        ar4.EdgeColor= 'none';
    
        ar5=area(transpose(year),transpose(T_MC_CO2(n*.9,1:year_n)), 'FaceColor' , [0.9 1 0.9],'basevalue',-1);
        ar5.EdgeColor= 'none';
        area(transpose(year),transpose(T_MC_CO2(n*.8,1:year_n)), 'FaceColor' , [0.7 0.96 0.7],'basevalue',-1)
        area(transpose(year),transpose(T_MC_CO2(n*.7,1:year_n)), 'FaceColor' , [0.7 0.90 0.7],'basevalue',-1)
        area(transpose(year),transpose(T_MC_CO2(n*.6,1:year_n)), 'FaceColor' , [0.7 0.80 0.7],'basevalue',-1)
        area(transpose(year),transpose(T_MC_CO2(n*.4,1:year_n)), 'FaceColor' , [0.7 0.90 0.7],'basevalue',-1)
        area(transpose(year),transpose(T_MC_CO2(n*.3,1:year_n)), 'FaceColor' , [0.7 0.96 0.7],'basevalue',-1)
        area(transpose(year),transpose(T_MC_CO2(n*.2,1:year_n)), 'FaceColor' , [0.9 1 0.9],'basevalue',-1)
        ar6=area(transpose(year),transpose(T_MC_CO2(n*.1,1:year_n)), 'FaceColor' , [1 1 1],'basevalue',-1);
        ar6.EdgeColor= 'none';
        p2=plot(year, T_IRF_CO2(1,1:year_n) ,'color', [0.3, 0.5, 0.3],'linewidth',1);%[1, 0.8, 0]
        p3=plot(year, T_IRF_CH4(1,1:year_n),'color', [0 0 1],'linewidth',1);
        p4=plot(year, T_IRF_CO2(1,1:year_n)  + T_IRF_CH4(1,1:year_n),'color', [1 0 0],'linewidth',1);
        yline(0,'linewidth',0.1)
        hold off
        if Projectlength(i) == 20 || Projectlength(i) == 30; lgd=legend([p3 p2 p4 ],'Methane Emission','CO2 Removal','Net Effect','Location','SouthEast');end
        lgd.FontSize=8;
        %xlabel('years')
        ylabel('Temperature effect (°C)', 'FontSize', 8);
        axis([2015 2100 -inf inf]);
        set(gca, 'FontSize', 8);
        pbaspect([3 1 1])
        orient landscape
        %print(['MethaneOffsetTemp',num2str(Projectlength(i)),'y'],'-fillpage', '-dpdf' ); 
        exportgraphics(gcf, ['MethaneOffsetTemp',num2str(Projectlength(i)),'y.png'], 'Resolution', 300);
    
    %figure damages
    figure() 
        p2=plot(year, GDP0*Damage_IRF_CO2(1,1:year_n),'color', [1 0.8 0],'linewidth',1);
        hold on
        p3=plot(year, GDP0*Damage_IRF_CH4(1,1:year_n),'color', [0 0 1],'linewidth',1);
        p4=plot(year, GDP0*Damage_net(1,1:year_n),'color', [1 0 0],'linewidth',1);
        yline(0,'linewidth',0.1)
        hold off
        if Projectlength(i) == 20 || Projectlength(i) == 30; lgd=legend([p2 p3 p4],'Carbon Offsets','Methane','Net Effect','Location','SouthEast');end
        lgd.FontSize=8;
        xlabel('years')
        ylabel('Damages ($)', 'FontSize', 8)
        set(gca,'yticklabel',[])
        axis([2015 2100 -inf inf]);
        set(gca, 'FontSize', 8);
        pbaspect([3 1 1])
        orient landscape
        %print(['MethaneOffsetDamage',num2str(Projectlength(i)),'y'],'-fillpage', '-dpdf' ); 
        exportgraphics(gcf, ['MethaneOffsetDamage',num2str(Projectlength(i)),'y.png'], 'Resolution', 300);
    
    % Present value damages  
    figure()
        p2=plot(year, GDP0*PV_Damage_IRF_CO2(1,1:year_n),'color', [1 0.8 0],'linewidth',1);
        hold on
        p3=plot(year, GDP0*PV_Damage_IRF_CH4(1,1:year_n),'color', [0 0 1],'linewidth',1);
        p4=plot(year, GDP0*PV_Damage_net(1,1:year_n),'color', [1 0 0],'linewidth',1);
        yline(0,'linewidth',0.1)
        hold off
        if Projectlength(i) == 20 || Projectlength(i) == 30; lgd=legend([p2 p3 p4],'Carbon Offsets','Methane','Net Effect','Location','SouthEast');end
        lgd.FontSize=8;
        xlabel('years')
        ylabel('PV of Damages ($)', 'FontSize', 8)
        %set(gca,'yticklabel',[])
        axis([2015 2100 -inf inf]);
        set(gca, 'FontSize', 8);
        pbaspect([3 1 1])
        orient landscape
        %print(['MethaneOffsetDamagePV',num2str(Projectlength(i)),'y'],'-fillpage', '-dpdf' ); 
        exportgraphics(gcf, ['MethaneOffsetDamagePV',num2str(Projectlength(i)),'y.png'], 'Resolution', 300);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TABLE
sc=[2,3,6];%RCP2.6, 4.5 and 6
discount=[0.025 0.03 0.035];%the discount rate
phi=[0 0.005 0.01];%failure risk
Projectlength=[20 25 30 35 40 100 498];%
Table_eq=zeros(28,10);
Table_net=zeros(28,10);

for j=1:3 %RCP
    for k=1:3 %discount rate
        [~, ~, ~,PV_Damage_IRF_CH4, ~, SCM]=IRF_FAIR2(E_base_CH4(sc(j),:), E_project_CH4(sc(j),:), F_ext_TOT(sc(j),:)*0.9, R0(:,2), S0, paramR(2,:), paramT, gamma, discount(k),g); 
        if j==1;[~, ~, ~,~, GWP_CH4,~]     =IRF_FAIR2(E_base_CH4(sc(j),1:100), E_project_CH4(sc(j),1:100), F_ext_TOT(sc(j),1:100)*0.9, R0(:,2), S0, paramR(2,:), paramT, gamma, discount(k),g);end
        for l=1:3 %phi+varphi
            for m=1:7 %endtime
                Shock=zeros(1,horizon);
                Shock(:,projectstart)=-Projectsize;  
                Shock(:,projectstart+Projectlength(m))=Projectsize;
                E_project_CO2=E_base_CO2+Shock;
                [~, ~, ~,PV_Damage_IRF_CO2, GWP_CO2, SVO]=IRF_FAIR2(E_base_CO2(sc(j),:), E_project_CO2(sc(j),:), F_ext_NONCO2(sc(j),:), R0(:,1), S0, paramR(1,:), paramT, gamma, discount(k)+phi(l),g);
                %Equivalence table
                Equivalence=-SCM/SVO;%nature paper has inverse SCO./SCCorM
                Table_eq(9*(j-1)+3*(k-1)+l,1)=sc(j);
                Table_eq(9*(j-1)+3*(k-1)+l,2)=discount(k);
                Table_eq(9*(j-1)+3*(k-1)+l,3)=phi(l);
                Table_eq(9*(j-1)+3*(k-1)+l,3+m)=Equivalence; %column 4,to 10 %Table_eq(9*(j-1)+3*(k-1)+l,11)=-SVO*GDP0; %column 11 SCC (at max project length)
                Table_eq(28,3+m)=-GWP_CH4/GWP_CO2; % last row column 4 to 10
                %Table with net damages.
                PV_Damage_net=PV_Damage_IRF_CO2*Equivalence + PV_Damage_IRF_CH4;
                Table_net(9*(j-1)+3*(k-1)+l,1)=sc(j);
                Table_net(9*(j-1)+3*(k-1)+l,2)=discount(k);
                Table_net(9*(j-1)+3*(k-1)+l,3)=phi(l);
                Table_net(9*(j-1)+3*(k-1)+l,3+m)=sum(abs(PV_Damage_net))/SCM; 
            end
        end
    end
end
[~,idx]=min(Table_net(:,4:10),[],2);
Table_minNetDamage=zeros(27,1);
for i=1:27; Table_minNetDamage(i)=(Projectlength(idx(i)));end
writematrix(Table_eq, 'Table_Equivalence.xlsx')

writematrix(Table_net, 'Table_NetDamages.xlsx')

