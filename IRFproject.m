%Graphs for a CMIP5 emission scearios corresponding to RCP's
%ppm and GtC is converted to GtCO2 (above steadystate) unless indicated by the name of the matrix
close all; clear;
cd(fileparts(mfilename('fullpath')))
%load data
[RCP,txt_RCP,Table_RCP] = xlsread('iamc_db_Emissions.xlsx',1,'E1:O8','basic'); % does not work without 'basic'
RCP(2:8,:)=RCP(2:8,:)/1000; % convert Mt into GtCO2
[J,txt_J,Table_J] = xlsread('decay and inertia parameters.xlsx',1,'A1:H19','basic');
[G,txt_G,Table_G] = xlsread('decay and inertia parameters.xlsx',3,'A2:G20','basic');
%[CONC] = xlsread('iamc_db_CO2ppm.xlsx',1,'F1:P8','basic');
[TEMP] = xlsread('iamc_db_Temp.xlsx',1,'F1:P8','basic');
[F_TOT] = xlsread('iamc_db_TotalF.xlsx',1,'F1:P8','basic');
[F_CO2] = xlsread('iamc_db_CO2F.xlsx',1,'F1:P8','basic');
F_NONCO2=F_TOT-F_CO2;
load('ALPHA_RCP.mat')
[IRFMethane]=xlsread('IRF_Methane.xlsx',1,'A2:B26','basic');
%[ForcingMethane]=xlsread('IRF_Methane.xlsx',1,'C2:D26','basic');%fill in the excell file in the correct cells
IRFMethane=interp1(IRFMethane(:,1),IRFMethane(:,2),[0:100],'makima');
IRFMethane=IRFMethane'; %is now a column vector
%ForcingMethane=interp1(ForcingMethane(:,1),ForcingMethane(:,2),[0:100],'makima');
%ForcingMethane=ForcingMethane'; %is now a column vector
%Parameters
T0atm=0.85;
T0ocean=0.28; % Initial Value from FAIR
MequilibriumGtCO2=588*3.664; %background GtCO2 stock in 2015 
Climatesensitivity=3.1; %value in DICE
Ecum2015=571*3.664; %Fair excel file has571. Richard Millar in Nature has 545GtC for 2015
firstyear=2015;
lastyear=2100;
year=firstyear:1:lastyear ;
Decennia=RCP(1,:);
Decennia2005=[2005 2010:10:2100];
%txt_RCP(9,1)={'Constant 2015 emissions'};RCP(9,:)=39.15;F_NONCO2(9,:)=0.01811;CONC(9,:)=missing;TEMP(9,:)=missing; 
%scenario1-7= offset starting in 2020 ending in 2050 with emissions background RCP1.9 to RCP8.5
%scenario8  = offset starting in 2020 ending in 2050  with constant emissions background
%scenario9  = increase of emission in 2020 and absorption in 2050 with RCP2.6 background
%scenario10= permanent increase emissions with RCP2.6 background
%%%%%%%%%%%%%%%%% KEY PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    projectlength=20; %(for projectlenght 80, comment out line 51 (no reemission)
    CarbonPerMethane=136; %this is the number of CO2 offsets for one tonne of methane emission.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Baseline emission path
    Epath=zeros(10,86); 
    for i=1:7 ;Epath(i,:)=interp1(Decennia,RCP(i+1,:),year, 'spline'); end
    %Epath(8,:) =39.15 ;  
    %Epath(9,:)=Epath(2,:);
    %Epath(10,:)=Epath(2,:);
    
    %emissions path with the project
    Eproject=10; % this is irrelevant, pulse is in this case 10GtCO2 to avoid too large rounding errors, scaled back to 1GtCO2 later

    Epath2=Epath;
    Epath2(1:8,6) =Epath2(1:8 ,6)-Eproject; %absorption in 2020
    Epath2(1:8,6+projectlength)=Epath2(1:8,6+projectlength)+Eproject;% reemission at end of project
    %Epath2(9:10,6)=Epath2(9:10,6)+Eproject; Epath2(9  ,56)=Epath2(9  ,56)-Eproject;
    %F_NONCO2(9,:)=0.01811;F_NONCO2(10,:)=F_NONCO2(3,:);F_NONCO2(11,:)=F_NONCO2(3,:);
    %TEMP(10,:)=TEMP(3,:);TEMP(11,:)=TEMP(3,:);
    %txt_RCP(9,1)={'Constant 2015 emissions'};
    %txt_RCP(10,:)={[txt_RCP{3,1} 'TemporaryIncrease']};
    %txt_RCP(11,:)={[txt_RCP{3,1} 'PermanentIncrease']};
Mi0=[0.473184713605214;0.359104481245797;0.143649591297395;0.0240612138515941]*ones(1,18)*263*3.664;% FAIR, calculated from historic emissions in file FitDiceDecay.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options= odeset('RelTol',1e-7);%1e-6 earlier version
for scenario=2:2% Scenario 2 is RCP2.6, change here for other background emissions
%JOOS GEOFFROY  
F_nonco2=interp1(Decennia2005,F_NONCO2(scenario+1,:),year);
F_nonco2_5y=interp1(Decennia2005,F_NONCO2(scenario+1,:),2015:5:2100);
M_Joos=zeros(18,86);M_Joos2=zeros(18,86);
for i=1:18
    [t,M]=ode45(@(t,M) odefcn_decay(t,M,Epath(scenario,:),year,J(i,:),1),[firstyear lastyear],Mi0(:,i),options);
    M=sum(transpose(M));
    Mtemporary=interp1(t,M(1,:),year);
    M_Joos(i,:)=transpose(Mtemporary);
    %with project
    [t,M2]=ode45(@(t,M) odefcn_decay(t,M,Epath2(scenario,:),year,J(i,:),1),[firstyear lastyear],Mi0(:,i),options);
    M2=sum(transpose(M2));
    Mtemporary2=interp1(t,M2(1,:),year);
    M_Joos2(i,:)=transpose(Mtemporary2);
end
M_Joos_diff=(M_Joos2-M_Joos)/Eproject;
T_JoosGeof=zeros(16*16,86);
T_JoosGeof2=zeros(16*16,86);
for k=1:16 %loop runs over carbon decay models
  for i=1:16 %loop runs over thermal inertia models
    Forcing=Climatesensitivity*G(i,3)*log((MequilibriumGtCO2+M_Joos(k,:))/MequilibriumGtCO2)/log(2)+F_nonco2;
    [t,T]=ode45(@(t,T) odefcn_inertia(t,T,Forcing,year,G(i,:)),[firstyear lastyear],[T0atm;T0ocean]);
    Ttemporary=interp1(t,T(:,1),year);
    T_JoosGeof(16*(k-1)+i,:)=transpose(Ttemporary);
    %with project
    Forcing2=Climatesensitivity*G(i,3)*log((MequilibriumGtCO2+M_Joos2(k,:))/MequilibriumGtCO2)/log(2)+F_nonco2;
    [t,T2]=ode45(@(t,T) odefcn_inertia(t,T,Forcing2,year,G(i,:)),[firstyear lastyear],[T0atm;T0ocean]);
    Ttemporary=interp1(t,T2(:,1),year);
    T_JoosGeof2(16*(k-1)+i,:)=transpose(Ttemporary);
  end 
end
T_JoosGeof_diff=(T_JoosGeof2-T_JoosGeof)/Eproject;
T_sorted=sort(T_JoosGeof_diff,1);
T_deciles=zeros(9,86);
for i=1:9
   T_deciles(i,:)=T_sorted(round(i*25.6),:);
end
%bestfit trajectory
Forcing=Climatesensitivity*G(17,3)*log((MequilibriumGtCO2+M_Joos(18,:))/MequilibriumGtCO2)/log(2)+F_nonco2;
[t,T]=ode45(@(t,T) odefcn_inertia(t,T,Forcing,year,G(17,:)),[firstyear lastyear],[T0atm;T0ocean],options);
Ttemporary=interp1(t,T(:,1),year);
T_JoosGeofbestfit(scenario,:)=transpose(Ttemporary);
%with project
Forcing2=Climatesensitivity*G(17,3)*log((MequilibriumGtCO2+M_Joos2(18,:))/MequilibriumGtCO2)/log(2)+F_nonco2;
[t,T2]=ode45(@(t,T) odefcn_inertia(t,T,Forcing2,year,G(17,:)),[firstyear lastyear],[T0atm;T0ocean],options);
Ttemporary=interp1(t,T2(:,1),year);
T_JoosGeofbestfit2(scenario,:)=transpose(Ttemporary);
T_JoosGeofbestfit_diff(scenario,:)=(T_JoosGeofbestfit2(scenario,:)-T_JoosGeofbestfit(scenario,:))/Eproject;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FAIR
[t,X]=ode45(@(t,X) odefcn_endogAlpha(t,X,year,Epath(scenario,:),F_nonco2,J(18,:),G(17,:),Climatesensitivity),...
    [firstyear lastyear],[Mi0(:,1);T0atm;T0ocean;Ecum2015],options);
Ttemporary=interp1(t,X(:,5),year);
T_FAIR(scenario,:)=transpose(Ttemporary);
Mtemporary=sum(X(:,1:4),2);
M_FAIR(scenario,:)=transpose(interp1(t,Mtemporary,year));
F_FAIR(scenario,:)=Climatesensitivity*G(17,3)*log(1+M_FAIR(scenario,:)/(3.663*588))/log(2);
%with project
[t,X2]=ode45(@(t,X) odefcn_endogAlpha(t,X,year,Epath2(scenario,:),F_nonco2,J(18,:),G(17,:),Climatesensitivity),...
    [firstyear lastyear],[Mi0(:,1);T0atm;T0ocean;Ecum2015],options);
Ttemporary=interp1(t,X2(:,5),year);
T_FAIR2(scenario,:)=transpose(Ttemporary);
Mtemporary2=sum(X2(:,1:4),2);
M_FAIR2(scenario,:)=transpose(interp1(t,Mtemporary2,year));
F_FAIR2(scenario,:)=Climatesensitivity*G(17,3)*log(1+M_FAIR2(scenario,:)/(3.663*588))/log(2);
M_FAIR_diff(scenario,:)=(M_FAIR2(scenario,:)-M_FAIR(scenario,:))/Eproject;
T_FAIR_diff(scenario,:)=(T_FAIR2(scenario,:)-T_FAIR(scenario,:))/Eproject;
F_FAIR_diff(scenario,:)=(F_FAIR2(scenario,:)-F_FAIR(scenario,:))/Eproject;
GWP=sum(F_FAIR_diff(scenario,:)) %report this manually to fill in in the equivalence Table
figure()
plot(year,F_FAIR_diff(scenario,:))

%--------------------------------------------------------------------------------------------------------------------------------
%scale offset impact with CarbonPerMethane
T_FAIR_offset=CarbonPerMethane.*T_FAIR_diff(scenario,:);
T_deciles_offset=CarbonPerMethane.*T_deciles;
T_JoosGeof_offset=CarbonPerMethane*T_JoosGeofbestfit_diff(scenario,:);
T_Methane=[zeros(1,5) IRFMethane(1:81)'];
%Damages    
kappa=0.0077; % output loss for 1°C
gamma=2*kappa; %total damages=exp(-gamma/2^*T²)
Damage_FAIR=gamma.*T_FAIR(scenario,:).*T_FAIR_offset;
Damage_Methane=gamma.*T_FAIR(scenario,:).*T_Methane;
discount=0.03-0.02; %r-g to account for damages being proportional to production
PV_Damage_FAIR=Damage_FAIR.*exp(-discount*(0:85));
PV_Damage_Methane=Damage_Methane.*exp(-discount*(0:85));

%PLOTS
%subplot(3,1,1);
figure()
    ar=area(transpose(year),transpose(T_deciles_offset(9,:)), 'FaceColor' , [0.5 0.98 0.6],'basevalue',-1);
    hold on
    area(transpose(year),transpose(T_deciles_offset(8,:)), 'FaceColor' , [0.5 0.96 0.6],'basevalue',-1)
    area(transpose(year),transpose(T_deciles_offset(7,:)), 'FaceColor' , [0.5 0.90 0.6],'basevalue',-1)
    area(transpose(year),transpose(T_deciles_offset(6,:)), 'FaceColor' , [0.5 0.80 0.6],'basevalue',-1)
    area(transpose(year),transpose(T_deciles_offset(4,:)), 'FaceColor' , [0.5 0.90 0.6],'basevalue',-1)
    area(transpose(year),transpose(T_deciles_offset(3,:)), 'FaceColor' , [0.5 0.96 0.6],'basevalue',-1)
    area(transpose(year),transpose(T_deciles_offset(2,:)), 'FaceColor' , [0.5 0.98 0.6],'basevalue',-1)
    area(transpose(year),transpose(T_deciles_offset(1,:)), 'FaceColor' , [1 1 1],'basevalue',-1)
    p1=plot(year,T_JoosGeof_offset,'color',[0 0.2 0],'linewidth',1.5);
    p2=plot(year, T_FAIR_offset,'color', [1 0.8 0],'linewidth',1);
    %p3=plot(year, [zeros(1,5+delay) mean(T_FAIR_diff(scenario,5+delay:5+projectlength+delay))*ones(1,projectlength) zeros(1,81-projectlength-delay)],'linewidth',1);
    p3=plot(year, T_Methane,'color', [0 0 1],'linewidth',1);
    p4=plot(year, T_Methane+T_FAIR_offset,'color', [1 0 0],'linewidth',1);
    yline(0,'linewidth',0.1)
    hold off
    lgd=legend([ar p1 p2 p3 p4],'Deciles Joos-Geoffroy combinations','Best fit Joos-Geoffroy ensemble','FAIR','Methane','Net Effect','Location','SouthEast');
    lgd.FontSize=8;
    xlabel('years')
    ylabel('Effect on Temperature (°C)')
    axis([2015 2100 -inf inf]);
    pbaspect([3 1 1])
    orient landscape
    print(['MethaneOffsetTemp',num2str(projectlength),'y'],'-fillpage', '-dpdf' ); 

%subplot(3,1,2);%damages
figure()
    p2=plot(year, Damage_FAIR,'color', [1 0.8 0],'linewidth',1);
    hold on
    p3=plot(year, Damage_Methane,'color', [0 0 1],'linewidth',1);
    p4=plot(year, Damage_Methane+Damage_FAIR,'color', [1 0 0],'linewidth',1);
    yline(0,'linewidth',0.1)
    hold off
    lgd=legend([p2 p3 p4],'Carbon Offsets','Methane','Net Effect','Location','SouthEast');
    lgd.FontSize=8;
    xlabel('years')
    ylabel('Damages (relative to GDP)')
    set(gca,'yticklabel',[])
    axis([2015 2100 -inf inf]);
    pbaspect([3 1 1])
    orient landscape
    print(['MethaneOffsetDamage',num2str(projectlength),'y'],'-fillpage', '-dpdf' ); 
%subplot(3,1,3);  % Present value damages  
figure()
    p2=plot(year, PV_Damage_FAIR,'color', [1 0.8 0],'linewidth',1);
    hold on
    p3=plot(year, PV_Damage_Methane,'color', [0 0 1],'linewidth',1);
    p4=plot(year, PV_Damage_Methane+PV_Damage_FAIR,'color', [1 0 0],'linewidth',1);
    yline(0,'linewidth',0.1)
    hold off
    lgd=legend([p2 p3 p4],'Carbon Offsets','Methane','Net Effect','Location','SouthEast');
    lgd.FontSize=8;
    xlabel('years')
    ylabel('PV of Damages')
    set(gca,'yticklabel',[])
    axis([2015 2100 -inf inf]);
    pbaspect([3 1 1])
    orient landscape
    print(['MethaneOffsetDamagePV',num2str(projectlength),'y'],'-fillpage', '-dpdf' ); 

end
