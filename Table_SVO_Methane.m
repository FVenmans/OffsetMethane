%Make table on adjustment factors
clear; close all;
%constant parameters
    g=0.02;
    %eta=1.35;
    %delta=0.005;
    %discount=delta+eta*g;
    discount=[0.025 0.03 0.035];
    GDP0=85000;
    T0=1.2; %°C
    zeta=0.0006;
    tau1=3;
%parameters which will vary
    Endtime=[20 40 395];
    tau2=Endtime+tau1;
    phivarphi=[0 0.005 0.01];
    varphi2=1000;
    kappa=0.0077; % output loss for 1°C
    gamma=2*kappa; %total damages=exp(-gamma/2^*T²)
%Load Temperature path and interpolate.
    [TEMP] = xlsread('iamc_db_Temp.xlsx',1,'H1:P8','basic');%starts in 2020
    TEMP=TEMP([3 4 7 8],:); %select scenario SSP1-26, SSP2-45, SSP4-60 and SSP5-85(Baseline)
    for i=1:4;    Tpath(i,:)=interp1([0:10:80],TEMP(i,:),[0:80], 'spline');end
    Tpath(:,[82:401])=Tpath(:,81)*ones(1,320); %constant temperature after 2100
    Tpath=Tpath';
%Load methane and CO2 IRF and interpolate
    [IRFMethane]=xlsread('IRF_Methane.xlsx',1,'A2:B26','basic');
    IRFMethane=interp1(IRFMethane(:,1),IRFMethane(:,2),[0:100],'spline');
    IRFMethane=IRFMethane';
    %IRFmethane(:,[102:401])=Tpath(:,81)*ones(1,320); %constant temperature after 2100
    [IRFCO2]=xlsread('IRF_CO2.xlsx',1,'A2:B19','basic');
    IRFCO2=interp1(IRFCO2(:,1),IRFCO2(:,2),[0:100],'spline');
    zeta=mean([IRFCO2])
%Paths
    year=[0:400]';
    GDP=GDP0*(exp(g*year));
    Discount=exp(-year*discount);
    Table=zeros(36,11);



for j=1:3 %RCP
    for k=1:3 %discount rate
        for l=1:3 %phi+varphi
            for m=1:3 %lifetime
                for i=1:2 %SCC and SCM
    RiskFactor=exp(-phivarphi(l)*year).*(1-exp(-varphi2.*year));
    if i==1; SCCorM=sum(gamma.*IRFMethane.*Tpath(1:101,j).*GDP(1:101).*Discount(1:101,k)); %Social cost of methane
    else;    SCCorM=sum(gamma.*zeta.*Tpath(tau1:end,j).*GDP(tau1:end).*Discount(tau1:end,k));end  %SCC=gamma*zeta²*S=gamma*zeta*Temp
    SCO=sum(gamma.*zeta.*Tpath(tau1:tau2(m),j).*GDP(tau1:tau2(m)).*Discount(tau1:tau2(m),k).*RiskFactor(tau1:tau2(m)));
    AdjFactor=transpose(SCCorM./SCO);%nature paper has inverse SCO./SCCorM
    Table(9*(j-1)+3*(k-1)+l,1)=j;
    Table(9*(j-1)+3*(k-1)+l,2)=discount(k);
    Table(9*(j-1)+3*(k-1)+l,3)=phivarphi(l);
    Table(9*(j-1)+3*(k-1)+l,3+(i-1)*4+m)=AdjFactor; %column 4,5,6 and 8,9,10
    Table(9*(j-1)+3*(k-1)+l,7+4*(i-1))=SCCorM;%column 7 and 11 
                end
            end
        end
    end
end

%Cubic damages
gamma=gamma/2;
TableCubic=zeros(36,11);
for j=1:3 %RCP
    for k=1:3 %discount rate
        for l=1:3 %phi+varphi
            for m=1:3 %lifetime
                for i=1:2 %SCC and SCM
    RiskFactor=exp(-phivarphi(l)*year).*(1-exp(-varphi2.*year));
    if i==1; SCCorM=sum(gamma.*IRFMethane.*Tpath(1:101,j).^2.*GDP(1:101).*Discount(1:101,k)); %Social cost of methane
    else;    SCCorM=sum(gamma.*zeta.*Tpath(tau1:end,j).^2.*GDP(tau1:end).*Discount(tau1:end,k));end  %SCC=gamma*zeta²*S=gamma*zeta*Temp
    SCO=sum(gamma.*zeta.*Tpath(tau1:tau2(m),j).^2.*GDP(tau1:tau2(m)).*Discount(tau1:tau2(m),k).*RiskFactor(tau1:tau2(m)));
    AdjFactor=transpose(SCCorM./SCO);%nature paper has inverse SCO./SCCorM
    TableCubic(9*(j-1)+3*(k-1)+l,1)=j;
    TableCubic(9*(j-1)+3*(k-1)+l,2)=discount(k);
    TableCubic(9*(j-1)+3*(k-1)+l,3)=phivarphi(l);
    TableCubic(9*(j-1)+3*(k-1)+l,3+(i-1)*4+m)=AdjFactor; %column 4,5,6 and 8,9,10
    TableCubic(9*(j-1)+3*(k-1)+l,7+4*(i-1))=SCCorM;%column 7 and 11 
                end
            end
        end
    end
end