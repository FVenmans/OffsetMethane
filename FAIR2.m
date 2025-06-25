%FAIR 2.0.0
%OUTPUTS
%T is atmospheric warming (in °C)
%F is total forcing
%S is a 3xt matrix with the 3 temperature boxes (adding up to the atmospheric temperature) over t years 
%R is an 7xtxp matrix related to the gaz absorption cycle, t years, p gases (if only one gas is modelled, 7xt matrix works too). The 7 outputs are
    % Reservoir1 in PgC or TgCH4
    % Reservoir2
    % Reservoir3
    % Reservoir4
    % Cumulative emissions in PgC or TgCH4
    % Alpha
    % Forcing of the particular gaz
%INPUTS
%E: pxt matrix, emissions in PgC or TgCH4
%F_ext: pxt matrix, other forcing
%R0 is 5xp (initial alpha and forcing will be calculated)
%S0 is 3x1
%paramR: px17 matrix, table S2= a1:4, tau1:4,r0,ru,rT,ra,PI_conc,emis2conc,f1:3
%paramT: 1x6 matrix, d1,d2,d3,q1,q2,q3, (standard values 0.903, 7.92, 355, 0.180, 0.297, 0.386)
%units: CO2 ppm and PgC, methane ppb and TgCH4 

function [T,F,S,R]=FAIR2(E, F_ext, R0, S0, paramR, paramT)
%number of years t and gases p
t_end=size(E,2); 
p_end=size(E,1); 


%preallocate output matrices
T=zeros(1,t_end);
F=zeros(1,t_end);
S=zeros(3,t_end);
R=zeros(7,t_end,p_end); 

%inititial values
T(:,1)=sum(S0);
S(:,1)=S0;
for p=1:p_end
    R(1:size(R0,1),1,p)=R0(:,p);
end

%F2xCO2= paramR(15) * log(2) + paramR(17) * (sqrt(2*paramR(13))-sqrt(paramR(13)));

%parameters g1 and g0 for alpha
g1=zeros(1,p_end);
g0=zeros(1,p_end);
iIRF=zeros(1,p_end);
for p=1:p_end
    for i=1:4
        g1(p)=g1(p) + paramR(p,i) * paramR(p,i+4) * (1 - (1+100/paramR(p,i+4)) * exp(-100/paramR(p,i+4)));
        iIRF(p)=iIRF(p) + paramR(p,i) * paramR(p,i+4) * (1-exp(-100/paramR(p,i+4))) ;
    end 
    g0(p)=exp(-iIRF(p)/g1(p)); %These values do not give alpha=1 in 2015, but rather 0.41 in 2025, they do give alpha=1 for methane and N2O
end

for t=2:t_end %for each period of 1 year
    for p=1:p_end %for each gas
        
        %calculate alpha for preceding period
        Ga=sum(R(1:4,t-1,p),1); %total atmospheric burden since preindustrial era in GtC
        Gu=R(5,t-1,p)-Ga; %total uptake of gas in sinks (or decayed) in GtC        
        R(6,t,p)= g0(p) * exp((paramR(p,9) + paramR(p,10)*Gu + paramR(p,11)*T(t-1) + paramR(p,12)*Ga)/g1(p)); %alpha at t is based on values from t-1 in python code (SM9)
        
        %calculate carbon (or other gas) in each of the 4 reservoirs (or methane in 1 reservoir)
        for k=1:4 
            R(k,t,p)= paramR(p,k)*E(p,t) + R(k,t-1,p)*(1-1/paramR(p,k+4)/R(6,t,p));  %approximation Rt=a_i*E + R_t-1*(1-1/tau_i/alpha) reservoir 1 to 4 in GtC
            %R(k,t,p)= E(p,t) * paramR(p,k) * R(6,t,p) * paramR(p,k+4) *(1-exp(1/R(6,t,p) / paramR(p,k+4)))  + R(k,t-1,p) * exp(1/R(6,t,p)/paramR(p,k+4)); %precise equation from SM9 gives negative values, probably wrong.
        end
        R(k,t,:)=R(k,t,:).*(R(k,t,:)>0); % sometimes alpha is below 0.25 which leads to negative values in fastest box with lifetime below 1 year.
        %cumulative emissions in GtC
        R(5,t,p)= R(5,t-1,p)+E(p,t); 
        
        %forcing of gaz p 
        %(element t has the average forcing between t and t-1)
        Conc = paramR(p,13) + (sum(R(1:4,t-1,p),1) + sum(R(1:4,t,p),1)) / 2 * paramR(p,14); %C0 + (Ga_t + Ga_t-1)/2*emis2conc
        R(7,t,p) = paramR(p,15) * log(Conc/paramR(p,13)) + paramR(p,16) * (Conc-paramR(p,13)) + paramR(p,17) * (sqrt(Conc)-sqrt(paramR(p,13))); %f1 * log(C/C0) + f2 * (C-C0) + f3 * (sqrtC-sqrtC0)
    end 
    
    %total forcing of all gases     
    F(t)=sum(R(7,t,:),3) + F_ext(t);
    %temperature in each of the 3 boxes
    for k=1:3
        S(k,t)=paramT(k+3) * F(t) * (1-exp(-1/paramT(k))) + S(k,t-1) * exp(-1/paramT(k)); %equation from SM9
    end 
    %temperature
    T(t)= sum(S(:,t),1); % SM9 takes the average over St and St-1 but i think that is wrong.
    if imag(T(t))~=0;disp(['Complex matrix detected at iteration' ]); keyboard;end
end
