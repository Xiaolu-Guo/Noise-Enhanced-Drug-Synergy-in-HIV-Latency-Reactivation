
%% Set Parameter of LTR-8-State Model
%SetParameterLTR4State_V7
koff_ver='/koff_Hill_';%'_koff_Sigmoid'

%% Set Parameter of Tat positive feedback SSA Model
total_cell=10000;%20000*100000
time_steps=100000;
if ~exist('on_threshold')
    on_threshold=75;
end

k_mbasal=0.01;
k_Tat=10;
d_m=1;
d_Tat=0.125;
k_trs1=5;
k_trs2=1;
observ_time=100;


version=strcat(version,koff_ver);
switch koff_ver
    case '/koff_Hill_'
        koff_ratio=@(Tat,on_threshold) (on_threshold^3+1e-2*Tat.^3)./(on_threshold^3+Tat.^3);
        
    case '/koff_Sigmoid_'
        koff_ratio=@(Tat,on_threshold) (1+1e-3*exp(Tat-on_threshold))./(1+0.1*exp(Tat-on_threshold));
        
end
