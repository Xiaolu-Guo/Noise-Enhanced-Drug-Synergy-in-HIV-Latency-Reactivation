ndb = 0;
d_Tat_ornot = 1;
d_Tat0 = 0;

cell_cycle = 4;


%% addpath so all the functions in the filefolder can be called
%clear all

%lambda on and off run figs4_draw_lambda_onoff_V7

% change the path
% cd /Users/admin/Documents/MATLAB/HIV/EssayCodesV7
% addpath('/Users/admin/Documents/MATLAB/HIV/EssayCodesV7')
% cd /Users/Admin/.nutstore/.nutstore/MATLAB/HIV/EssayCodesV7
% addpath('/Users/Admin/.nutstore/.nutstore/MATLAB/HIV/EssayCodesV7')
% clear all
if ~exist('reactivation_or_not')
    reactivation_or_not=1;%whether calculate Reactivation or not; Recativation calculation will cost more time
    %reactivation_or_not=1: calculate reactivation and Pon.
    %reactivation_or_not=0: calculate only Pon.
end

if ~exist('rundata')
    rundata=1;
    % rundata=1: run and save data.
    % rundata=0: load saved data.
end

if ~exist('ndb')
    ndb=1;
    % ndb=1: calculate non-detailed balance results in fig6.
    % ndb=0: calculate detailed balance results in fig4.
end

%% Set Paramters of the model and drawing fig
run('../EssayCodesV8/SetParameterLTR4State_V7.m');
run('../EssayCodesV8/SetParameterTatPosFeed_V7.m');

%% initialise variables
gamma_vec0=[1,(1:25)*10^8];
Pon=zeros(2,alpha_vec_len);%first row represents without AC; second row represents with AC
MDT_off=zeros(2,alpha_vec_len,alpha_NM_len);
MDT_on=zeros(2,alpha_vec_len,alpha_NM_len);
Pon_bar=zeros(5,1);
Reac_bar=zeros(5,1);
total_cell = 10000;

if d_Tat_ornot
    d_Tat = d_Tat0;
end

d_Tat_str = replace(num2str(d_Tat),'.','p');

if ndb
    cas='R*P-P&R*P-R*';
    beta=10;
    fig4n6_data_filename=strcat('./data_mat_form/fig6_data_reac_ndb_cellcycle_',num2str(cell_cycle),'hours_dTat_',d_Tat_str);
else
    cas='DetailBalance';
    beta=0;
    fig4n6_data_filename=strcat('./data_mat_form/fig4_data_reac_db_cellcycle_',num2str(cell_cycle),'hours_dTat_',d_Tat_str);
end

%% Figure3n4 Pon-AC-N

i_vec=[1,2,1,2,2];
j_vec=[1,1,alpha_NE_NS_pos,alpha_NE_NS_pos,alpha_NE_NS_pos];
i_NM_vec=[2,2,1,1,3];


for i_bar=4 %1:5% i_bar=1:5 correspondind to untreated, AC only, NE only, AC+NE, AC+NS
    
    gamma=gamma_vec(i_vec(i_bar));
    
    alpha=alpha_NM(i_NM_vec(i_bar))*alpha_vec(j_vec(i_bar));%set parameter for NE or NK or ND
    total_cell = 20;
    
    LTR4State_TatPosFeed_SSA_ReacRatio_MFPT_CellCycle_v2_traj...
        (k_act,k_unact,k_bindp,k_unbindp,gamma,omega,alpha,beta,cas,...
        total_cell,time_steps,observ_time,on_threshold,k_mbasal,k_Tat,d_Tat,k_trs1,k_trs2,d_m,koff_ratio,cell_cycle);
LTR4State_TatPosFeed_SSA_ReacRatio_MFPT_v2_traj...
        (k_act,k_unact,k_bindp,k_unbindp,gamma,omega,alpha,beta,cas,...
        total_cell,time_steps,observ_time,on_threshold,k_mbasal,k_Tat,d_Tat,k_trs1,k_trs2,d_m,koff_ratio,cell_cycle);

end


%


