function [] = LTR4State_TatPosFeed_SSA_ReacRatio_MFPT_CellCycle_v2_traj...
    (k_act,k_unact,k_bindp,k_unbindp,gamma,omega,alpha,beta,cas,...
    total_cell,time_steps,observ_time,on_threshold,k_mbasal,k_Tat,d_Tat,k_trs1,k_trs2,d_m,koff_ratio,cell_cycle)
%LTR8StatePon Calculate the mean first passage time of Polymerase binding to or unbinding from LTR in HIV
%LTR-8-State Markov model and the Probability of Polymerase binding to LTR
%   A: Activator
%   N: Noise modulator
%   P: Polymerase

run('../EssayCOdesV8/SetJumpRateMetrix_V7.m');

Q0_nonzero=Q0_nonzero';
Q0=Q0';

%% Set value for Reaction Vector and Q rate????
Q_len=sum(sum(Q0_nonzero));
num_variables=6;
num_molecular_reac=4;
num_reactions=Q_len+num_molecular_reac;
[Q_raw,Q_col,Q_rate]=find(Q0_nonzero);
reaction_vector=zeros(num_reactions,num_variables);
for i_Q=1:Q_len
    reaction_vector(i_Q,Q_raw(i_Q))=1;
    reaction_vector(i_Q,Q_col(i_Q))=-1;
    Q_rate(i_Q)=Q0(Q_raw(i_Q),Q_col(i_Q));
end
reaction_vector((Q_len+1):num_reactions,end-1:end)=...
    [1,0
    0,1
    -1,0
    0,-1];

%% Set initial values of variables
%total_cell = 100;
Tat0=0;
zero_vector=zeros(total_cell,1); % initialize the vector
one_vector=ones(total_cell,1);
one_array=ones(1,num_reactions);

LTR=1*ones(total_cell,time_steps);
LTR_RStar=zeros(total_cell,time_steps);
LTR_RStarP=zeros(total_cell,time_steps);
LTR_P=0*ones(total_cell,time_steps);

mRNA=0*ones(total_cell,time_steps);
Tat=Tat0*ones(total_cell,time_steps);
T=zeros(total_cell,time_steps+1); % the reaction time matrix
koff_decrease=ones(total_cell,1);

%% variables monitor to stop
index_t=1:total_cell;
index_t =index_t';
first_passage_time=-ones(total_cell,1);
all_reac=0;%record whether all cells have evoluted to the observ time

% cell_cycle = 8;
time_num =10;
t_division = ones(total_cell,1) * (cell_cycle : cell_cycle : (time_num * observ_time+cell_cycle));
% d_Tat = 0;
t_division_index = zeros(total_cell,1);

for t=1:time_steps-1
    k_m=k_mbasal+k_trs1*(Tat(index_t,t)./k_trs2)./(1+Tat(index_t,t)./k_trs2);
    koff_decrease(index_t)=koff_ratio(Tat(index_t,t),on_threshold);%if Tat>thresh then reduce koff
    Rate_matrix=[ Q_rate(1)*LTR(index_t,t),     Q_rate(2)*LTR(index_t,t), ...
        Q_rate(3)*LTR_RStar(index_t,t),   Q_rate(4)*LTR_RStar(index_t,t), ...
        Q_rate(5)*LTR_P(index_t,t).*koff_decrease(index_t),  Q_rate(6)*LTR_P(index_t,t), ...
        Q_rate(7)*LTR_RStarP(index_t,t).*koff_decrease(index_t),  Q_rate(8)*LTR_RStarP(index_t,t),...
        k_m.*(LTR_P(index_t,t)+LTR_RStarP(index_t,t)),...
        k_Tat*mRNA(index_t,t),...
        d_m.*mRNA(index_t,t),...
        d_Tat*Tat(index_t,t)]; % calculate the reaction rate of each reaction at present
    Rate_matrix_sum=cumsum(Rate_matrix,2);
    %rng(14);
    tau=exprnd(one_vector(index_t)./Rate_matrix_sum(:,num_reactions)); % generate the random value of next reation time
    
    T(index_t,t+1)=T(index_t,t)+tau; % update the time of each cell going through
    
    
    fire_reaction=unifrnd(zero_vector(index_t),Rate_matrix_sum(:,num_reactions))*one_array; % generate random value to decide which reaction will be fired in tau
    K=sum(cumsum(fire_reaction<Rate_matrix_sum,2)==0,2)+1; %decide which reaction will be fired in tau
    LTR(index_t,t+1)=LTR(index_t,t)+reaction_vector(K,1);
    LTR_RStar(index_t,t+1)=LTR_RStar(index_t,t)+reaction_vector(K,2);
    LTR_P(index_t,t+1)=LTR_P(index_t,t)+reaction_vector(K,3);
    LTR_RStarP(index_t,t+1)=LTR_RStarP(index_t,t)+reaction_vector(K,4);
    mRNA(index_t,t+1)=mRNA(index_t,t)+reaction_vector(K,5);
    Tat(index_t,t+1)=Tat(index_t,t)+reaction_vector(K,6);
    
    index_cell_div = sum( T(index_t,t+1)>= t_division(index_t,:) ,2) ~= sum( T(index_t,t)>= t_division(index_t,:) ,2) ;
    t_division_index(index_t) = sum( T(index_t,t)>= t_division(index_t,:),2) +1;        
    while max(t_division_index)>size(t_division,2)
        t_division = ones(total_cell,1) * (cell_cycle : cell_cycle : (time_num * observ_time+cell_cycle));
        time_num = time_num+10;
        index_cell_div = sum( T(index_t,t+1)>= t_division(index_t,:) ,2) ~= sum( T(index_t,t)>= t_division(index_t,:) ,2) ;
        t_division_index(index_t) = sum( T(index_t,t)> t_division(index_t,:),2) +1;
    end
    T(index_t(index_cell_div),t+1)=t_division(index_t(index_cell_div)+(t_division_index(index_t(index_cell_div))-1)*total_cell);
    LTR(index_t(index_cell_div),t+1)=LTR(index_t(index_cell_div),t);
    LTR_RStar(index_t(index_cell_div),t+1)=LTR_RStar(index_t(index_cell_div),t);
    LTR_P(index_t(index_cell_div),t+1)=LTR_P(index_t(index_cell_div),t);
    LTR_RStarP(index_t(index_cell_div),t+1)=LTR_RStarP(index_t(index_cell_div),t);
    mRNA(index_t(index_cell_div),t+1)=mRNA(index_t(index_cell_div),t);
    Tat(index_t(index_cell_div),t+1)=ceil(Tat(index_t(index_cell_div),t)/2);
    
    index=Tat(index_t,t+1)<on_threshold;
    first_passage_time(index_t(~index))=T(index_t(~index),t+1);
    % index_t=index_t(index);
%     if isempty(index_t)
%         all_reac=1;
%         break
%     end
end

%% Calculate MFPT and ReacRatio
% if all_reac
%     MFPT_GFPmoreThresh=mean(first_passage_time);
%     ReacRatio=sum(first_passage_time<=observ_time)/total_cell;
% else
%     index_narr=first_passage_time<0;
%     if T(index_narr,time_steps)<=0
%         error('time calculating wrong!')
%     end
%     first_passage_time(index_narr)=T(index_narr,time_steps)*2;
%     MFPT_GFPmoreThresh=mean(first_passage_time);
%     ReacRatio=sum(first_passage_time<=observ_time)/total_cell;
% end

% debug
for i_cell = 1:total_cell
    plot(T(i_cell,Tat(i_cell,:)>0),Tat(i_cell,Tat(i_cell,:)>0),'LineWidth',2); hold on
end
xlim([0,100])
xlabel('Time (hour)');
ylabel('Tat');

set(gca,'FontSize',14,'FontWeight','b')

end






