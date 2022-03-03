%set parameter for LTR4state model
version='210425';
fig_path1=strcat('./',version);
mkdir(fig_path1);
fig_path=strcat(fig_path1,'/');

%% Set Parameter of LTR-8-State Model
%model parameter
alpha_NE=1;
alpha_NK=0;
alpha_NS=-1;
NM={'NE','NK','NS'};
omega=100;
gamma_AC=2.5e9;
gamma_untreated=1;
k_act=1e-11;
k_unact=0.1;
k_bindp=0.001;
k_unbindp=0.1;
beta=10;

DetailBalenceBreak={'DetailBalance','R-R*',...
    'R*-R','R*-R*P','R*P-R*'...
    'R*P-P','P-R*P','P-R','R-P',...
    'R*P-P&R*P-R*'};
i_DB=1;%i_DB=3 landscape for AC and NE
i_nDB=10;
cas=DetailBalenceBreak{i_DB};

gamma_vec=[gamma_untreated,gamma_AC];
gamma_vec_len=length(gamma_vec);
alpha_NM=[alpha_NE,alpha_NK,alpha_NS];
alpha_NM_len=length(alpha_NM);
alpha_vec=0:0.1:2;%%
alpha_vec_len=length(alpha_vec);
alpha_NE_NS_pos=11;
NM_color={'r','c','b'};
