% Polynmomial regression of partial faceted data

function Polynomial_main(job_identifier)
%% Data preprocessing & model parameters

DataType_all=[1 2 3 4];
Lbin_all=[0.01 0.02 0.035 0.05 0.1 0.2];
nmin_all=[3 4 5];
dc_all=[0.01 0.025 0.05];
[sub1,sub2,sub3,sub4]=ind2sub([length(DataType_all),length(Lbin_all),length(nmin_all),length(dc_all)],job_identifier);
DataType=DataType_all(sub1);
Lbin=Lbin_all(sub2);
nmin=nmin_all(sub3);
dc=dc_all(sub4);

%DataType=3;% 1:Spring netwrok 2:AcH3 3:Scenescence

if DataType==1
    dataset_label='SpringNetwork_Sig_0p02';
    load("Spring_sig_0p02.mat");
    X_prcs=dF_training;
    Y_prcs=dX_training;
elseif DataType==2
    dataset_label='SpringNetwork_Sig_0p1';
    load("Spring_sig_0p1.mat");
    X_prcs=dF_training;
    Y_prcs=dX_training;
elseif DataType==3
    dataset_label='Scenescence_cond1';
    load('Senescence_cellcondition1_prcs.mat');
    X_prcs=X_train_prcs;
    Y_prcs=Y_train_prcs;
else
    dataset_label='Scenescence_cond3';
    load('Senescence_cellcondition3_prcs.mat');
    X_prcs=X_train_prcs;
    Y_prcs=Y_train_prcs;
end

% [X_prcs,Y_prcs]=DataPreprocessing(DataType);

[N,l]=size(X_prcs);
[~,q]=size(Y_prcs);
s=cov(X_prcs);

mu=zeros(1,l);
x_diff_measure=cell(1,l);
y_diff_measure=cell(1,l);
x_bin=cell(1,l);
N_sub_data=floor(N/l);

% nmin=5; %minimum number of points required in each bin for statistical analysis 
% Lbin=5e-2; %bin size

%% conditional mean and variance of experiemtal data
for i=1:l
    x_diff_measure{i}=X_prcs((1+(i-1)*N_sub_data):i*N_sub_data,:);
    y_diff_measure{i}=Y_prcs((1+(i-1)*N_sub_data):i*N_sub_data,:);
    x_bin{i}=min(x_diff_measure{i}(:,i)):Lbin:max(x_diff_measure{i}(:,i));%divide the data into several consecutive bins
    k=1;
    for j=1:length(x_bin{i})-1 %j=1:length; later use i without empty set case
        Bt{i}{j}=find(x_diff_measure{i}(:,i)>x_bin{i}(j) & x_diff_measure{i}(:,i)<=x_bin{i}(j+1));
        if length(Bt{i}{j})>nmin
            B1{i}{k}=Bt{i}{j};
            xt{i}(k)=(x_bin{i}(j)+x_bin{i}(j+1))/2;
            M_Exp{i}(k)=mean(y_diff_measure{i}(B1{i}{k}));
            V_Exp{i}(k)=var(y_diff_measure{i}(B1{i}{k}));
            k=k+1;
        end
    end
end

%% optimization using simulated annealing method
%dc=0.05;

temp=2*rand(l,l);
c3t=tril(temp,-1)+triu(temp',0);

c={2*rand,2*rand(l,1),c3t};
T(1)=1e5*1;
alpha=0.95;

imax=50000;
err=inf;
p=2;
ct=cell(1,3);
for i=1:imax
    i/imax
    for f=1:2
        ct{f}=c{end,f}+dc*(rand(size(c{end,f}))-0.5);
    end
    if DataType==1
        ct{1}=0;
    end
    temp=rand(l,l)-0.5;
    ct{3}=c{end,3}+dc*(tril(temp,-1)+triu(temp',0));
    Mf1=Mean_nD(xt,s,ct);
    Vf2=Variance_nD(xt,s,ct);
    delta=0;
    for q=1:l
        deltaM=Mf1{q}-M_Exp{q};
        deltaV=Vf2{q}-V_Exp{q};
        delta=delta+norm(deltaM)^2+norm(deltaV)^2;
    end
    if delta<err(p-1)
        Err_all_step(i)=delta;
        err(p)=delta;
        c{p,1}=ct{1};
        c{p,2}=ct{2};
        c{p,3}=ct{3};
        p=p+1;
    else
        Err_all_step(i)=Err_all_step(i-1);
        if rand()<exp(-(delta-err(p-1))/T(i))
            err(p)=delta;
            Err_all_step(i)=delta;
            c{p,1}=ct{1};
            c{p,2}=ct{2};
            c{p,3}=ct{3};
            p=p+1;
        end
    end
    T(i+1)=T(i)*alpha;
    save(strcat('Poly_',dataset_label,...
        '_Lbin=',num2str(Lbin),...
        '_nmin=',num2str(nmin),...
        '_Step_SA=',num2str(dc),...
        '.mat'));
end