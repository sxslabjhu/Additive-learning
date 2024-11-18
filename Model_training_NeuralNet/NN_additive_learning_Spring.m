function NN_additive_learning_Spring2(job_identifier)
% ---------------- Parameter for NN ------------------
%clc;clear;
DataType=1;
node=[20,20];

d_sub_all=[1 2 4];
Lbin_all=[0.02 0.05 0.1 0.2 0.3 0.5];
nmin_all=[1 2 3 4];

% d_sub_all=[1 2 4];
% Lbin_all=[0.05 0.5];
% nmin_all=[1 2 4];

[sub1,sub2,sub3]=ind2sub([length(d_sub_all),length(Lbin_all),length(nmin_all)],job_identifier);
d_sub=d_sub_all(sub1);
Lbin=Lbin_all(sub2);
n_min=nmin_all(sub3);

%parameters for SA
dl_MC=0.05;
imax=5e4;
n_sample=60000;
Err=inf;
i_ErrUpdate=1;
T(1)=1e5*1;
alpha=0.95;

%d_sub=1; %dimension of each subgroup
%Lbin=5e-2; %bin size for calculating mean values of experimental data
%n_min=4; %minimum # of data pts for a bin

%fprintf('Job Identifier: %s\n', job_identifier);
%fprintf('d_sub: %s\n', d_sub);

strcat('Job Identifier = ',num2str(job_identifier))
strcat('d_sub = ',num2str(d_sub))
strcat('L_bin = ',num2str(Lbin))
strcat('nmin = ',num2str(n_min))

load('Spring_sig_0p1.mat');
X_prcs=dF_training;
Y_prcs=dX_training;

if DataType==1
    dataset_label='SpringNetwork';
elseif DataType==2
    dataset_label='Nucleus';
else
    dataset_label='Scenescence';
end

[N,d]=size(X_prcs);
[~,q]=size(Y_prcs);
S=cov(X_prcs);

%Divide the complete dataset into several subgroups
n_partial_set=ceil(d/d_sub); %# of partial data sets
i_partial=zeros(n_partial_set,d_sub); %which variables are selected for each subgroup
for i=1:n_partial_set
    i_partial(i,:)=(1+(i-1)*d_sub):i*d_sub;
end
if mod(d,d_sub)~=0
    i_partial(n_partial_set,:)=(d-d_sub+1):d;
end

%getting the partial data and statistics
N_sub=floor(N/n_partial_set); %number of data points in each set of measurement (floor is used in case it's out of range)
X_partial=zeros(N_sub*n_partial_set,d_sub); %All partial independent variables list together
Y_partial=zeros(N_sub*n_partial_set,1); %All partial dependent variables
A=zeros(d-d_sub,d-d_sub,n_partial_set);
B=zeros(d-d_sub,d_sub,n_partial_set);
D=zeros(d_sub,d_sub,n_partial_set);
mu1_cf=zeros(d-d_sub,d_sub,n_partial_set);
s1=zeros(d-d_sub,d-d_sub,n_partial_set); %Coefficients for calculating mean and variance

%Calculating mean and variance for exp. data AND conditional distribution (s1) for
%theoretical prediction
for i=1:n_partial_set
    X_partial((i-1)*N_sub+1:i*N_sub,:)=X_prcs((i-1)*N_sub+1:i*N_sub,i_partial(i,:));
    Y_partial((i-1)*N_sub+1:i*N_sub,1)=Y_prcs((i-1)*N_sub+1:i*N_sub,:);

    X_temp=X_partial((i-1)*N_sub+1:i*N_sub,:);
    Y_temp=Y_partial((i-1)*N_sub+1:i*N_sub,1);

    x_minmax=[min(X_temp);max(X_temp)]'; %[min,max;min,max;...]

    position{i}=ceil((X_temp-min(X_temp))/Lbin);
    N_all_dim{i}=ceil((max(X_temp)-min(X_temp))/Lbin);

    [position_net{i},ia{i},ic{i}]=unique(position{i},'row');% EXPERIMETNAL: Divide exp. data into bins

    k=1;
    for j=1:length(ia{i}) %j=1:length; later use i without empty set case
        if sum(ic{i}==j)>n_min
            xt{i}(k,:)=min(X_temp)+(position_net{i}(j,:)-0.5)*Lbin;
            M_Exp{i}(k)=mean(Y_temp(ic{i}==j));
            V_Exp{i}(k)=var(Y_temp(ic{i}==j));
            k=k+1;
        end
    end

    S_temp=S;
    col=S_temp(:,i_partial(i,:));
    S_temp(i_partial(i,:),:)=[];
    S_temp(:,i_partial(i,:))=[];
    A(:,:,i)=S_temp;
    D(:,:,i)=col(i_partial(i,:),:);
    col(i_partial(i,:),:)=[];
    B(:,:,i)=col;
    mu1_cf(:,:,i)=B(:,:,i)/D(:,:,i);
    s1(:,:,i)=A(:,:,i)-B(:,:,i)/D(:,:,i)*B(:,:,i)'; %Coefficients for estimating mean and variance
end

%network info
n_HL=length(node); % # of hidden layer
node_all=[d,node,q];
W=cell(1,n_HL+1); 
for i=1:n_HL+1
    W{i}=zeros(node_all(i)+1,node_all(i+1));%rand(node_all(i)+1,node_all(i+1));
end

%main program
for i=1:imax
    W_temp=W;
    for j=1:n_HL+1
        W_temp{j}=W_temp{j}+dl_MC*(rand(size(W_temp{j}))-0.5);
    end

    %mean value:
    delta=0;
    for j=1:n_partial_set
        for k=1:length(xt{j})
            X_sample=mvnrnd(mu1_cf(:,:,j)*xt{j}(k,:)',s1(:,:,j),n_sample);%sample for MC estimation of the mean value
            X_temp=X_sample;
            X_temp=[X_temp,ones(n_sample,1).*xt{j}(k,:)];
            if j~=n_partial_set
                X_temp(:,i_partial(j,1):end)=[ones(n_sample,1).*xt{j}(k,:),X_temp(:,i_partial(j,1):end-d_sub)];%%%%%%%%%
            end
            Y_PrdctBySmpl=MyNN_Prediction(X_temp,W_temp);
            Mean_PrdctBySmpl{j}(k)=mean(Y_PrdctBySmpl);
            Var_PrdctBySmpl{j}(k)=mean((Y_PrdctBySmpl-Mean_PrdctBySmpl{j}(k)).^2);
        end
        deltaM=Mean_PrdctBySmpl{j}-M_Exp{j};
        deltaV=Var_PrdctBySmpl{j}-V_Exp{j};
        delta=delta+norm(deltaM)^2+norm(deltaV)^2; %The weight should be adjusted propoerly !!!!!!!!!!!!!!!!!!!!!!!!!!
    end
    delta=delta/length(xt);

    strcat('step',num2str(i),' of ',num2str(imax),' - Process:',num2str(i/imax*100),'%')

    if delta<Err
        W=W_temp;
        Err=delta;
        Err_total(i_ErrUpdate)=delta;
        i_ErrUpdate=i_ErrUpdate+1;
    else
        if exp((Err-delta)/T(i))>rand
            W=W_temp;
            Err=delta;
            Err_total(i_ErrUpdate)=delta;
            i_ErrUpdate=i_ErrUpdate+1;
        end
    end
    T(i+1)=T(i)*alpha;
    Err_all_step(i)=Err;
    save(strcat(dataset_label,'_NN_para_[',num2str(node),']_dsub=',num2str(d_sub),'_Lbin=',num2str(Lbin),'_nmin=',num2str(n_min),...
        '_Step_SA=',num2str(dl_MC),'.mat'));
end