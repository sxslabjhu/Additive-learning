%input: F_i
%output: delta x_i (vector)

clc;clear;

%% generate normaly distributed forces for ML

%
%First case with 12 DOF
X=[-2,0;-1,0.5;0,0.5;1,0.5;2,0;1,-0.5;0,-0.5;-1,-0.5]; %initial configuration
Node_fixed=[1,5];
output_displacement_node=[2,2];
%}

%{
%High dimensional case with 100 DOF
X=zeros(52,2); %initial configuration
aa=linspace(-5,5,10); 
bb=linspace(-2.5,2.5,5);
[xxx,yyy]=meshgrid(aa,bb);
X(2:51,:)=[xxx(:),yyy(:)];
Node_fixed=[1,52];
output_displacement_node=[5,2];
%}

%{
%2 dimensional case: testing the metastcble testing by clustering
X=[-1,0;1,0;0,1]; %initial configuration
Node_fixed=[1,2];
output_displacement_node=[3,2];
%}

%
%4 DOF
% X=[-1,0;0,1;1,0;0,-1];
% F=[0,0;0,0.1;0,0;0,0];
% Node_fixed=[1,3];
% output_displacement_node=[2,2];
%}

N=length(X);  
k_spring=ones(N,N)-eye(N); %spring stiffness
l0=zeros(N,N); 
for i=1:N
    for j=1:N
        l0(i,j)=norm(X(i,:)-X(j,:));
    end
end   %rest length corresponds to the initial configuration
Nf=N-2; %number of free nodes (except 1 and 5)
MU=zeros(Nf*2,1); 
SIGMA=0.1*eye(Nf*2);
N_sample=5000;%Nf*2*1500; %statistical parameters of random force

dX_training=zeros(N_sample,1);
dF_training=mvnrnd(MU,SIGMA,N_sample);
for i=1:N_sample
    disp(strcat('PART I - generating original data -',num2str(i/N_sample*100),'%'))
    F=zeros(N,2);
    F(setdiff(1:N,Node_fixed),:)=reshape(dF_training(i,:),[Nf,2]);
    dX=calc_def(X,Node_fixed,l0,k_spring,F);
    dX_training(i)=dX(output_displacement_node(1),output_displacement_node(2)); %here the output displacement is Y displacement at node 2
end

%% generate parameter scan plot (e.g. F1_scan: Fx2, F2_scan: Fy2, dX_scan: dY2)
% F1_scan=-0.5:0.025:0.5;
% F2_scan=-0.5:0.025:0.5;
% dX_scan=zeros(length(F1_scan),length(F2_scan));
% 
% for i=1:length(F1_scan)
%     for j=1:length(F2_scan)
%         F_total_scan=zeros(N,2);
%         F_total_scan(output_displacement_node(1),1)=F1_scan(i);
%         F_total_scan(output_displacement_node(1),2)=F2_scan(j);
%         dX_total_scan=calc_def(X,Node_fixed,l0,k_spring,F_total_scan);
%         dX_scan(i,j)=dX_total_scan(output_displacement_node(1),output_displacement_node(2));
%     end
%     disp(strcat('PART II - generating scan plot -',num2str(i/length(F1_scan)*100),'%'))
% end