%% prediction vs. true - testing set

% load('Spring_sig_0p02.mat');
for i=1:length(Y_train_prcs)
    Y_train_prcs_predict(i,1)=MyNN_Prediction(X_train_prcs(i,:),W); % need to be modified to be suitable for new normalization
end

for i=1:length(Y_test_prcs)
    Y_test_prcs_predict(i,1)=MyNN_Prediction(X_test_prcs(i,:),W); % need to be modified to be suitable for new normalization
end

scatter(Y_test_prcs,Y_test_prcs_predict,'.');
hold on;
plot([-1.5,1.5],[-1.5,1.5],'linewidth',2);
xlabel('True');
ylabel('Predicted');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
box on;

xlim([-1.5,1.5]);
ylim([-1.5,1.5]);

err_test=abs(Y_test_prcs-Y_test_prcs_predict);
mean(err_test)

%% Complete NN plot

node=[20,20];
Mdl = fitrnet(X_train_prcs,Y_train_prcs,"Standardize",true, ...
    "LayerSizes",node);

Y_train_prcs_predict=predict(Mdl,X_train_prcs);
Y_test_prcs_predict=predict(Mdl,X_test_prcs);

scatter(Y_test_prcs,Y_test_prcs_predict,'.');
hold on;

xlabel("True")
ylabel("Predicted")
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1)
xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
box on;

err_test=abs(Y_test_prcs-Y_test_prcs_predict);
mean(err_test)

%% projection test of the result

f=7; %variable#

figure(1);
scatter(xt{f},M_Exp{f},'x');
hold on;
plot(xt{f},Mean_PrdctBySmpl{f},'r','linewidth',2);
%xlabel(strcat('$x_',num2str(f),'$'),'interpreter','latex');
xlabel('$P_{22}$','interpreter','latex');
ylabel('$\mu (\delta y_2)$','interpreter','latex');
l1=legend('True','Predicted');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
set(l1,'box','off');
box on;
    
figure(2);
scatter(xt{f},V_Exp{f},'d');
hold on;
plot(xt{f},Var_PrdctBySmpl{f},'r','linewidth',2);
%xlabel(strcat('$x_',num2str(f),'$'),'interpreter','latex');
xlabel('$P_{22}$','interpreter','latex');
ylabel('$\sigma ^2 (\delta y_2)$','interpreter','latex');
l2=legend('True','Predicted');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
set(l2,'box','off');
box on;

%% histogram comparison
for i=1:length(Y_train_prcs)
    Y_train_prcs_predict(i,1)=MyNN_Prediction(X_train_prcs(i,:),W); % need to be modified to be suitable for new normalization
end

for i=1:length(Y_test_prcs)
    Y_test_prcs_predict(i,1)=MyNN_Prediction(X_test_prcs(i,:),W); % need to be modified to be suitable for new normalization
end

histogram(Y_test_prcs,'normalization','pdf');
hold on;
histogram(Y_test_prcs_predict,'normalization','pdf');
xlabel('$\delta y_2$','interpreter','latex');
ylabel('P.D.F');
% xlim([0,2e5])
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
l1=legend('True','Predicted');
set(l1,'box','off');
box on;

%% err_plot
plot(Err_all_step,'linewidth',2);
xlabel('Step');
ylabel('Loss function');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);

%% 3D test of the result - normalized data
x_range=-0.5:0.01:0.5;
y_range=x_range;
[X,Y]=meshgrid(x_range,y_range);
XX=x_range;
YY=y_range;
L=length(XX);
for i=1:L
    for j=1:L
        X_temp=zeros(d,1);
        X_temp(1)=XX(i);
        X_temp(7)=YY(j);
        Z(i,j)=MyNN_Prediction(X_temp',W);
    end
end

surface_pred=mesh(Y,X,Z);
surface_pred.FaceColor='flat';
hold on;
mesh(F2_scan,F1_scan,dX_scan');
l1=legend('Predicted','True');
xlabel('$P_{21}$','interpreter','latex');
%xlabel('${\mu}$','interpreter','latex', 'FontWeight','bold')
ylabel('$P_{22}$','interpreter','latex');
zlabel('$\delta y_2$','interpreter','latex');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
set(l1,'box','off');
box on;

%% Raw data visulization
clc;clear;
load('Fsig=0p02_SpringNetwork_NN_para_[20  20]_dsub=1_Lbin=0.02_nmin=2_Step_SA=0.05.mat');
c=histogram2(X_train_prcs(:,7),Y_train_prcs,'Normalization','pdf','DisplayStyle','tile');
set(c,'EdgeColor','none');
grid off;
xlabel('$P_{22}$','interpreter','latex');
ylabel('$\delta y_2$','interpreter','latex');
colorbar;
set(gca,'FontSize',20,'Fontname', 'Arial','linewidth',1);