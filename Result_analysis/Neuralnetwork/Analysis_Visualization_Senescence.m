%% prediction vs. true - testing set

for i=1:length(Y_train_prcs)
    Y_train_prcs_predict(i,1)=MyNN_Prediction(X_train_prcs(i,:),W);
end

for i=1:length(Y_test_prcs)
    Y_test_prcs_predict(i,1)=MyNN_Prediction(X_test_prcs(i,:),W);
end

Y_train_ori_predict=Y_train_prcs_predict*std(Y_train)+mean(Y_train);
Y_test_ori_predict=Y_test_prcs_predict*std(Y_train)+mean(Y_train);

scatter(Y_test_prcs,Y_test_prcs_predict,'.');
hold on;
plot([-2.5 8.5],[-2.5 8.5],'linewidth',2);
xlabel('True');
ylabel('Predicted');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
box on;
% axis equal;
xlim([-2.5 8.5]);
ylim([-2.5 8.5]);

err_test=abs(Y_test_prcs-Y_test_prcs_predict);
mean(err_test)

%% Same model - prediction on different cell conditions
load('Senescence_cellcondition2_prcs.mat');
Y_train_prcs_predict=[];
for i=1:length(Y_train_prcs)
    Y_train_prcs_predict(i,1)=MyNN_Prediction(X_train_prcs(i,:),W); % need to be modified to be suitable for new normalization
end
scatter(Y_train_prcs,Y_train_prcs_predict,'.');
hold on;

load('Senescence_cellcondition3_prcs.mat');
Y_train_prcs_predict=[];
for i=1:length(Y_train_prcs)
    Y_train_prcs_predict(i,1)=MyNN_Prediction(X_train_prcs(i,:),W); % need to be modified to be suitable for new normalization
end
scatter(Y_train_prcs,Y_train_prcs_predict,'.');
hold on;

load('Senescence_cellcondition4_prcs.mat');
Y_train_prcs_predict=[];
for i=1:length(Y_train_prcs)
    Y_train_prcs_predict(i,1)=MyNN_Prediction(X_train_prcs(i,:),W); % need to be modified to be suitable for new normalization
end
scatter(Y_train_prcs,Y_train_prcs_predict,'.');
hold on;

%% Comparison between cross predict model and full NN
load('Senescence_cellcondition4_prcs.mat');
Y_train_prcs_predict=[];
for i=1:length(Y_train_prcs)
    Y_train_prcs_predict(i,1)=MyNN_Prediction(X_train_prcs(i,:),W); % need to be modified to be suitable for new normalization
end
scatter(Y_train_prcs,Y_train_prcs_predict,'.');
hold on;

node=[20,20];
Mdl = fitrnet(X_train_prcs,Y_train_prcs,"Standardize",true, ...
    "LayerSizes",node);
Y_train_prcs_predict=[];
Y_train_prcs_predict=predict(Mdl,X_train_prcs);
scatter(Y_train_prcs,Y_train_prcs_predict,'.');
hold on;

plot([-2.5 8.5],[-2.5 8.5],'linewidth',2);
%plot(Y_test,Y_test,'linewidth',2);
xlabel('True');
ylabel('Predicted');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
box on;
% axis equal;
xlim([-2.5 8.5]);
ylim([-2.5 8.5]);

l1=legend('FL, trained by condition 1','Full NN, trained by condition 4');
set(l1,'box','off');

%% Complete NN plot

node=[20,20];
Mdl = fitrnet(X_train_prcs,Y_train_prcs,"Standardize",true, ...
    "LayerSizes",node);

Y_train_prcs_predict=predict(Mdl,X_train_prcs);
Y_test_prcs_predict=predict(Mdl,X_test_prcs);

scatter(Y_test_prcs,Y_test_prcs_predict,'.');
hold on;
plot([-2.5 8.5],[-2.5 8.5],'linewidth',2);
hold off

xlabel("True")
ylabel("Predicted")
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1)
xlim([-2.5 8.5]);
ylim([-2.5 8.5]);
box on;

err=abs(Y_test_prcs-Y_test_prcs_predict);
mean(err)

%% projection test of the result

f=2; %variable#

figure(1);
scatter(xt{f},M_Exp{f},'x');
hold on;
plot(xt{f},Mean_PrdctBySmpl{f},'r','linewidth',2);
xlabel(Mintoutlabel{inputnum(f)});
ylabel(strcat('$\mu $(',Mintoutlabel{outputnum},')'),'interpreter','latex');
l1=legend('Exp.','Theory');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
set(l1,'box','off');
box on;
    
figure(2);
scatter(xt{f},V_Exp{f},'d');
hold on;
plot(xt{f},Var_PrdctBySmpl{f},'r','linewidth',2);
xlabel(Mintoutlabel{inputnum(f)});
ylabel(strcat('$\sigma ^2$(',Mintoutlabel{outputnum},')'),'interpreter','latex');
l2=legend('Exp.','Theory');
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
xlabel('P53');
ylabel('P.D.F');
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
f_plt=[1,2]; %var #
min=-2;
max=2;
Z_exp=Y_test_prcs;

X_fix=0;
[X1,X2]=meshgrid(min:0.1:max,min:0.1:max);
X_temp_plt=X_fix*ones(length(X1(:)),d);
X_temp_plt(:,f_plt)=[X1(:),X2(:)];
scatter3(X_test_prcs(:,f_plt(1)),X_test_prcs(:,f_plt(2)),Z_exp,'x');
hold on;
[H,Width]=size(X1);
Z_predict=MyNN_Prediction(X_temp_plt,W);
Z_predict=reshape(Z_predict,H,Width);
mesh(X1,X2,Z_predict);

xlim([min,max]);
ylim([min,max]);
l1=legend('True','Predicted');
xlabel(num2str(Mintoutlabel{inputnum(f_plt(1))}));
ylabel(num2str(Mintoutlabel{inputnum(f_plt(2))}));
zlabel('P53');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
set(l1,'box','off');
box on;

%% Raw data visulization
clc;clear;
load('Senescence_cellcondition1_prcs.mat');
c=histogram2(X_train(:,2),Y_train,'Normalization','pdf','DisplayStyle','tile');
set(c,'EdgeColor','none');
grid off;
xlabel(Mintoutlabel{inputnum(2)});
ylabel('P53');
colorbar;
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);

%% GMM fitting performance
clc;clear;
load('Senescence_cellcondition1_prcs.mat');
nG=2;
[~,d]=size(X_train_prcs);
X_prcs=X_train_prcs;
npair=[1,2];
for k=1:d
    gmm=fitgmdist(X_prcs(:,k),nG);
    mu_temp(k,:)=gmm.mu;
    var_temp0=reshape(gmm.Sigma,[nG,1,1]);
    var_temp(k,:)=var_temp0';
    Coef_temp(k,:)=gmm.ComponentProportion;
end

i1=npair(1);
i2=npair(2);

ct=1;
for i=1:nG
    for j=1:nG
        MU(ct,:)=[mu_temp(i1,i),mu_temp(i2,j)];
        SIGMA(1,:,ct)=[var_temp(i1,i),var_temp(i2,j)];
        P(ct,1)=Coef_temp(i1,i)*Coef_temp(i2,j);
        ct=ct+1;
    end
end

gm_pred = gmdistribution(MU,SIGMA,P);
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm_pred,[x0 y0]),x,y);
h = histogram2(X_train_prcs(:,1),X_train_prcs(:,2),'DisplayStyle','tile','ShowEmptyBins','on');
set(h,'EdgeColor','None');
hold on;
fcontour(gmPDF,'LineWidth',1.5);
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
xlabel(Mintoutlabel{inputnum(i1)});
ylabel(Mintoutlabel{inputnum(i2)});
xlim([-2,2]);
ylim([-2,2]);

figure;
h2 = histogram2(X_train(:,1),X_train(:,2),'DisplayStyle','tile','ShowEmptyBins','on');
set(h2,'EdgeColor','None');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
xlabel(Mintoutlabel{inputnum(i1)});
ylabel(Mintoutlabel{inputnum(i2)});
xlim([0,2.5e5]);
ylim([0,10e5]);

%% Different distribution in different conditions
index=8;
figure;
load('Senescence_cellcondition1_prcs.mat')
h1=histogram(X_train(:,index),'Normalization','pdf');hold on;
set(h1,'edgecolor','none')
load('Senescence_cellcondition2_prcs.mat')
h2=histogram(X_train(:,index),'Normalization','pdf');hold on;
set(h2,'edgecolor','none')
load('Senescence_cellcondition3_prcs.mat')
h3=histogram(X_train(:,index),'Normalization','pdf');hold on;
set(h3,'edgecolor','none')
load('Senescence_cellcondition4_prcs.mat')
h4=histogram(X_train(:,index),'Normalization','pdf');hold on;
set(h4,'edgecolor','none')
xlabel(Mintoutlabel{inputnum(index)});
ylabel('P.D.F.');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
legend('control','quiescent','50uM Bleo','200 nM Doxo');