%% GMM fitting performance
clc;clear;
load('Senescence_cellcondition1_prcs.mat');
nG=2;
[~,d]=size(X_train_prcs);
X_prcs=X_train_prcs;
npair=[1,2];
for k=1:d
%   [X_prcs,Y_prcs]=DataPreprocessing(DataType);
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