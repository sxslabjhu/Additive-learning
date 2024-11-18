% %plot the results
% 
% Total_data=[X_prcs,Y_prcs];%Prcs_Log_All_data;
% Input_var=[1,2];
% % Dif_Input_var=setdiff(1:3,Input_var);
% Output_var=9; %These numbers correspond to "picked-out" variables (X_prcs,Y_prcs)
% load('condition1.mat');
% X_train=Log_All_data(1:20900,[2,4:6,8:11]);
% Y_train=Log_All_data(1:20900,7);

%% true vs. predicted plot

c_temp={c{end,:}};
for i=1:length(Y_test_prcs)
    X_temp=X_test_prcs(i,:);
    Y_test_prcs_predict(i,1)=c_temp{1}+c_temp{2}'*X_temp'+1/2*X_temp*c_temp{3}*X_temp';
end

scatter(Y_test_prcs,Y_test_prcs_predict,'.');
hold on;
plot([-2,8],[-2,8],'linewidth',2);
xlim([-2,8]);
ylim([-2,8]);
xlabel('True');
ylabel('Predicted');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
box on;
% l1=legend('Data point','y=x');
% set(l1,'box', 'off')

Err=abs(Y_test_prcs_predict-Y_test_prcs);
mean(Err)

%% Complete NN plot

node=[20,20];
Mdl = fitrnet(X_train_prcs,Y_train_prcs,"Standardize",true, ...
    "LayerSizes",node);

Y_train_prcs_predict=predict(Mdl,X_train_prcs);
Y_test_prcs_predict=predict(Mdl,X_test_prcs);
% Y_train_ori_predict=Y_train_prcs_predict*std(Y_train)+mean(Y_train);
% Y_test_ori_predict=Y_test_prcs_predict*std(Y_train)+mean(Y_train);

scatter(Y_test_prcs,Y_test_prcs_predict,'.');
hold on;
% plot([-1.5,1.5],[-1.5,1.5],'linewidth',2);
%plot(Y_test,Y_test,'linewidth',2)
hold off

xlabel("True")
ylabel("Predicted")
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1)
% xlim([-1.5,1.5]);
% ylim([-1.5,1.5]);
box on;

err_test=abs(Y_test_prcs-Y_test_prcs_predict);
mean(err_test)

%% 3D test of the result

Input_var=[1,2];
c_temp={c{end,:}};
x_range=-2.5:0.01:4.5;
y_range=x_range;
[X,Y]=meshgrid(x_range,y_range);
XX=x_range;
YY=y_range;
L=length(XX);
for i=1:L
    for j=1:L
        X_temp=zeros(l,1);
        X_temp(Input_var(1))=XX(i);
        X_temp(Input_var(2))=YY(j);
        Z(i,j)=c_temp{1}+c_temp{2}'*X_temp+1/2*X_temp'*c_temp{3}*X_temp;
    end
end
%Zero_projection_label=find(Total_data(:,Dif_Input_var)<0.1 & Total_data(:,Dif_Input_var)>-0.1);
%scatter3(Total_data(Zero_projection_label,Input_var(1)),Total_data(Zero_projection_label,Input_var(2)),Total_data(Zero_projection_label,Output_var));
scatter3(X_test_prcs(:,Input_var(1)),X_test_prcs(:,Input_var(2)),Y_test_prcs,'x');
hold on;
mesh(Y,X,Z);
xlabel(Mintoutlabel{inputnum(Input_var(1))});
ylabel(Mintoutlabel{inputnum(Input_var(2))});
zlabel('P53');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
l1=legend('True','Predicted');
set(l1,'box','off');
box on;

%% projection test of the result

for i=1:l
xx{i}=-2:0.01:4;
end

k_test=2;

% c_real={[0],[1;1;1],zeros(3,3)};
c_theory={c{end,:}};
% L_real=Mean_nD(xx,s,c_real);
L_theory=Mean_nD(xx,s,c_theory);
% var_real=Variance_nD(xx,s,c_real);
var_theory=Variance_nD(xx,s,c_theory);

xx1=xx{k_test};
% L_real1=L_real{1};
L_theory1=L_theory{k_test};
% var_real1=var_real{1};
var_theory1=var_theory{k_test};

figure(1);
scatter(xt{k_test},M_Exp{k_test},'x');
hold on;
% plot(xx1,L_real1,'b','linewidth',1.5);
% hold on;
plot(xx1,L_theory1,'r','linewidth',1.5);
xlabel(Mintoutlabel{inputnum(k_test)});
ylabel('$\mu (P_{53})$','interpreter','latex');
l1=legend('True','Predicted');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
set(l1,'box','off');
box on;
    
figure(2);
scatter(xt{k_test}(V_Exp{k_test}~=0),V_Exp{k_test}(V_Exp{k_test}~=0),'d');
hold on;
% plot(xx1(var_real1~=0),var_real1(var_real1~=0),'b','linewidth',1.5);
% hold on;
plot(xx1(var_theory1~=0),var_theory1(var_theory1~=0),'r','linewidth',1.5);
xlabel(Mintoutlabel{inputnum(k_test)});
ylabel('$\sigma ^2 (P_{53})$','interpreter','latex');
l2=legend('True','Predicted');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
set(l2,'box','off');
box on;

%% histogram comparison

c_temp={c{end,:}};
for i=1:length(Y_test_prcs)
    X_temp=X_test_prcs(i,:);
    Y_test_prcs_predict(i,1)=c_temp{1}+c_temp{2}'*X_temp'+1/2*X_temp*c_temp{3}*X_temp';
end

histogram(Y_test_prcs,'normalization','pdf');
hold on;
histogram(Y_test_prcs_predict,'normalization','pdf');
xlim([-2 8]);
xlabel('P53');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ylabel('P.D.F');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
l1=legend('True','Predicted');
set(l1,'box','off');
box on;

%% Error plot

plot(Err_all_step,'linewidth',2);
xlabel('Step');
ylabel('Loss function');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);

%% Linear and quadratic coefficients of the polynomial regression model of the P53
plot(c{end,2},'linewidth',2);
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
xlabel('variable number');
ylabel('$F_{\alpha}^{\prime}$','Interpreter','latex');

figure;
Fab=c{end,3};
Fab=[Fab,zeros(length(Fab),1);
    zeros(1,length(Fab)+1)];
h=pcolor(Fab);
% set(h, 'EdgeColor', 'none');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
xlabel('variable number');
ylabel('variable number');
colorbar;