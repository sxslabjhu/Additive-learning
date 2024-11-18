%plot the results

%% 3D test of the result
c_temp={c{end,:}};
x_range=-0.5:0.01:0.5;
y_range=x_range;
[X,Y]=meshgrid(x_range,y_range);
XX=x_range;
YY=y_range;
L=length(XX);
for i=1:L
    for j=1:L
        X_temp=zeros(l,1);
        X_temp(1)=XX(i);
        X_temp(7)=YY(j);
        Z(i,j)=c_temp{1}+c_temp{2}'*X_temp+1/2*X_temp'*c_temp{3}*X_temp;
    end
end
% scatter3(dF_training(:,1),dF_training(:,7),dX_training);
% hold on;
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

%% projection test of the result

for i=1:l
xx{i}=-0.5:0.01:0.5;
end

k_test=7;

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
xlabel('$P_{22}$','interpreter','latex');
ylabel('$\mu (\delta y_2)$','interpreter','latex');
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
xlabel('$P_{22}$','interpreter','latex');
ylabel('$\sigma ^2 (\delta y_2)$','interpreter','latex');
l2=legend('True','Predicted');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
set(l2,'box','off');
box on;

%% prediction vs. true

c_temp={c{end,:}};

for i=1:length(Y_test_prcs)
    X_temp=X_test_prcs(i,:)';
    Y_test_prcs_predict(i,1)=c_temp{1}+c_temp{2}'*X_temp+1/2*X_temp'*c_temp{3}*X_temp;
end

scatter(Y_test_prcs,Y_test_prcs_predict,'.');
hold on;
% plot(dX_training,dX_training,'linewidth',2);
plot([-1.5,1.5],[-1.5,1.5],'linewidth',2);
xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
xlabel('True');
ylabel('Predicted');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
box on;
% l1=legend('Data point','y=x');
% set(l1,'box','off');
Err2=abs(Y_test_prcs-Y_test_prcs_predict);
mean(Err2)

%% histogram comparison
c_temp={c{end,:}};

for i=1:length(Y_test_prcs)
    X_temp=X_test_prcs(i,:)';
    Y_test_prcs_predict(i,1)=c_temp{1}+c_temp{2}'*X_temp+1/2*X_temp'*c_temp{3}*X_temp;
end

histogram(Y_test_prcs,'normalization','pdf');
hold on;
histogram(Y_test_prcs_predict,'normalization','pdf');
xlabel('$\delta y_2$','interpreter','latex');
ylabel('P.D.F');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
l1=legend('True','Predicted');
set(l1,'box','off');
box on;

%% Error plot

plot(Err_all_step,'linewidth',2);
xlabel('Step');
ylabel('Loss function W');
set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);