function Y=MyNN_Prediction(X,W)

N_HL=length(W)-1;
X_temp=X;

for k=1:N_HL
    [M,~]=size(X_temp);
    X_temp=[ones(M,1),X_temp];
    X_temp=X_temp*W{k};
    X_temp=Act(X_temp);
end
[M,~]=size(X_temp);
X_temp=[ones(M,1),X_temp];
Y=X_temp*W{N_HL+1};
end