%calculate the network deformation given input F
function [dX,Energy_in_all_steps]=calc_def(X,Node_fixed,l0,k_spring,F)
N=length(X);
dX=0.4*(rand(N,2)-0.5)/4;%zeros(size(X));
dX(Node_fixed,:)=0;

eps=1e-6;
alpha=0.01/5;
k=1;
imax=1500;
Err=inf;
ddX=zeros(N,2);
Energy_in_all_steps=zeros(1,imax);

while ((Err>eps) & (k<imax))
    E0_total=Elastic_Energy(X+dX,k_spring,l0)-sum(F(:).*dX(:));
    Energy_in_all_steps(k)=E0_total;
    for i=1:N
        Norm_temp=vecnorm(X(i,:)-X+dX(i,:)-dX,2,2);
        Norm_temp(Norm_temp==0)=1;
        Grad=(sum(k_spring(:,i).*(Norm_temp-l0(:,i)).*(X(i,:)-X+dX(i,:)-dX)./Norm_temp)-F(i,:));
%       ddX(i,:)=-alpha*(sum(k_spring(:,i).*(Norm_temp-l0(:,i)).*(X(i,:)-X+dX(i,:)-dX)./Norm_temp)-F(i,:));
        ddX(i,:)=-alpha*Grad/norm(Grad);
        ddX(Node_fixed,:)=zeros(2,2);
    end
%     ddX
    dX=dX+ddX;
    k=k+1;
    E_total=Elastic_Energy(X+dX,k_spring,l0)-sum(F(:).*dX(:));
    Err=abs(E_total-E0_total);
end

end