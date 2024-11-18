%calculate the elastic energy
function E=Elastic_Energy(X,k,l0)
E=0;
[N,~]=size(X);
for i=1:N
    E=E+1/4*sum(k(:,i).*(vecnorm(X-X(i,:),2,2)-l0(:,i)).^2); %it's 1/4 because there are repeated sum
end
end