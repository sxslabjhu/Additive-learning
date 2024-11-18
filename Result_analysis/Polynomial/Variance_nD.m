function V=Variance_nD(x_bin,s0,c) %notice: elements in x_bin should be row vector

n=length(s0);
L=Mean_nD(x_bin,s0,c);

A=zeros(n-1,n-1,n);
B=zeros(n-1,n);
D=zeros(1,n);
mu1_cf=zeros(n-1,n);
s1=zeros(n-1,n-1,n);
C=zeros(n,n,n);
M_cf=zeros(n,n);

M=cell(1,n);
V=cell(1,n);

for i=1:n
    st=s0;
    col=st(:,i);
    st(i,:)=[];
    st(:,i)=[];
    A(:,:,i)=st;
    D(i)=col(i);
    col(i)=[];
    B(:,i)=col;
    mu1_cf(:,i)=B(:,i)/D(i);
    s1(:,:,i)=A(:,:,i)-B(:,i)*B(:,i)'/D(i);
    M_cf(:,i)=[mu1_cf(1:i-1,i);1;mu1_cf(i:end,i)];
    C(1:(i-1),1:(i-1),i)=s1(1:(i-1),1:(i-1),i);
    C((i+1):end,(i+1):end,i)=s1(i:end,i:end,i);
    C(1:(i-1),(i+1):end,i)=s1(1:(i-1),i:end,i);
    C((i+1):end,1:(i-1),i)=s1(i:end,1:(i-1),i);
end

for i=1:length(x_bin) %notice: elements in x_bin should be row vector
    M{i}=M_cf(:,i)*x_bin{i}; %n×m matrix, m is the data amount in xbin{i}
    V{i}=zeros(1,length(x_bin{i}));
    for j=1:length(x_bin{i})
        E1{i}{j}=M{i}(:,j);
        E2{i}{j}=C(:,:,i)+M{i}(:,j)*M{i}(:,j)';
        V{i}(j)=c{1}^2+2*c{1}*(c{2}'*E1{i}{1,j}+0.5*sum(sum(c{3}.*E2{i}{1,j})))+c{2}'*E2{i}{1,j}*c{2}-L{i}(j)^2;
        for k1=1:n
            E3{i}{j}(k1,:,:)=M{i}(k1,j)*C(:,:,i)+M{i}(:,j)*C(k1,:,i)+C(k1,:,i)'*M{i}(:,j)'+M{i}(k1,j)*M{i}(:,j)*M{i}(:,j)';
            for k2=1:n
                E4{i}{j}{k1,k2}=C(k1,k2,i)*C(:,:,i)+C(k1,:,i)'*C(k2,:,i)+C(k2,:,i)'*C(k1,:,i)+...
                                        M{i}(k1,j)*M{i}(k2,j)*C(:,:,i)+M{i}(k1,j)*M{i}(:,j)*C(k2,:,i)+...
                                        M{i}(k1,j)*C(k2,:,i)'*M{i}(:,j)'+M{i}(k2,j)*M{i}(:,j)*C(k1,:,i)+...
                                        M{i}(k2,j)*C(k1,:,i)'*M{i}(:,j)'+M{i}(:,j)*M{i}(:,j)'*C(k1,k2,i)+...
                                        M{i}(k1,j)*M{i}(k2,j)*M{i}(:,j)*M{i}(:,j)';
                 V{i}(j)=V{i}(j)+0.5*(c{2}(k2)*c{3}(k1,:)+c{2}(k1)*c{3}(k2,:))*E3{i}{1,j}(k1,:,k2)'+0.25*c{3}(k1,:)*E4{i}{1,j}{k1,k2}*c{3}(k2,:)';
            end
        end
    end
end