function L=Mean_nD(x_bin,s0,c) %notice: elements in x_bin should be row vector; c{2} is column vector

n=length(s0);

A=zeros(n-1,n-1,n);
B=zeros(n-1,n);
D=zeros(1,n);
mu1_cf=zeros(n-1,n);
s1=zeros(n-1,n-1,n);
C=zeros(n,n,n);
M_cf=zeros(n,n);

M=cell(1,n);
L=cell(1,n);

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
    M2=zeros(n,n,length(x_bin{i}));
    M{i}=M_cf(:,i)*x_bin{i}; %n×m matrix, m is the data amount in xbin{i}
    for j=1:length(x_bin{i})
        M2(:,:,j)=(M{i}(:,j)*M{i}(:,j)');
        L{i}(1,j)=c{1}+c{2}'*M{i}(:,j)+1/2*sum(sum(c{3}.*(C(:,:,i)+M2(:,:,j))));
    end
end