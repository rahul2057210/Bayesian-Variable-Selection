
% Example 4.1

clear all
clc

n = 60;
X = normrnd(0,1,60,5);
X0=ones(60,1);
X4 = X(:,4);
X5 = X(:,5);
%X(:,3)=X(:,5)+0.15*normrnd(0,1,1,60)';
sigma2 = 6.25;
e = normrnd(0,2.5,60,1);

beta4 = 1;
beta5 = 1.2;

Y = X4*beta4 + X5*beta5 + e;
Z=[X0];

Y=(eye(60)-Z*((Z'*Z)\eye(1))*Z')*Y;      % Frisch Waugh Lowell property


X=(eye(60)-Z*((Z'*Z)\eye(1))*Z')*X;       % Frisch Waugh Lowell property


y=[1,1,1,1,1];
bols = (X'*X)\(X'*Y)
Bols=bols;
B1=Bols;
l=((Y-X*Bols)'*(Y-X*Bols))/(55);
sig=l;
k8=(diag((sig)*(X'*X\eye(5)))).^0.5;

Tau=k8';
C=[10,10,10,10,10];
R=eye(5);

V=0;
lamb=0;
jj=Tau.*C;
D=diag(jj); % D matrix



Bstar=Bols';
Ystar=y;
je=1:5;
sigstar=sig;
ZZ=X;
for k = 1:5000
    
    
    A=((1/sig)*(X'*X) + (D*eye(size(D))*D)\eye(size(D)) );
    A=A\eye(size(A));
    
    B=mvnrnd((((1/sig)*A)*(X'*X)*bols)', A,1);
    
    
    je=1:5;
    
    sig=1/(gamrnd((V/2+30),2/((Y-X*B')'*(Y-X*B')),1));
    
    
    n=5;
    kk=n;
    sigstar=[sigstar;sig];
    
    
    
    while kk>1     % generating random permutation
        u1=rand(1);
        l=floor(kk*u1)+1;
        b=l;
        l=kk;
        kk=b;
        je(l)=kk;
        je(kk)=l;
        kk=kk-1;
        
    end
    
    for j=je
        
        D(j,j)=C(j)*Tau(j);
       
  
        a=normpdf(B(j),0,D(j,j)); % copying the sq of diag(D) in identity matrix
        
        
        D(j,j)=Tau(j);
        
        
        b=normpdf(B(j),0,D(j,j)); % copying the sq of diag(D) in identity matrix
        
   
        
        y(j)=binornd(1,a/(a+b),1);
        
        if y(j)==1
            D(j,j)=C(j)*Tau(j);
          
            
            
        else
            D(j,j)=Tau(j);
            
        end
        
    end
    
    %Bstar=[Bstar;B1];
    Ystar=[Ystar;y];
    
end
[ii,jj,kk]=unique(Ystar,'rows');
f=histc(kk,1:numel(jj)); % Frequency
g=sort(f);
r=length(f);
aa=[];
for i=1:5         % print the top 5 models
    j=0;
    k=true;
    while(k)
        j=j+1;
        if f(j)==g(r)
            aa=[aa;ii(j,:) f(j)];
            k=false;
        end
    end
    r=r-1;
end
[1 10]
result=aa


