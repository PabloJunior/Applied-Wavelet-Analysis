function Cooley_Tukey_Factorization ()
% Факторизация Кули-Тьюки
clc;clear all;
n=8;
 
el=w(n);
F=zeros(n);
for i=1:n
    for j=1:n
        F(i,j)=(el^(i-1))^(j-1);
    end
end

disp('Matrix F_8:');
disp(F);
 
t=log2(n);
A=eye(n);
B=eye(n);
 
for q=1:t
    L=2^q;    
    sub=eye(L);
 
    even=sub(1:2:end,:);
    odd=sub(2:2:end,:);
    PTL=[even;odd]; %PTL
    
    IL2=eye(L/2);%IL2
    
    r=w(L);
    lya=zeros(1,L/2);
    for i=1:L/2
        lya(i)=r^(i-1);
    end
    
    lyambda=diag(lya); %lyambda
    
    BL=[IL2,lyambda;IL2,(-1)*lyambda;]; %BL
    
    Aq=zeros(n);
    RTL=zeros(n);
    p=length(BL);
    g=n/p;    
    for k=1:g
        for i=1:p
            for j=1:p
                Aq((k-1)*p+i,(k-1)*p+j)=BL(i,j);
                RTL((k-1)*p+i,(k-1)*p+j)=PTL(i,j);
            end
        end
    end   
  
    A=Aq*A;
    fprintf('Matrix A_%.f: \n',q);
    disp(Aq);
    B=B*RTL; %RTL   
end

disp('Matrix P_8T:');
disp(B);

disp('Matrix O=F_8-A_3*A_2*A_1*P_8T:');
O=F-A*B;
disp(O)
end
 
 
function w = w(n)
pi=3.1415;
i=sqrt(-1);
w=exp(-2*pi*i/n);
end
