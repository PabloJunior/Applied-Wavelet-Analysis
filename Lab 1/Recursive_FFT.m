function Recursive_FFT()
% Рекурсивный алгоритм БПФ
clc; clear all;
T1=[]; T2=[]; K=[];
B=[]; E=[];
 
for k=1:16
    K=[K,k];
    n=2^k;
    x = randi(100,n,1);
    
    tic;
    y=FFT(x);
    t1=toc;
    T1 = [T1, t1];  
         
    tic;
    z=truefft(x);
    t2=toc;
    T2 = [T2, t2]; 
    
    x_=IFFT(y);    
    En=norm(x-x_.');
    E=[E,En];
    
    Bn = norm(y-z.');
    B = [B, Bn];
    
    if k==3
        format bank;
% форматированный вывод, дабы не выводилось 50.0000 + 0.0000i 
        disp('x^8:'); disp(x.');
        disp('x~8:'); disp(x_);   
        disp('y^8:'); disp(y);
        disp('z^8:'); disp(z.');
        
        format;
% восстановление формата вывода по умолчанию для вывода невязок
        disp('e^8:'); disp(En);
        disp('b^8:'); disp(Bn);
    end 
end
 
ty = double(T1);
tz=double(T2);
 
figure('Color','w');
plot( K,ty, 'M*', K, tz, 'ro')
grid on
zoom on
xlabel('Степень вектора, k');
axis([0 17 0 2.5]);
ylabel('Время, c');
title('Точечные графики t_y^n и t_z^n');
legend({'t_y^n','t_z^n'},'Location','northwest','Orientation','vertical')
zoom on
 
disp('------------------------------------------------------');
disp('|  K |     En      |     Bn      |  t_y^n  |  t_n^z  |');
 
for i=1:16
    fprintf('------------------------------------------------------\n'); 
    formatSpec = '|%3.f | %0.5e | %0.5e | %0.5f | %0.5f |\n';
    fprintf(formatSpec,K(i),E(i),B(i),ty(i),tz(i));    
end

end
 
function y = sub(x,p)
N = length(x);
if N == 1
    y = x;
else
    M = N / 2;
    W = exp(p*2i*pi*(0:M-1)/N);
    odd = sub(x(1:2:N),p);
    even = sub(x(2:2:N),p);
    l = odd(1:M) + W.*even(1:M);
    r = odd(1:M) - W.*even(1:M);
    y = [l,r];
end
end
 
function y = FFT(x)
y=sub(x,-1);
end
 
function zn=truefft(x)
zn=fft(x);
end
 
function y = IFFT(x)
N = length(x);
y=1/N*sub(x,1);
end


