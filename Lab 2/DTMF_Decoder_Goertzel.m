function DTMF_Decoder_Goertzel ()
% Распознавание звукового сигнала
clear all;close all;clc;
set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Arial Cyr');
set(0,'DefaultTextFontSize',10,'DefaultTextFontName','Arial Cyr');
 

filename='s1.wav';
[data, freq]=audioread(filename);
%sound(wave, freq);       
time=(1:length(data))/freq;
figure('Color','w');
plot(time, data);
title('График звуковых волн файла');
xlabel('Время (с)')
ylabel('Амплитуда')

sample = 0.05; 
t = linspace(0,sample,sample*freq); 
width = length(t); 
height = length(data)/width; 
dataS = reshape(data(1:height*width),width,height);
 
 
fr = [697,770,852,941,1209,1336,1477,1633];
KEYS = [];     ncol = size(dataS,2); 
 
for i=1:4:ncol
    y = [dataS(:, i)',dataS(:, i+1)'];
    n = length(y);
    freq_indices = round(fr/freq*n) + 1;   
    dft_data = goertzel(y,freq_indices);
    t=abs(dft_data);
    %stem(fr,t);
    
    [val,ind] = sort(t,'descend');
    an=fr(ind(1:2));
    [val,ind] = sort(an,'descend');
    KEYS=[KEYS,detect(val)];
end
 
 
formatSpec = 'В файле %s были нажаты следующие клавиши:\n';
fprintf(formatSpec,filename);   
disp(KEYS);
end
  
function c=detect(val)
hifreq=val(1);
lofreq=val(2);
dialArray = ['1' '2' '3' 'A';'4' '5' '6' 'B';
    '7' '8' '9' 'C';'*' '0' '#' 'D'];
rowList = [1209 1336 1477 1633];
colList = [697 770 852 941];
r = find(rowList==hifreq);
k = find(colList==lofreq);
c=dialArray(k,r);
end

