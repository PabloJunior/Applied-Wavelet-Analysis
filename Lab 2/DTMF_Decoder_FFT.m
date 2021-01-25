function DTMF_Decoder_FFT ()
% Распознавание звукового сигнала
FFT('s1.wav')
FFT('s2.wav')
FFT('s3.wav')
end



function FFT(filename)
 
[data, freq] = audioread(filename);
 
time=(1:length(data)) / freq;   
plot(time, data);
 
answer_symbols = [];
n = 1600;
 
for i=1:length(data)/n
    start = (i - 1) * n + 1;
    finish = i * n;
    y = data(start:finish);
    data_fft = unique(abs(fft(y, n)), 'stable');
    if data_fft~=0  %Если звука вообще нет, то ничего не делаем
        all_freqs = (freq / n) * (0:n-1);
        [~, ind] = sort(data_fft, 'descend');
        freqs_for_detect = all_freqs(ind(1:2));
        answer_symbols = [answer_symbols, detect(freqs_for_detect)];
    end
end
    
disp('File:');
disp(filename);
disp('Decoded:');
disp(answer_symbols);
disp('Size:');
disp(length(answer_symbols));
end
 
 
function c = detect(freqs)
 
required_freqs = [697 770 852 941 1209 1336 1477 1633];
[~, ind] = min(abs(required_freqs - freqs(1)));
hifreq = required_freqs(ind);
[~, ind] = min(abs(required_freqs - freqs(2)));
lofreq = required_freqs(ind);
 
if hifreq < lofreq
    tmp = hifreq;
    hifreq = lofreq;
    lofreq = tmp;
end
 
c='('; %bad symbol
 
if lofreq==697
    if hifreq==1209
        c='1';
    end
    
    if hifreq==1336
        c='2';
    end
    
    if hifreq==1477
        c='3';
    end
    
    if hifreq==1633
        c='A';
    end   
end
 
 
if lofreq==770
    if hifreq==1209
        c='4';
    end
    
    if hifreq==1336
        c='5';
    end
    
    if hifreq==1477
        c='6';
    end
    
    if hifreq==1633
        c='B';
    end   
end
            
if lofreq==852
    if hifreq==1209
        c='7';
    end
    
    if hifreq==1336
        c='8';
    end
    
    if hifreq==1477
        c='9';
    end
    
    if hifreq==1633
        c='C';
    end   
end
 
if lofreq==941
    if hifreq==1209
        c='*';
    end
    
    if hifreq==1336
        c='0';
    end
    
    if hifreq==1477
        c='#';
    end
    
    if hifreq==1633
        c='D';
    end   
end
 
end


