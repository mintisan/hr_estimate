function [Y, N_prev]= bpmFft3(E,n,Y0,srate, N_prev)
% Y = bpmFft2(E,n)
% E is the noise removed PPG Signal of length 1000
% n is the factor by which FFT points will be calculated
% n should be equal to 4 always in current implementation
% srate is sampling rate
% This computes BPM using FFT MEthod
% Written for SP CUP 2015
    E = fft(E,n*1024);
    M = srate*60/4096;
    NEX=1:round(242/M);
    EX = E(NEX);   
    EX = abs(EX);
    [~,I2]=max(EX);
    Y=(I2-1)*M;
    DEL = 10;
    T0 = 10;
    
    [~,I]=max(EX);
    [a,~] = findpeaks(EX);
     p = (a>.25*max(EX));
     if(sum(p)<=1)
             
        %I = HR_PPG(PPGX, I, M, 2);
        Y = (I(1)-1)*M;
        N_prev = I;
        return;
     end
    
    if(Y0>0)
        I = N_prev;
        EX(NEX<I-DEL | NEX>I+DEL)=0; 
        [a,ai]=findpeaks(EX);
        z = ai>I-DEL & ai<I+DEL; %DELS = 10
        a = a(z);
        ai = ai(z);
        [a2,i] = max(a);
        I = ai(i);
        if isempty(I)
            [a,I] = max(EX);
        end
        Y2=(I(1)-1)*M;
        
        if(abs(Y-Y0)<T0)
            N_prev = I2;
            return;
        else
            N_prev = I(1);
            Y=Y2;
        end
    end
end