function [Y, DOM, N_prev]= bpmFft3(E,n,Y0,srate, N_prev1)
    E = fft(E,n*1024);
    M = srate*60/4096;
    NEX=1:round(242/M);
    EX = E(NEX);   
    EX = abs(EX);
    [~,I]=max(EX);
    
    Y=(I(1)-1)*M;
    DOM = Y;
    N_prev  = I(1);
    
    [a,ai]=findpeaks(EX);
     p = (a>.25*max(EX));
     if(sum(p)>0)
         l = 0;
     else 
         l=1;
         return;
     end
    
    
    DEL = 10;
    
    if(Y0>0)
        I = N_prev1;
        EX(NEX<I-DEL | NEX>I+DEL)=0; 
        [a,ai]=findpeaks(EX);
        z = ai>I-DEL & ai<I+DEL; %DELS = 10
        a = a(z);
        ai = ai(z);
        [~,i] = max(a);
        I = ai(i);
        if isempty(I)
            [~,I] = max(EX);
        end
        Y2=(I(1)-1)*M;
        
        if(abs(Y-Y0)<10)
            return;
        else
            N_prev  = I(1);
            Y=Y2;
        end
    end
end