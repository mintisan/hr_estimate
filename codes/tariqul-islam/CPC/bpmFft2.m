function Y = bpmFft2(E,n,Y0,srate)
    E = fft(E,n*1024);
    M = srate*60/4096;
    NEX=1:round(242/M);
    EX = E(NEX);   
    EX = abs(EX);
    [a,I]=max(EX);
    Y=(I(1)-1)*M;
    
    if(Y0>0)
        I = Y0/M+1;
        EX(NEX<I-10 | NEX>I+10)=0; %DELS = 10
        [a,I]=max(EX);
        Y2=(I(1)-1)*M;
        
        if(abs(Y-Y0)<10)
            return;
        else
            Y=Y2;
        end
    end
end