function [Y, N_prev] = initial_estimate(PPGCX,srate)
    E = fft(PPGCX,4*1024);
    M = srate*60/4096;
    NEX=1:round(242/M);
    EX = E(NEX);   
    EX = abs(EX);
    [a,I]=max(EX);
    N_prev = I(1);
    Y=(N_prev-1)*M;

end
