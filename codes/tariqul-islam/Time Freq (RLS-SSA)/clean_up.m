function PPG = clean_up(PPG1,srate)
%
%PPG = clean_up(PPG1,srate)
%
%Bandpass Filter: 0.4 to 3.5Hz
%
    ff = fft(PPG1,1024);
    n = 0:1023;
    np = n*2*pi/1024;
    
    ff(np<2*pi*0.4/srate | np>2*pi*(1-0.4/srate)) = 0;
    ff(np>2*pi*3.5/srate & np<2*pi*(1-3.5/srate)) = 0;
    PPG = ifft(ff);
    PPG = real(PPG(1:1000));
end