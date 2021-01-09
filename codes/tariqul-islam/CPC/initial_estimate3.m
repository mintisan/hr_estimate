function y = initial_estimate3(PPGCX,X,Y,Z,srate)
        P0=10*eye(55);
        ha=adaptfilt.rls(55,0.999,P0);
        [~,E1] = filter(ha,X,clean_up(PPGCX,srate));
        ha=adaptfilt.rls(55,0.999,P0);
        [~,E1] = filter(ha,Y,E1);
        ha=adaptfilt.rls(55,0.999,P0);
        [~,E1] = filter(ha,Z,E1);
        
        y=bpmFft2(E1,4,0,srate);    
end
