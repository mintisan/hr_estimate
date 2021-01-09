function [BPM, DOM_PEAK, N_prev]= CPC(sig,BPM_00,srate,p,idnb, N_prev, N_LMS, N_RLS,see,BPM_orig, BPM_DOM)


   %%
    %taking individual array
    PPG = .5*(clean_up(sig(1,:),srate)+clean_up(sig(2,:),srate)); %adding both channel
    XData = clean_up(sig(3,:),srate);
    YData = clean_up(sig(4,:),srate);
    ZData = clean_up(sig(5,:),srate);
    L0 = length(BPM_00);
    

    %Some constants
    M = 4096/(srate*60);
    %%
    %LMS FILTER DATA
    lms_order=26; %length
    mux = 0.0001; %step size
    muy = 0.0001;
    muz = 0.0001;
    
    %%
    %rls filter data
    rls_lambda=0.999;
    rls_order = 55;
    P0 = 10*eye(rls_order);
    
    %%
    %Affine Projection with Recursive Update
    ap_order = 26;
    muapx = 0.008;
    muapy = 0.008;
    muapz = 0.008;
    
    %%
    if L0==0
        Y0 = initial_estimate3(PPG,XData,YData,ZData,srate);
        N_prev = Y0*M+1; 
    else
        Y0=BPM_00(end);
    end
    
    %%

    
    %Perfoming LMS
    
    E1 = zeros(1,1000);
    for i = 1:N_LMS
        E1 = E1 + lms_step(PPG,XData,YData,ZData,mux,muy,muz,lms_order);
    end
    E1 = E1/norm(E1); %Normalizing
    

%     %Performing RLS
    
     E2 = zeros(1,1000);
        for i = 1:N_RLS
           E2 = E2 + rls_step(PPG,XData,YData,ZData,rls_lambda,rls_order,P0);
        end
        E2 = E2/norm(E2); %normalizing
     


    lambda =  1/2;
    
    %EX = lambda* (E1(1:end-3)+ E2(4:end) + E3(4:end));
    EX = lambda* (E1+ E2);

    
    
    %EX = E1+E2;
    EX = clean_up(EX,srate);
    [BPM, DOM_PEAK, N_prev]= bpmFft3(EX,4,Y0, srate, N_prev); %%Applying Trivial BPM Calculation
    
    
    
 
    
    if L0>3
        
        BPM = 0.9*BPM+0.05*BPM_00(end)+0.05*BPM_00(end-1);
        
    end
end