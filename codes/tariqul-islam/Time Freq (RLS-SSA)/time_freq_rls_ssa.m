function [BPM, N_prev, d]= time_freq_rls_ssa(sig,BPM_00,srate,i,idnb, N_prev, ABPM ,d)

   
    %%
    %taking individual array
    PPG1 = clean_up(sig(1,:),srate);
    PPG1 = PPG1/norm(PPG1);
    PPG2 = clean_up(sig(2,:),srate);
    PPG2 = PPG2/norm(PPG2);
    PPG = .5*(PPG1+PPG2);
    
    XData = clean_up(sig(3,:),srate);
    XData = XData/norm(XData);
    YData = clean_up(sig(4,:),srate);
    YData = YData/norm(YData);

    ZData = clean_up(sig(5,:),srate);
    ZData = ZData/norm(ZData);
    
    %figure, plot(ACC(1:180));
    L0 = length(BPM_00);

    %Some constants
    M = srate*60/4096;

    
    %%
    %rls filter data
    rls_lambda=0.999;
    rls_order = 55;
    P0 = 10*eye(rls_order);
    
    E2 = rls_step(PPG,XData,YData,ZData,rls_lambda,rls_order,P0, ABPM);
    %%
    if L0<3
        [Y0, N_prev]= initial_estimate(PPG, srate);
        %Nprev = round(Y0*M+1); 
        
    else
        Y0=BPM_00(end);
        %Nprev = round(Y0*M+1);
    end
    
    BPM = Y0;
    
    if L0>3
        
        E2 = rls_step(PPG,XData,YData,ZData,rls_lambda,rls_order,P0, ABPM);
        E3 = pre_ssa(PPG,XData,YData,ZData,4096,N_prev,srate);

       [BS,~]= bpmFft3(E2,4,Y0,srate, N_prev);
    
    if abs(BS- Y0)> 15
        EX = E2;
    else
        EX = E2+E3; %Adding
    end
    
        
    [BPM N_prev]= bpmFft3(EX,4,Y0,srate, N_prev); %%Applying Trivial BPM Calculation
        BPM = 0.9*BPM+0.05*BPM_00(end)+0.05*BPM_00(end-1);
        if((BPM-BPM_00(end))>=5)
         BPM=BPM_00(end)+5;
         N_prev = BPM/M+1;
        elseif((BPM-BPM_00(end))<=-3 )
         BPM=BPM_00(end)-3;
         N_prev = BPM/M+1;
        end
    end

end
