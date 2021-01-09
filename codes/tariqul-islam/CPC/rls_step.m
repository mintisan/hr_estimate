function y = rls_step(PPGCX,X,Y,Z,lambda,filter_length,P0)
%%
    y = rls_noise_canceller(PPGCX,X,lambda,filter_length,P0);
    y = rls_noise_canceller(y,Y,lambda,filter_length,P0);
    y = rls_noise_canceller(y,Z,lambda,filter_length,P0);
    
end

function y = rls_noise_canceller(ppg_data,acc_channel,lambda,order,P0)
%%
%RLS Filter Data

%Handle of the filter/filter object
ha = adaptfilt.rls(order,lambda,P0);%,randn(1,order)/10);

%%
%Filtering using adaptive noise canceller
[a,y] = filter(ha,acc_channel,ppg_data);

end