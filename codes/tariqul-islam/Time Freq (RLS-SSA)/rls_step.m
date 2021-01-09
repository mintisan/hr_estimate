function y = rls_step(PPGCX,X,Y,Z,lambda,filter_length,P0, ABPM)
%%
    y = rls_noise_canceller(PPGCX,X,lambda,filter_length,P0);
    
%     figure();
%     %subplot(4,1,1)
%     xlabel('(a)');
%     subplot(4,1,2), peri(y, 4, ABPM)
%     set(gca,'XMinorTick','on','YMinorTick','off');
%     set(gca,'yticklabel',[]);
%     xlabel('(b)');
%     title('Periodogram of PPG1 (after RLS Block 1)');
    
     y = rls_noise_canceller(y,Y,lambda,filter_length,P0);
%     subplot(4,1,3), peri(y, 4, ABPM)
%     set(gca,'XMinorTick','on','YMinorTick','off');
%     set(gca,'yticklabel',[]);
%     xlabel('(c)');
%     title('Periodogram of PPG2 (after RLS Block 2)');
    
    y = rls_noise_canceller(y,Z,lambda,filter_length,P0);
    
%     subplot(4,1,4), peri(y, 4, ABPM)
%     set(gca,'XMinorTick','on','YMinorTick','off');
%     set(gca,'yticklabel',[]);
%     xlabel('(d)');
%     title('Periodogram of PPG Final Signal (after RLS Block 3)');
    

    
end

function y = rls_noise_canceller(ppg_data,acc_channel,lambda,order,P0)
%%
%RLS Filter Data

%Handle of the filter/filter object
ha = dsp.RLSFilter('Length',order,'ForgettingFactor',lambda,'InitialInverseCovariance',P0);%,randn(1,order)/10);

%%
%Filtering using adaptive noise canceller
[a,y] = ha(acc_channel,ppg_data);

end