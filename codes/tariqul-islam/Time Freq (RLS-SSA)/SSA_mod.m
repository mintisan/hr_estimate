function yn =  SSA_mod( L, K, nheart, loca, nfft, tol, tol_s, tol_ds)

% in this code if a certain percentage(tol) of dominant frequencies of
% a time series are in Facc, it won't be added


%   nheart = 1000 prefiltered ppg sample
%   L = number of rows in hankel matrix
%   K = number of colums in hankel matrix
%   loca = location of dominant frequencies without locr-10:locr+10
%   nfft = number of points in FFT
%   tol = percentage of dominant frequeny in Facc to be MA component 

%   yn = SSA filtered ppg signal
    
    N = 1000 ;
    yn = zeros(1,N) ;
    %d = min(L,K) ;
    
    % Embedding Step
    XX = hankel(nheart(1:L),nheart(L:end)) ;
    
    % SVD Step and reconstruction
    [U, S, V] = svd(XX) ;
          s = diag(S) ;
%     s = s(s>0.01*S(1,1)) ;
    s = s(s >= tol_s*s) ;

    ds = s(1:end-1) - s(2:end);
    ds = [0; ds];

    [~ , loc_ds] = findpeaks(ds);
    loc_ds=loc_ds-1;
    
    for k1 = 1:length(loc_ds)
        Yip = zeros(size(XX));
        if ds(loc_ds(k1)+1) >= tol_ds*max(ds) && ds(loc_ds(k1)) <= 0.5*ds(loc_ds(k1)+1)
%             flag1 = 1;
%         else
%             flag1 = 0 ;
%         end
%             
        if(k1 == 1 ) 
            start_k1 = 1;
            end_k1 = loc_ds(1);
%         elseif(k1==length(loc_ds))
%             start_k1 = end_k1+1;
%             end_k1 = d;
        else
            start_k1 = loc_ds(k1)-1 ;
            end_k1 = loc_ds(k1) ;
        end
            
        for j = start_k1:end_k1
            Yi = S(j,j)*U(:,j)*V(:,j)' ;
            Yip = Yip + Yi ;
         end
            fYi = fliplr(Yip) ;

            ydn = zeros(1,N) ;
            k = K-1 ;
            for j1 = 1:N 
                ydn(j1) = mean(diag(fYi,k)) ;
                k = k-1 ;
            end
            Ydn = abs(fft(ydn,nfft)).^2 ;  
 %           M_Ydn = max(abs(Ydn)) ;
 
            Ydn = [0 Ydn(1:end/2)] ;
            [~,locm] = max(Ydn(2:end)) ;
            
 %           locrr = [locr-Del:locr+Del 2*(locr-1)+1-Dels:2*(locr-1)+1+Dels] ;
%            for jj = 1:length(locrr)
%                      if( any(locrr == locm) )
%                          flag2 = 1 ;
%                      else
%                          flag2 = 0 ;
%                      end
% %            end
%             if flag1 == 1 || flag2 == 1

            [~,loc_pk] = findpeaks(Ydn) ;
            loc = loc_pk(Ydn(loc_pk) >= 0.5*max(Ydn) );
            loc = loc-1 ;
      
            l_loc = length(loc) ;
            count = 0 ;
            for j2 = 1:l_loc
                if( any(loc(j2)==loca) )
%                     flag = 1 ;
%                     break ;
                    count = count+1 ;
                end
            end
            
            if tol == 0
               if count == 0
                   yn = yn+ydn ;
               end
            else    
                if(count < tol*l_loc)
                    yn = yn+ydn ;
                end
            end
        end     
   
  
    end

end
