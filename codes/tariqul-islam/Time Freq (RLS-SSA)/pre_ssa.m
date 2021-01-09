function yn = pre_ssa(nheart, xdata, ydata, zdata, nfft, locr,srate)
   N = length(nheart);
   L = 400;
   K = N-L+1;
   tol_s = 0.1 ;
   tol_ds = 0.01 ;
   Del = 5 ;
        X = abs(fft(xdata,nfft)).^2 ;
%         X = X.*index ;
        Y = abs(fft(ydata,nfft)).^2 ;
%         Y = Y.*index ;
        Z = abs(fft(zdata,nfft)).^2 ;
%         Z = Z.*index ;
        
%         X = [0 X(1:end/2)] ;
%         Y = [0 Y(1:end/2)] ;
%         Z = [0 Z(1:end/2)] ;
        
        [~,loc_pk] = findpeaks(X) ;
        locx = loc_pk(X(loc_pk) >= 0.5*max(X) );
        
        [~,loc_pk] = findpeaks(Y) ;
        locy = loc_pk(Y(loc_pk) >= 0.5*max(Y) );
        
        [~,loc_pk] = findpeaks(Z) ;
        locz = loc_pk(Z(loc_pk) >= 0.5*max(Z) );
        loca = [locx locy locz] ;

%         locx = find(X >= 0.25*max(X)) ;
%         locy = find(Y >= 0.25*max(Y)) ;
%         locz = find(Z >= 0.25*max(Z)) ;
%         loca = [locx locy locz] ;
% 
%          loca = loca-1 ;
        
         locrr = [locr-Del:locr+Del 2*(locr-1)+1-Del:2*(locr-1)+1+Del] ;
         for jj = 1:length(locrr)
             loca = loca(loca ~= locrr(jj)) ;
         end
    
    
        yn = clean_up(SSA_mod(L, K, clean_up(nheart,srate),loca, nfft, 1, tol_s, tol_ds),srate) ;    
%       yn = SSA_tn(L, K, nheart, loca, nfft) ;
end


