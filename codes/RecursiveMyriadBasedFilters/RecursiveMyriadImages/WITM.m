% ITM filter; only allow positive weights
function [opt itetime] = WITM (data,w)
n       = length(data);
itetime = 0;
maxw    = max(w);

while(1)
    itetime = itetime+1;
    % get the mu and tau
    mu  = sum(data.*w)/sum(w);        
    tau = (sum((abs(data-mu)).*w))/sum(w);
    % truncate the data    
    bh    = mu+tau;
    bl    = mu-tau;
    data(data>bh) = bh;
    data(data<bl) = bl;
    % stop criteria 1
    wu       = sum(w(data>mu));
    wl       = sum(w(data<=mu));
    if abs(wu-wl)<=maxw
        break
    end
    % stop criteria 2    
    if itetime>=2*n^0.5;
        break
    end           
end

opt  = sum(data.*w)/sum(w);






