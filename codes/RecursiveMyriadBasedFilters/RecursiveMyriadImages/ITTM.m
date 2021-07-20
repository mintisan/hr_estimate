
function [opt] = ITTM (data,opt_type)
%% ITTM filter
%input: data      //input signal
%       opt_type  // type of the filter ouput: ITTM1, ITTM2, ITTM3

%output: output  // filter output

% created by Miao Zhenwei.
% Nanyang Technological Univerisity, ROSE center
% Webpage: https://sites.google.com/site/miaozhenwei/


% Please report bugs and/or send comments to Miao Zhenwei.
% zwmiao@ntu.edu.sg

%  Reference: Z. W. Miao and X. D. Jiang, "Additive and Exclusive Noise Suppression byIterative Trimmed and Truncated Mean Algorithm,” Signal Processing, vol. 99, pp. 147-158, June, 2014.

% Relevant papers

%             Z. W. Miao and X. D. Jiang, “Weighted Iterative Truncated Mean Filter,” IEEE Transactions on Signal Processing, Vol. 61, no. 16, pp. 4149-4160, August, 2013.
%
%             Z. W. Miao and X. D. Jiang, “Further Properties and a Fast Realization
%             of the Iterative Truncated Arithmetic Mean Filter” IEEE Transactions on Circuits and Systems-II, 
%             vol. 59, no. 11, pp. 810-814, November 2012.
%  
%             X.D. Jiang, "Iterative Truncated Arithmetic Mean Filter And Its Properties," IEEE Transactions 
%             on Image Processing, vol. 21, no. 4, pp. 1537-1547, April 2012.
%%
n     = length(data);
ntauh = 0;
ntaul = 0;
bh    = 0;
bl    = 0;
rn    = n;
itetime = 0;
S3minus1 = -1;

while(1)
    itetime = itetime+1;
    
    
    if(ntauh>=ntaul)
        mu  = (sum(data)+(ntauh-ntaul)*bh)/(n-2*ntaul);
        tau = (sum(abs(data-mu))+(ntauh-ntaul)*(bh-mu))/(n-2*ntaul);
    else
        mu  = (sum(data)+(ntaul-ntauh)*bl)/(n-2*ntauh);
        tau = (sum(abs(data-mu))+(ntaul-ntauh)*(mu-bl))/(n-2*ntauh);
    end
    
    bh    = mu+tau;
    bl    = mu-tau;
    ntauh = sum(data>bh)+ntauh;
    ntaul = sum(data<bl)+ntaul;
    data  = data(data>=bl&data<=bh);
    rn    = length(data);
    nh    = sum(data>mu)+ntauh;
    nl    = sum(data<=mu)+ntaul;
    
    
    ditn    = abs(nh-nl);
    if ditn<=1
        break
    end
    % stop criteria 2
    
    if itetime>=2*n^0.5;
        break
    end
    
    % stop criteria 3
    S3 = abs(ntauh-ntaul);
    if S3>=n^0.5
        break
    end
    %         % stop criteria 4
    if (S3>=(n-n^0.5))&&(S3 == S3minus1)
        break
    end
    S3minus1 = S3;
    
    
    
end


if opt_type == 1 % output 1
    opt  = (sum(data)+ntauh*bh+ntaul*bl)/n;
else
    if opt_type == 2       % output 2
        if rn>(n/4)
            opt = mean(data);
        else
            opt  = (sum(data)+ntauh*bh+ntaul*bl)/n;
        end
    else
        if opt_type == 3       % output 3
            if(ntauh>ntaul)
                if ( 2*ntaul<n*0.75)
                    opt = (sum(data)+(ntauh-ntaul)*bh)/(n-2*ntaul);
                else
                    opt  = (sum(data)+ntauh*bh+ntaul*bl)/n;
                end
            else
                if ( 2*ntauh<n*0.75)
                    opt = (sum(data)+(ntaul-ntauh)*bl)/(n-2*ntauh);
                else
                    opt  = (sum(data)+ntauh*bh+ntaul*bl)/n;
                end
            end
        end
    end
end





