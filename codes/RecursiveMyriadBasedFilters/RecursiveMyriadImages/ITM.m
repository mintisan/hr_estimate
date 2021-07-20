function [output] = ITM (data,opt_type)

%% ITM filter 

%input: data      //input signal
%       opt_type  // type of the filter ouput: ITM1, ITM2

%output: output // filter output

% created by Miao Zhenwei.
% Nanyang Technological Univerisity, ROSE center
% Webpage: https://sites.google.com/site/miaozhenwei/


% Please report bugs and/or send comments to Miao Zhenwei.
% zwmiao@ntu.edu.sg

% Reference 1: Z. W. Miao and X. D. Jiang, Further Properties and a Fast Realization
%             of the Iterative Truncated Arithmetic Mean Filter IEEE Transactions on Circuits and Systems-II, 
%             vol. 59, no. 11, pp. 810-814, November 2012.
% Reference 2:
%             X.D. Jiang, "Iterative Truncated Arithmetic Mean Filter And Its Properties," IEEE Transactions 
%             on Image Processing, vol. 21, no. 4, pp. 1537-1547, April 2012.
%%

n            = length(data);
c2           = 2*n^0.5;
c3           = (n-n^0.5)/2;
c4           = n^0.5;

iteration_times  = 0;
S3minus1         = -1;

while(1)
    iteration_times = iteration_times+1;
    
    mu       = mean(data);
    tao3      = mean(abs(data-mu));
    
    bh       = mu+tao3;
    bl       = mu-tao3;
    
    nh       = sum(data>mu);
    nl       = n-nh;
    
    phtao   = data>mu+tao3;
    pltao   = data<mu-tao3;
    ntaoh   = sum(phtao);
    ntaol   = sum(pltao);
    
    data(phtao) = mu+tao3;
    data(pltao) = mu-tao3;
    
    % stop criteria 1
    ditn    = abs(nh-nl);
    if ditn<=1
        break
    end
    % stop criteria 2
    if iteration_times>=c2
        break
    end
    % stop criteria 3
    S3 = abs(ntaoh-ntaol);
    if S3>=c3
        break
    end
    % stop criteria 4
    if (S3>=c4)&&(S3 == S3minus1)
        break
    end
    S3minus1 = S3;
end

% output 1
if opt_type==1
    output  = mean(data);
else if opt_type==2 % output 2

        idxr    = data>bl&data<bh;
        if sum(idxr)>(n/4)
            output = mean(data(idxr));
        else
            output  = mean(data);
        end
    end
end
