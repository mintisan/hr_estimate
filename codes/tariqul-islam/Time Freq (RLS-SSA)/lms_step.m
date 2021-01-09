function y = lms_step(PPGCX,X,Y,Z,mux,muy,muz,filter_length)
%%
    y = lms_noise_canceller(PPGCX,X,mux,filter_length);
    y = lms_noise_canceller(y,Y,muy,filter_length);
    y = lms_noise_canceller(y,Z,muz,filter_length);
    
end

function y = lms_noise_canceller(ppg_data,acc_channel,mu_step,filter_length)
%%
%Function uses Adaptive LMS Filter in NOISE CANCELLING SCHEME
%
%
%       INPUTS
%
%ppg_data:      PPG Data of length 1000
%acc_data:      Accelerometer data of length 1000, time synchronized with
%               PPG Data, only one channel
%mu_step:       Step Size of the filter, usually of order 1e-3
%               Defaults to 0.003
%filter_length: Length of the filter. Defaults to 20
%
%       OUTPUTS
%y:             Noise Removed data
%
%       Example:
%
%   y=lms_noise_canceller(PPGC(z),XData(z),0.003,20);
%   z is a preferred time window of length 1000
%

%%
%LMS Filter Data
if nargin<4
    filter_length=20;
end
if nargin<3
    mu_step = 0.003;
end
%Handle of the filter/filter object
ha = adaptfilt.lms(filter_length,mu_step);

%%
%Filtering using adaptive noise canceller
[a,y] = filter(ha,acc_channel,ppg_data);

end