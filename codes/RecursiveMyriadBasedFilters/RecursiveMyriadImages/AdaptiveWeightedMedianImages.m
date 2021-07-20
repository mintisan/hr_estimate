function [Weights, e] = AdaptiveWeightedMedianImages(InputImage, DesiredImage, Weights, WindowSize, u)
%% Adaptive Weighted Median on Images
%
% Inputs:   InputImage      Noisy image
%           DesiredImage    Clean image
%           Weights         Initial weight vector
%           WindowSize      Window size
%           u               Step size
%
% Outputs:  Weihgts     `   Output weight vector
%           e               Error 
%
% Created by Juan Marcos Ramirez (December, 2015)
% Please report bugs to Juan Ramirez (juanra@ula.ve)
% Universidad de Los Andes, Facultad de Ingenieria, Escuela de Electrica
% Merida, Venezuela
%
% References    (1) ARCE, Gonzalo R. A general weighted median filter 
%                   structure admitting negative weights. Signal 
%                   Processing, IEEE Transactions on, 1998, vol. 46, 
%                   no 12, p. 3195-3205.

%%
[Row, Col]  = size(InputImage);
N1          = Col - WindowSize + 1;
N2          = Row - WindowSize + 1;
e               = zeros(N1*N2,1);
m               = 1;
for i = 1 : N1
    for j = 1 : N2
        Ii      = InputImage(j : j + WindowSize - 1, i : i + WindowSize -1);
        window  = Ii(:);
        d       = DesiredImage(j + floor(WindowSize/2),i + floor(WindowSize/2));
        
        x   = sign(Weights).*window;
        aw  = abs(Weights);
        y   = wmedian(x, aw);
        
        e(m)    = d - y;
        Weights = Weights + u * e(m) * sign(Weights) .* sign(x - y);
        m       = m + 1;
    end
end