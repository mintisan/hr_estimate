function [Y, r1, r2]=astable(M, N, alpha, beta, dispersion, location)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASTABLE (function)
%       Function for generating alpha stable random vector 'y'.
%       inproved version of a previous FORTRAN function written by
%       John M. Chambers and John Nolan in which a standardized variable 
%       with characteristic exponent, alpha, and symmetric parameter, beta, 
%       is outputted.
%
% INPUT ARGUMENTS:      M -> Length of matrix to be returned
%                       N -> Width of matrix to be returned
%                   alpha -> Characteristic exponent
%                    beta -> Skewness in revised parameterization
%                    disp -> dispersion
%                     loc -> location parameter
%                      
% OUTPUT ARGUMENTS:  y -> Vector containing alpha stable random numbers
%                   
%
%
%
% function [y, r1, r2]=astable(M, N, alpha, beta, dispersion, location)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin == 0) N=1; M=1; end;
if (nargin == 1) N=1; end;
if (nargin <= 2) alpha = 1; end;
if (nargin <= 3) beta = 0; end;   
if (nargin <= 4) dispersion = 1; end;
if (nargin <= 5) location = 0; end;

u=(rand(M,N)-0.5)*pi;
w=-log(rand(M,N));

if (alpha ~= 1)
        X=sin(alpha.*u)./(cos(u)).^(1/alpha).*(cos((1-alpha).*u)./w).^((1-alpha)/alpha);
        Y=X-beta.*tan(pi/2*alpha);
else
        X=sin(alpha*u)./(cos(u)).^(1/alpha).*(cos((1-alpha).*u)./w).^((1-alpha)/alpha);
        Y=X;
end;

Y = Y*dispersion^(1/alpha)+location;
return;


