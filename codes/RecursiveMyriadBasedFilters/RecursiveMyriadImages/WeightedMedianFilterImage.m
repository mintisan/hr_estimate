function [OutputImage] = WeightedMedianFilterImage(InputImage, Weights, WindowSize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[Row, Col]      = size(InputImage);
OutputImage     = zeros(Row,Col);

for i = 1 : Col - WindowSize + 1
    for j = 1 : Row - WindowSize + 1
        Ii = InputImage(j : j + WindowSize - 1, i : i + WindowSize -1);
        window = Ii(:);
        
        x   = sign(Weights).*window;
        aw  = abs(Weights);
        OutputImage(j + floor(WindowSize/2),i + floor(WindowSize/2)) = wmedian(x, aw);
    end
end

end

