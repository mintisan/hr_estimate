function [OutputImage] = ITMFilterImage(InputImage, WindowSize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[Row, Col]      = size(InputImage);
OutputImage     = zeros(Row,Col);

for i = 1 : Col - WindowSize + 1
    for j = 1 : Row - WindowSize + 1
        Ii = InputImage(j : j + WindowSize - 1, i : i + WindowSize -1);
        OutputImage(j + floor(WindowSize/2),i + floor(WindowSize/2)) = ITM(Ii(:), 1);
    end
end

end

