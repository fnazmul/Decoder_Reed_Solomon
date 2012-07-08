function diff = gfdifferentiate(poly)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differentiate a polynomial with respect to x                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:                                                              %
%       poly:   matrix representing polynomial format in field       %
%               i.e. [1 0 2 3] for A + x + A^2x^2 + A^3x^3           %
%               where A is primitive element in the field            %
%Output:                                                             %
%       diff:   differentiated polynomial with respect to x          %
%               in matrix representing polynomial format in field    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calcualte the length of the polynomial
len = length(poly);

for pow = 1:len-1
        %all the even powers are zero
        if mod(pow,2) == 0 
            diff(pow) = -Inf; 
        %coefficient of x^i will be the 
        %coefficient of x^i+1
        else
            diff(pow) = poly(pow+1);
        end
end