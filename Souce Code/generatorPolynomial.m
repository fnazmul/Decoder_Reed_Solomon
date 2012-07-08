function g = generatorPolynomial(twoT, field)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the generator polynomial of the twoT/2 error correcting  %
% Reed Solomon code from the product of (x+alpha^i), where i=1..twoT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:                                                              %
%       twoT:   if the code is t error correcting then 2*t           %
%       field:  list of all elements in the field p^m                %
%               generated using gftuple                              %
%Output:                                                             %    
%       g:      the generator polynomial in matrix                   %
%               representing polynomail format                       %
%               i.e. [1 0 1 1 1 0 0 0 1]  for field 2^8              %    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %generate the generator polynomial
    
    %represents (alpha + X )
    gen = [1 0]; 
    temp(1) = gen(1);

    for i = 1:twoT-1
        %temp represents (alpha^i + X )
        temp(1) = gfmul(temp(1),1,field);
        temp(2) = 0;
        %muliplies ((alpha + X )...(alpha^i-1 + X )) and (alpha^i + X )
        gen = gfconv(gen,temp,field);
    end
    g = gen;

end