function code = RSencoder(info)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Encode a (255,239) Reed Solomon code using the systematic encoding  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:                                                              %
%       info:     information bits represented in field              %
%Output:                                                             %
%       code:     an encoded codeword in matrix representing         %
%                 polynomial format in field                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('Receiving information bits..');
    %prime
    p = 2;
    % Reed Solomon code over GF(2^m)
    m = 8; %8
    % Length of codeword
    n = 2^m -1; 

    % number of errors can be corrected
    t = 8; %8;
    twoT = 2*t;

    % Dimension of codeword
    k = n - twoT; %239

    %default primitive polynomial
    prim_poly = gfprimdf(m,p);
    
    %generate a list of elements of GF(2^m)
    field = gftuple([-1:p^m-2]',m,p);
    
    %generate the generator polynomial
    g = generatorPolynomial(twoT, field);
    
    disp('Encoding information bits..');
    
    %Systematic encryption
    %parity bits are calculated by (X^(n-k).i(X)) / g(X)
    %codeword = (X^(n-k).i(X)) + parity bits
        
    %a polynomial representing X^(n-k)
    shiftPoly(1:n-k) = -Inf;
    shiftPoly(n-k+1) = 0;
    %multiplying it with the info to shift
    shiftInfo = gfconv(info,shiftPoly,field);

    %divide shifted info by g(x)
    [quot, parity] = gfdeconv(shiftInfo, g, field);

    %if padding is needed for the parity bits
    temp = -Inf;
    while length(parity) < n-k
        parity = [parity temp];
    end

    %concatenate the parity bits to the data
    code = [parity info];

end