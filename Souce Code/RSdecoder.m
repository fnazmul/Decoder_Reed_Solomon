function decoded = RSdecoder(recWord)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Decode a (255,239) Reed Solomon code using the Euclidean algorithm  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:                                                              %
%       recWord:  a received codeword in matrix representing         %
%                 polynomial format in field                         %
%                 i.e. [1 0 2 3] for A + x + A^2x^2 + A^3x^3         %
%                 where A is primitive element in the field          %
%Output:                                                             %
%       decoded:  a decoded codeword in matrix representing          %
%                 polynomial format in field                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('Receiving the code..');
    %recWord = [7 14 12 14 9 10 14 6 11 1 5 13 11 14 9];
    
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

    %convert any negative value to -Inf
    for i = 1:n
       if (recWord(i) < 0)
           recWord(i) = -Inf;
       end
    end
    
    
    disp('Decoding the code. Please wait...');
    %calculate syndromes
    S = [];
    %to get S(point) evaluating received polynomial at alpha^point
    for point = 1:twoT
        S(point)= gfpolyval(recWord,point,n,field);
    end
    
    %check if there is no error
    emptyPoly(1:twoT) = -Inf ;
    if(isequal(S, emptyPoly))
        decoded = recWord;
    else
    
        %generate the polynomial S(x)
        for i = 1:twoT
            Sx(i) = S(twoT-i+1);
        end

        %generate a polynomial for x^2t
        xPow2t(1: twoT)= -Inf;
        xPow2t(twoT+1) = 0;

        %run the Euclid alg on x^2t, S(x)
        [gj, rj] = gfeuclidRS(xPow2t, Sx, t, field);

        %find the roots of gj
        roots = gfallroots(gj,n,field);
        
        %differentiate gj
        gjDiff = gfdifferentiate(gj);
        
        %check the number of roots equals the degree of gj
        if not(length(roots)==(length(gj)-1))
            disp(sprintf('   There must have occurred more than %d errors.',t));
            decoded = recWord;
            %return
        else
            %find the error polynomial
            e(1:n) = -Inf;        
            for r = 1 : length(roots)
                %(B^i)^-(2t + 1)
                powB = mod((roots(r)*mod(-(twoT + 1),n)),n);
                %rj(B^i)
                val1 = gfpolyval(rj ,roots(r) ,n ,field);
                %gj'(B^i)
                val2 = gfpolyval(gjDiff ,roots(r) ,n ,field);
                %rj(B^i)/gj'(B^i)
                division = gfdeconv(val1, val2, field);
                %(B^i)^-(2t + 1) * (rj(B^i)/gj'(B^i))
                %v = gfconv(powB,division,field);
                e(roots(r)+1) = gfconv(powB,division,field);
            end

            %calculate the decoded word c = r+e
            decoded = gfadd(recWord,e, field);
                       
        end
    end

    %check if its a codeword
    [quot,remd] = gfdeconv(decoded,g,field);
    if~(remd == -Inf)
        disp('Decoded word is not a valid codeword!');
    end 
end


