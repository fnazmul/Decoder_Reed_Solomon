%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file checks the performance of the (255,239) Reed-Solomon decoder%
% It takes one frame of randomly generated codeword                     %
% Randomly creates some errors in the codeword                          %
% Then uses the decoder to decode the code and                          %
% checks if it is decoded successfully                                  %
% or cannot be decoded or decoded into wrong codeword                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;

%generate a list of elements of GF(2^m)
field = gftuple([-1:2^8-2]',8,2);

%generate the generator polynomial
g = generatorPolynomial(8, field);

%generate random data
info = randint(1,239,[-1 255-1]);

for i = 1:239
    if (info(i) < 0)
        info(i) = -Inf;
    end
end

%encoding information bits
encoded = RSencoder(info); 

disp('Sending the code');
send = encoded;

field = gftuple([-1:2^8-2]',8,2);

%creating random errors
 send(3) = gfadd(encoded(3),randint(1,1,[-1 255-1]),field);
 send(5) = gfadd(encoded(3),randint(1,1,[-1 255-1]),field);
 send(15) = gfadd(encoded(15),randint(1,1,[-1 255-1]),field);
 send(67) = gfadd(encoded(3),randint(1,1,[-1 255-1]),field);
 send(122) = gfadd(encoded(122),randint(1,1,[-1 255-1]),field);
 send(141) = gfadd(encoded(141),randint(1,1,[-1 255-1]),field);
 send(167) = gfadd(encoded(167),randint(1,1,[-1 255-1]),field);
 send(188) = gfadd(encoded(188),randint(1,1,[-1 255-1]),field);
 send(207) = gfadd(encoded(207),randint(1,1,[-1 255-1]),field);
 send(247) = gfadd(encoded(247),randint(1,1,[-1 255-1]),field);
  
 send(79) = gfadd(encoded(79),randint(1,1,[-1 255-1]),field);
 send(139) = gfadd(encoded(139),randint(1,1,[-1 255-1]),field);
 send(193) = gfadd(encoded(193),randint(1,1,[-1 255-1]),field);
 send(235) = gfadd(encoded(235),randint(1,1,[-1 255-1]),field);
%uses decoder to decode
DECODED = RSdecoder(send);

%checks the result of decoding
if (isequal(DECODED,encoded))
    disp('Succesful Decoding')
elseif (isequal(DECODED,send))
    disp('No Change')
else
    disp('Decoding Error')
end
