%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file checks the performance of the (255,239) Reed-Solomon decoder %
% It takes 100 frames of zero codewords                                  %
% Randomly creates same number of errors in each frame                   %   
% Then decodes each frame and checks if it is decoded successfully       %
% or cannot be decoded or decoded into wrong codeword                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc

% taking the parameters
n = 255;
k = 239;
errorNum = 9;

%generate a list of elements of GF(2^m)
field = gftuple([-1:2^8-2]',8,2);

%generate the generator polynomial
g = generatorPolynomial(8, field);

%a zero codeword
allEmpty(1:n) = -Inf;

%generating random errors in each frame of length n
recFrame = randerr(100,n,errorNum);

%for each frame
for(frame = 1:100)
    
    disp(' ');
    disp(sprintf('Processing frame.. %d',frame));
    
    %change the format to field format
    for (i = 1:n)
        if (recFrame(frame,i)== 0)
            recFrame(frame,i) = -Inf;
        else
            recFrame(frame,i) = 0;
        end
    end
    
    disp('Received word:')
    for bit = 0:16
        disp(sprintf('%s',num2str(recFrame(frame,1+(bit*15):15+(bit*15)))));
    end
       
    %disp('Sending the code');
    send = recFrame(frame,:);
    %decode the word with errors
    DECODED = RSdecoder(send);
    
    disp('Decoded word:')
    for bit = 0:16
        disp(sprintf('%s',num2str(DECODED(1+(bit*15):15+(bit*15)))));
    end
    
    %if it is decoded to the zero word
    if (isequal(DECODED,allEmpty))
        disp('      Succesful Decoding')
    %if it cannot be decoded and returned as it is
    elseif (isequal(DECODED,send))
        disp('      No Change')
    %if it is decoded into another codeword
    else
        disp('      Decoding Error')
    end
    
end

