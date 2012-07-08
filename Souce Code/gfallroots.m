function roots = gfallroots(poly, n, field)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the list of zeros or roots of a polynomial               %
% in the field 'field'                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:                                                              %
%       poly:   matrix representing polynomial format in field       %
%               i.e. [1 0 2 3] for A + x + A^2x^2 + A^3x^3           %
%               where A is primitive element in the field            %
%       n:      p^m-1 in the field p^m                               %
%       field:  list of all elements in the field p^m                %
%               generated using gftuple                              %
%Output:                                                             %
%       roots:  the list of roots of the polynomial                  %
%               representing polynomail format                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    roots = [];
    %for each alpha^num in the field
    for num = 0:n-1
       %evaluate the polynomial
       value = gfpolyval(poly,num,n,field);
       %if evaluated to zero than its a root
       if(value == -Inf)
          roots = [roots num] ;
       end
    end    
end
