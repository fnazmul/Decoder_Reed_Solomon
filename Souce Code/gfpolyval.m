function val = gfpolyval(poly,point,n,field)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate a polynomial at given point in the field 'field'          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:                                                              %
%       poly:   matrix representing polynomial format in field       %
%               i.e. [1 0 2 3] for A + x + A^2x^2 + A^3x^3           %
%               where A is primitive element in the field            %
%       point:  i means alpha^i                                      %
%       n:      p^m-1 in the field p^m                               %
%       field:  list of all elements in the field p^m                %
%               generated using gftuple                              %
%Output:                                                             %
%       val:    the result of evaluation as the power of alpha       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    len = length(poly);
    val = -Inf;
    for position = 0:len-1
        %(alpha^point)^postion
        xPow = mod(point*position,n);
        %r(position+1).(alpha^point)^postion
        xPos = gfmul(poly(position+1),xPow,field);
        %summation of all products
        val = gfadd(val,xPos,field);
    end
end