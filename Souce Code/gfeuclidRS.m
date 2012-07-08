function [g,r] = gfeuclidRS( a, b, t, field)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs the Euclid Algorithm on two polynomials a and b           %
% in the field until the degree of the remainder                     %
% polynomial is less than t                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:                                                              %
%       a, b:   matrix representing polynomial format in field       %
%               i.e. [1 0 2 3] for A + x + A^2x^2 + A^3x^3           %
%               where A is primitive element in the field            %
%       t:      the degree of the remainder should be <t             %
%       field:  list of all elements in the field p^m                %
%               generated using gftuple                              %
%Output:                                                             %
%       g, r:   matrix representing polynomial format in field       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %initialize the values of r0,f0,g0,r1,f1,g1 
    r0 = a;
    f0 = 0;     %represents 1;
    g0 = -Inf;  %represents 0;

    r1 = b;
    f1 = -Inf;   %represents 0;
    g1 = 0;      %represents 1;

    %divide r1 by r0
    [q2, r2] = gfdeconv(r0,r1,field);
    %update values of f2,g2
    f2 = gfsub(f0, gfconv(q2,f1,field), field);
    g2 = gfsub(g0, gfconv(q2,g1,field), field);
    
    %rearrange the values
    r0 = r1;
    r1 = r2;
    g0 = g1;
    g1 = g2;
    f0 = f1;
    f1 = f2;        
    
    %keep on dividing, updating and rearranging
    %until deg(r2)< t
    while  not( length(r1) < t+1 )%|| r1(t+1)== -Inf) 
        [q2, r2] = gfdeconv(r0,r1,field);
        f2 = gfsub(f0, gfconv(q2,f1,field), field);
        g2 = gfsub(g0, gfconv(q2,g1,field), field);
       
        r0 = r1;
        r1 = r2;
        g0 = g1;
        g1 = g2;
        f0 = f1;
        f1 = f2;        
    end

    g = g2;
    r = r2;
end
