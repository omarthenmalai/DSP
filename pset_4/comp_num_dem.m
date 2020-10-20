function [b,a] = comp_num_dem(b1,a1,b2,a2)
    b3 = conv(b1, a2); 
    b4 = conv(b2, a1);
    a = conv(a1, a2).*2; % common denominator
    b = (b3 - b4);
end