function [R] = resistor_lookup(val1, val2, Rtumor, Rtissue, Rvessel)
    rs = [Rtumor, Rtissue, Rvessel, Rvessel];
%     val2
    R = (rs(val1)+rs(val2))/2;
end