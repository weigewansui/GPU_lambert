function [a,b] = linear_coe_coeffs(t1,t2,c1,c2)
A = (t1 - t2);
a = (c1 - c2)/A;
b = -(c1*t2 - c2*t1)/A;
end
