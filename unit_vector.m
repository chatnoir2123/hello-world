function [v] = unit_vector(l,k)
% Produce a unit column vector which dimension is l, and k-th element is 1.
% All elements except k-th element are 0.

v=zeros(l,1);
v(k)=1;

end

