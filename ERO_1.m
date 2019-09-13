function [A] = ERO_1(A,i,j)
% Interchanging the ith row and the jth row

r_i=A(i,:);
A(i,:)=A(j,:);
A(j,:)=r_i;

end

