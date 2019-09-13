function [A] = ERO_2(A,i,c)
% Multiply ith row by nonzero constant c

if c==0
    error('must be nonzero constant');
end
A(i,:)=c*A(i,:);

end

