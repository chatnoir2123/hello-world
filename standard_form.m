function [A_stan,b_stan,c_stan] = standard_form(A,eq,b,c)
% Convert constraints to standard form with b>=0
% Details regarding the column vector eq
% <= : -1
% = : 0
% >= : 1

[m n]=size(A);

k=1;
% Make b>=0
while k<=m
   if b(k)<0 % b must be a non-negative vector
       A=ERO_2(A,k,-1);
       eq(k)=(-1)*eq(k);
       b(k)=(-1)*b(k);
   end
   k=k+1;
end

k=1;
p=1;

% Introduce slack variables and extend A, c if necessary
while k<=m
    if eq(k)==-1 % <=
        A(:,n+p)=unit_vector(m,k);
        c(n+p)=0;
        p=p+1;
    elseif eq(k)==1 % >=
        A(:,n+p)=-unit_vector(m,k);
        c(n+p)=0;
        
        % If b(k)=0, multiply -1 to row k
        % By doing this we can avoid unnecessary complexity 
        % (i.e. introducing artificial variables)
        if b(k)==0
            A=ERO_2(A,k,-1);
        end
        p=p+1;
    end
    k=k+1;
end

A_stan=A;
b_stan=b;
c_stan=c;

end

