function [A,x_a,b,idx_B,dir] = artificial_var(A,b)
% Introduce artificial variables if necessary
% Input matrix A must be in standard form. 

[m,n]=size(A);
idx_B=zeros(1,m);
i=n;
k=0;
flag=0;
dir=0; % dir=0 : proceed to Phase 1, dir=1 : skip Phase 1 and proceed to Phase 2

% Preprocessing to find unit vectors as many as possible
while k<=(m-1) % Search unit column vectors from column 1 to m
    i=n;
    flag=0;
    k=k+1;
    v=unit_vector(m,k); 
    %a unit column vector which dimension is m, and k-th element is 1.
   
   % Find a unit column vector (or a constant multiple of it) 
   while (i>=1)&&(flag==0)
       % v=unit_vector(m,k);
      if  A(:,i)==v % If a unit column vector is found
          idx_B(k)=i;
          %k=k+1;
          flag=1; % A unit vector is found
       
      %elseif (v'*A(:,i))==(norm(A(:,i)))&&(A(:,i)~=zeros(m,1))
      elseif (v'*A(:,i))==(norm(A(:,i)))&&(norm(A(:,i))~=0)
          % e.g. If A(:,i)=[0 0 2 0]', v=[0 0 1 0]'
          % v'*A(:,i) = 2 = norm(A(:,i))
          % Avoid the case when A(:,i)=[0 0 0 0]'
          % In this case simply dividing A(;,i) into 2 is enough
          % No need to introduce an artificial variable
          b(k)=b(k)/(norm(A(:,i)));
          A=ERO_2(A,k,1/(norm(A(:,i))));
          % b(k)=b(k)/(norm(A(:,i)));
          idx_B(k)=i;
          %k=k+1;
          flag=1; % A unit vector is found
          
      else % If there is no unit column vector
       i=i-1;
      end 
   end
   %k=k+1;
end

idx_id=find(idx_B==0);
% The indices for unit column vectors that are not found 
% by the preprocessing 

if isempty(idx_id) % i.e. no artificial variable is needed
    txt = sprintf('Skip Phase 1 and proceed to Phase 2 directly')
    x_a=[];
    dir=1; % Proceed to Phase 2
else % Introduce artificial variables to take initial B: identity matrix
    for p=1:length(idx_id)
        idx_B(idx_id(p))=n+p;
        x_a(:,p)=unit_vector(m,idx_id(p));
    end
    txt = sprintf('Proceed to Phase 1')
end

end

