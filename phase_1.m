function [A,T_2,T,idx_B,idx_N,B_inv] = phase_1(A,x_a,b,idx_B)
% Conduct Phase 1

[m,n]=size(A);
[m_a,n_a]=size(x_a);
idx_B_org=idx_B; % Define it in order to get B_inv

e=10^(-9);  % Phase 1 test criterion

idx_art=(n+1):1:(n+n_a); % Indices of artificial variables

% Construct idx_N
idx=1:(n+n_a); % The entire endices including those of artificial variables

p=length(idx_B);
idx_temp=idx;
for k=1:p
    idx_temp(idx_B(k))=0;
end
idx_N=find(idx_temp); 
% find the indices of nonzero elements in idx_temp
% i.e. find the indices of variables which are not basic variables


% Objective function for Phase 1
c=zeros((n+n_a),1);
for k=(n+1):(n+n_a)
    c(k)=1; 
end

% Construct simplex tableau
T=zeros((m+1),(n+n_a+1));
A_aug=[A x_a];
N=A_aug(:,idx_N);
T(2:(m+1),1:n)=A;
T(2:(m+1),(n+1):(n+n_a))=x_a;
T(2:(m+1),(n+n_a+1))=b;
T(1,idx_N)=(c(idx_B))'*N-(c(idx_N))'; % Reduced costs
T(1,(n+n_a+1))=(c(idx_B))'*b; % Current obj f value

reduced_cost=T(1,1:(n+n_a));

% Simplex
while(~isempty(find(reduced_cost>0))) % Proceed until are reduced costs are <=0
    % Find the maximum reduced cost and its index
    Max=max(reduced_cost(reduced_cost>0));
    i_Max=find(reduced_cost==Max);
    
    % Bland's rule
    if(length(i_Max)~=1)
       i_Max=min(i_Max);  
    end
    
    % Ratio test
    ratio_test=T(2:(m+1),(n+n_a+1))./T(2:(m+1),i_Max);
    y=T(2:(m+1),i_Max);
    % Find the minimum value when y>0
    idx_y=find(y>0);
    [Min,i_Min_temp]=min(ratio_test(idx_y));

    i_Min=idx_y(i_Min_temp); % Indicies of basic variables corresponding to minimum ratio
    
    
    % Bland's rule
    if(length(i_Min)~=1)
       i_B_Min=min(idx_B(i_Min));
       i_Min=find(idx_B==i_B_Min)
    end
    
    T=ERO_2(T,(i_Min+1),1/(T((i_Min+1),i_Max))); % Make the pivot to 1
    % Make all elements in that column to 0 except the pivot 
    for k=1:(m+1)
        if(k~=(i_Min+1))
            T=ERO_3(T,k,(i_Min+1),-(T(k,i_Max))/(T((i_Min+1),i_Max)));
        end
    end
    
    % Update new reduced costs, indices for basic and nonbasic variables
    reduced_cost=T(1,1:(n+n_a));
    idx_N(find(idx_N==i_Max))=idx_B(i_Min);
    idx_B(i_Min)=i_Max;
end

% Check whether it is good to proceed to Phase 2 or not
obj=T(1,(n+n_a+1))
if abs(obj)<e
    % If there are no redundant constraint, save B_inv again after
    % modifying T
    
    % Check whether there are artificial variables in the basis or not
    l_art=length(idx_art);
    for k=1:l_art
        B_art=find(idx_B==idx_art(k));
        if(~isempty(B_art))
             % if -> redundant constraint, eliminate it 
             % else -> change that artificial variable to legitimate variable
            if(norm(T((B_art+1),1:n))<10^(-6)) % All elements corresponding to legitimate variables are zeros (computational error)
                txt = sprintf('row %g is redundant. Eliminate this row.',(B_art+1))
                T((B_art+1),:)=[];
                m=m-1;
                idx_B_org(idx_B_org==idx_art(k))=[];
                idx_B(idx_B==idx_art(k))=[];
                A(B_art,:)=[];
            else % At least one of the elements is not zero.
               idx_legit_art=find(abs(T((B_art+1),1:n))>10^(-6)); % Index of legitimate var on art row which is not 0
               if(length(idx_legit_art)~=1) % Bland's Rule
                   idx_legit_art=min(idx_legit_art);    
               end
               
               % Elementary Row Opertations
               % Pivot element : (B_art+1), (idx_legit_art)
               T=ERO_2(T,(B_art+1),1/(T((B_art+1),idx_legit_art)));
               for k=1:(m+1)
                    if(k~=(B_art+1))
                        T=ERO_3(T,k,(B_art+1),-(T(k,idx_legit_art))/(T((B_art+1),idx_legit_art)));
                    end
               end
               
               % Change indices
               idx_N(find(idx_N==idx_legit_art))=idx_B(B_art);
               idx_B(B_art)=idx_legit_art;
            end
        end
    end
      
    % Construct a tableau for Phase 2
    T_2=T(:,1:n); % Phase 2 tableau
    T_2(:,(n+1))=T(:,(n+n_a+1)); % Shift current x_B
    T_2(1,:)=zeros(1,(n+1)); % initialize all reduced costs
    B_inv=T(2:(m+1),idx_B_org);
    
else % If objective function value is not zero
    error('The LP is infeasible. (Phase 1)')
end

end

