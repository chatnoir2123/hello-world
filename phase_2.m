function [T,x_opt,opt_f,rec_dir] = phase_2(A,T,c,idx_B,B_inv)
% Conduct Phase 2
 
[m,n]=size(A); % must be in standard form
[m_T,n_T]=size(T); % T : Tableau
rec_dir=[]; % If the LP is bounded, return empty value for a recession direction.

 
% Construct idx_N
idx=1:(n_T-1);
p=length(idx_B);
idx_temp=idx;
for k=1:p
    idx_temp(idx_B(k))=0;
end
idx_N=find(idx_temp)
% find the indices of nonzero elements in idx_temp
% i.e. find the indices of variables which are not basic variables
 
% Compute 0th row of the tableau
N=A(:,idx_N);
T(1,idx_N)=(c(idx_B))'*B_inv*N-(c(idx_N))'; 
T(1,n_T)=(c(idx_B))'*T(2:m_T,n_T); % Current obj f value
reduced_cost=T(1,1:(n_T-1));

% Computational error criterion
e=10^(-6);

% Iteration
iterations=0;

% Simplex
while(~isempty(find(reduced_cost>0))) % Proceed until are reduced costs are <=0
    
    % Find the maximum reduced cost and its index
    Max=max(reduced_cost(reduced_cost>0));
    i_Max=find(reduced_cost==Max); 
    
    % Bland's rule
    if(length(i_Max)~=1)
       i_Max=min(i_Max); 
    end
    
    % When the LP is unbounded.
    ubb=(T(2:(m_T),i_Max)<=0); % When y<=0
    % e.g. T(2:(m_T),i_Max)==[-1 -1 0]'
    % then ubb=(T(2:(m_T),i_Max)<=0)==[1 1 1]'
    % sum(ubb)== 3 == (m_T-1) == (4-1)

    if(sum(ubb)==(m_T-1))
        % Current BFS
        x_B_ubb=zeros(n,1);
        x_B_ubb(idx_B)=T(2:m_T,n_T);
        x_opt=x_B_ubb; % x_opt here is not an optimal point.

        
        % Recession direction
        rec_dir=zeros(n,1);
        rec_dir(idx_B)=-T(2:m_T,i_Max);
        rec_dir(i_Max)=1;

        opt_f=-inf; % Unbounded

        
        txt = sprintf('The LP is unbounded (Phase 2). The returned x_opt is not optimal point.')
        txt = sprintf(' x = x_opt + (rec_dir) * x_%g , x_%g>=0 ',i_Max,i_Max)
        return % Terminate the function if the LP is unbounded
    end
    
    % Ratio test
    ratio_test=T(2:(m_T),(n_T))./T(2:(m_T),i_Max);
    % Find the minimum value when y>0
    y=T(2:(m+1),i_Max);
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
    for k=1:m_T
        if(k~=(i_Min+1))
            T=ERO_3(T,k,(i_Min+1),-(T(k,i_Max))/(T((i_Min+1),i_Max)));
        end
    end
    
    % Update new reduced costs, indices for basic and nonbasic variables
    reduced_cost=T(1,1:(n_T-1)); 
    idx_N(find(idx_N==i_Max))=idx_B(i_Min); 
    idx_B(i_Min)=i_Max;
    
    % Computational error adjustion
    i_prev=find(abs(T)<e);
    if ~isempty(i_prev)
       T(i_prev)=0; 
    end
    
    T
    T(1,n_T)
    iterations=iterations+1
end
 
% Optimal solution
x_opt=zeros(n,1);
x_opt(idx_B)=T(2:m_T,n_T);
 
% Optimal objective function value
opt_f=T(1,n_T);
 
end
