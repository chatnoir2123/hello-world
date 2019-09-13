function [T,x_opt,opt_f,rec_dir] = my_simplex(A,eq,b,c)
% Convert constraints into standard form
[A_stan,b_stan,c_stan]=standard_form(A,eq,b,c)
 
% Compute the rank 
[A_ref,rank] = my_ref_rank(A_stan)
[m_stan,n_stan]=size(A_stan);
if rank==m_stan
    txt = sprintf('rank A=m, might be possible to skip phase 1')
else
    txt = sprintf('Phase 1 must be proceeded')
end
 
% Introduce artificial variables if necessary
[A,x_a,b,idx_B,dir] = artificial_var(A_stan,b_stan);
dir
if dir==1 % No artificial variable, skip Phase 1
    [m,n]=size(A);
    T=zeros((m+1),(n+1));
    T(2:(m+1),:)=[A b];
    [T,x_opt,opt_f,rec_dir]=phase_2(A,T,c_stan,idx_B,eye(m)); % Initial B_inv = Identity matrix
elseif dir==0 % Artificial variables are added, proceed to Phase 1
    [A,T_2,T,idx_B,idx_N,B_inv] = phase_1(A,x_a,b,idx_B)
    [T,x_opt,opt_f,rec_dir]=phase_2(A,T_2,c_stan,idx_B,B_inv);
end
end
