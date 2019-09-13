function [A] = ERO_3(A,i,j,c)
% Add a constant multiple of one row to another
% Mainâ�� a=~ �� �ϸ� a���� ��ȭ ����.

r_i=A(i,:);
r_j=A(j,:);

r_i_new=r_i+c*r_j;
A(i,:)=r_i_new;

end

