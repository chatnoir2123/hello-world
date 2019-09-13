function [A,rank] = my_ref_rank(A)
% Find the row echelon form of A and rank of A

[m n]=size(A);

x=1; % row
y=1; % col

% Main Process
if n>=m
    while x<m
        if A(x,y)~=0 % If there is a pivot in row x
            for l=(x+1):m
                % Make all elements below the pivot to zero
                A=ERO_3(A,l,x,-A(l,y)/A(x,y)); 
            end
            % Move to the element on the lower right side
            x=x+1;
            y=y+1;

        else
            nz=find(A(x+1:m,y))+x; % find non zero element in column y
            if isempty(nz) % No pivot in that column
                y=y+1; % Move to the element on the right side
            else % There is a pivot in column y
                A=ERO_1(A,nz(1),x); % Move up the row which has a pivot
                for l=(x+1):m
                    % Make all elements below the pivot to zero
                    A=ERO_3(A,l,x,-A(l,y)/A(x,y));
                end
                % Move to the element on the lower right side
                x=x+1;
                y=y+1;
            end
        end
    end

else % n<m, the process is the same as that of above.
    while y<=n
        if A(x,y)~=0
            for l=(x+1):m
                A=ERO_3(A,l,x,-A(l,y)/A(x,y));
            end
            x=x+1;
            y=y+1;

        else
            nz=find(A(x+1:m,y))+x; % find non zero element 
            if isempty(nz) % No pivot in that column
                y=y+1;
            else        
                A=ERO_1(A,nz(1),x);
                for l=(x+1):m
                    A=ERO_3(A,l,x,-A(l,y)/A(x,y));
                end
                x=x+1;
                y=y+1;
            end
        end
    end   
end

% Check the rank of the matrix
x=1; % row
y=1; % col
e=10^(-7); % error
rank=0;
while (x<=m)&(y<=n)
    if abs(A(x,y))<e % Take absolute value, due to computational error.
        y=y+1; % Move to the element on the right side
    else % If the row contains a pivot
        % Move to the element on the lower right side
        x=x+1;
        y=y+1;
        rank=rank+1;    
    end
end

end

