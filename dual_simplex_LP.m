function [ x,val,status] = dual_simplex_LP(f,A,b)
%primal problem:
%   min         f'*x
%   subject to  A*x <= b
%               0 <= x
%
    
    if ~all(f >= 0 )
        disp('It''s not dual feasible');
        status = 0;
        return;
    end
    % initialization
    % tableau
    %______________________________
    %   basis | Abar  | bbar    
    %         | cbar' | zeta
    %
    %
    [m,n] = size(A);
    nonbasis = 1:n;
    basis = n+1:n+m;    
    
    Abar = A;
    bbar = b;
    cbar = f;
    zeta0 = 0;
    xbar = zeros(m+n,1);     
    status = 1;
    
    iter = 0;
    while 1
        iter = iter + 1;
        if all(bbar >= 0 ) %optimial solution is found.
            xbar(basis) = bbar;        
            break;
        end
        % piovt selection  (i,j)
        ind = find(bbar +eps < 0);
        i = ind(1);
        if all(Abar(i,:)>=0) 
            disp('LP is infeasible, DP is unbounded');
            status = 0;    
            break;        
        end       
        idx = find(Abar(i,:)<0);
        [~,idx2] = max(cbar(idx)'./Abar(i,idx));
        j = idx(idx2);        

        %%%%%
        % piovt (i,j)
        % update
        tmp = basis(i);
        basis(i) = nonbasis(j);
        nonbasis(j) = tmp;
        
        Atmp = [Abar,bbar;cbar',zeta0];
        Atmp2 = Atmp - repmat(Atmp(:,j)./Atmp(i,j),1,n+1).*repmat(Atmp(i,:),m+1,1);
        Atmp2(i,:) = Atmp(i,:)/Atmp(i,j);
        Atmp2(:,j) = -Atmp(:,j)/Atmp(i,j);        
        Atmp2(i,j) = 1/Atmp(i,j);
        
        Abar = Atmp2(1:m,1:n);
        bbar = Atmp2(1:m,n+1);
        cbar = Atmp2(m+1,1:n)';
        zeta0 = Atmp2(m+1,n+1);       
    end
    
    if status==0        
        return;     
    end
    
    x = xbar(1:n);
    val = -zeta0;
    status = 1;
end

