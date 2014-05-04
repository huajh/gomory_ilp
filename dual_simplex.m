function [ x,val,status,model] = dual_simplex(model)
%primal problem:
%   max         f'*x
%   subject to  A*x <= b
%               0 <= x
%
    x = [];  val = []; status = 1;
    eps = 2^-24;
    if ~all(model.c + eps >= 0 )
        disp('It''s not dual feasible');
        status = 0;
        return;
    end
    % tableau
    %______________________________
    %   basis | Abar  | bbar    
    %         | cbar' | zeta
    %
    %
    [m,n] = size(model.A);
    Abar = model.A;
    bbar = model.b;
    cbar = model.c;
    zeta0 = model.zeta;
    basis = model.basis;
    nonbasis = model.nonbasis;    
    xbar = zeros(m+n,1);         
        
    iter = 0;
    maxiter = 1000;
    while iter < maxiter
        iter = iter + 1;
       % disp(bbar);
        if all(bbar +eps>= 0 ) %optimial solution is found.
            xbar(basis) = bbar;            
            break;
        end
        % piovt selection  (i,j)
        ind = find(bbar +eps < 0);
        i = ind(1);
        if all(Abar(i,:)+eps>=0) 
            disp('LP is infeasible, DP is unbounded');
            status = 0;   
            break;        
        end       
        idx = find( Abar(i,:) + eps < 0);
        [~,idx2] = max(cbar(idx)'./Abar(i,idx));
        j = idx(idx2);        
                
        %%%%%
        % piovt (i,j)
        % update
        tmp = basis(i);
        basis(i) = nonbasis(j);
        nonbasis(j) = tmp;
        
        Atmp = [Abar,bbar;cbar',zeta0];
        %[i,j]
        Atmp2 = Atmp - repmat(Atmp(:,j)./Atmp(i,j),1,n+1).*repmat(Atmp(i,:),m+1,1);
        Atmp2(i,:) = Atmp(i,:)/Atmp(i,j);
        Atmp2(:,j) = -Atmp(:,j)/Atmp(i,j);        
        Atmp2(i,j) = 1/Atmp(i,j);
                
        
        Abar = Atmp2(1:m,1:n);
        bbar = Atmp2(1:m,n+1);
        cbar = Atmp2(m+1,1:n)';
        zeta0 = Atmp2(m+1,n+1);
        %disp(Atmp2);
    end
    
    if (status==0 ||  iter == maxiter )        
        status = 0;
        disp('no feasiable solution');
        return;     
    end
    
    x = xbar(1:model.N);
    val = zeta0;
    model.A = Abar;
    model.b = bbar;
    model.c = cbar;
    model.zeta = zeta0;
    model.basis = basis;
    model.nonbasis = nonbasis;
    status = 1;
end
