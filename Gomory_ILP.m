function [x,val,status] = Gomory_ILP(f,A,b)
    
    
    % initialization
    % tableau
    %______________________________
    %   basis | Abar  | bbar    
    %         | cbar' | zeta
    %  
    
    % Use exact rational arithmetic 
    format rat;
    eps = 2^-24;
    
    [m,n] = size(A);
    model.A = A;
    model.b = b;
    model.c = f;
    model.zeta = 0;
    model.nonbasis = 1:n;
    model.basis = n+1:n+m;
    
    old_x = [];
    old_v = [];
    iter = 0;    
    %maxIter = 1000;
    while 1
        iter = iter + 1;        
        [x0,val0,status0,model] = dual_simplex(model);
        sf = [ '%d. x0 = [' repmat(' %.2f',1,n) '] val0=%.2f '];
        str = sprintf(sf,iter,x0,val0);
        disp(str);
        if status0<=0 % no feasiable solution
            x = old_x; 
            val = old_v; 
            status = status0;
            break;
        end        
        ind = find( abs(x0 - floor(x0+eps)) > eps);
        if isempty(ind)  % find feasiable solution
            x = round(x0);
            val = val0;
            status = 1;
            break;
        end
        [~,t] = max(abs(x0(ind)-round(x0(ind))));
        %t = 1;
        idx = find(model.basis == ind(t));
        % pivot selection / in which row ?
        % add a new constraint
        Abar = model.A;
        [m,n] = size(Abar);
        Anew = - ( Abar(idx,:) - floor(Abar(idx,:)+eps) );
        f0= - ( model.b(idx) - floor(model.b(idx)+eps) );
        model.A = [model.A;Anew];
        model.b = [model.b;f0];
        model.basis = [model.basis,m+n+1];
                
    end
    if status0 <=0 %|| iter > maxIter
        disp('no feasibale solution');
        x = old_x; 
        val = old_v; 
        status = 0;
        return;
    end
end

function [ x,val,status,model] = dual_simplex(model)
%primal problem:
%   min         f'*x
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
    %maxiter = 100;
    while 1
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
    
    if (status==0)
        disp('no feasiable solution');
        return;     
    end
    
    x = xbar(1:n);
    val = -zeta0;
    model.A = Abar;
    model.b = bbar;
    model.c = cbar;
    model.zeta = zeta0;
    model.basis = basis;
    model.nonbasis = nonbasis;
    status = 1;
end
