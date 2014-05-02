function [x,val,status] = Gomory_ILP(f,A,b,d,k)
    
    
    % initialization
    % tableau
    %______________________________
    %   basis | Abar  | bbar    
    %         | cbar' | zeta
    %  
    
    % Use exact rational arithmetic 
    format rat;
    eps = 2^-24;
    
%     [m,n] = size(A);
%     model.A = A;
%     model.b = b;
%     model.c = f;
%     model.zeta = 0;
%     model.nonbasis = 1:n;
%     model.basis = n+1:n+m;
  
    %initialization: the linear program solution
    [N,N] = size(A);
    model.A = 1/(2^(k-1))*(ones(N,N) - 2*(-A') );
    model.b = d/2^(k-1)*ones(N,1);
    model.c = 1/2^(k-1)*ones(N,1);
    model.zeta = -(2^k - 1)*d/(2^(k-1));
    model.basis = 1:N;
    model.nonbasis = N+1:N+N;
    x0 = model.b;
    status0 = 1;
    val0 = model.zeta;
    
    p = N;
    old_x = [];
    old_v = [];
    iter = 0;    
    %maxIter = 1000;

    while 1
        iter = iter + 1;                       
        sf = [ '%d. x0 = [' repmat(' %.2f',1,p) '] val0=%.2f '];
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
        %[~,t] = max(abs(x0(ind)-round(x0(ind))));
        t = 1;
        idx = find(model.basis == ind(t));
        % pivot selection cutting plane / in which row ?
        % add a new constraint
        Abar = model.A;
        [m,n] = size(Abar);
        Anew = - ( Abar(idx,:) - floor(Abar(idx,:)+eps) );
        f0= - ( model.b(idx) - floor(model.b(idx)+eps) );
        model.A = [model.A;Anew];
        model.b = [model.b;f0];
        model.basis = [model.basis,m+n+1];
        
        [x0,val0,status0,model] = dual_simplex(model,p);   %%remove: get the result of linear programming when first running.
    end
    if status0 <=0 %|| iter > maxIter
        disp('no feasibale solution');
        x = old_x; 
        val = old_v; 
        status = 0;
        return;
    end
end

function [ x,val,status,model] = dual_simplex(model,p)
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
    
    x = xbar(1:p);
    val = -zeta0;
    model.A = Abar;
    model.b = bbar;
    model.c = cbar;
    model.zeta = zeta0;
    model.basis = basis;
    model.nonbasis = nonbasis;
    status = 1;
end

% %     iter = 0;
% %     x0 = [zeros(n,1);b];
% %     A =  [A,eye(m)];
% %     B = [n+1:n+m]';
% %     c = [-f;zeros(n,1)];
% %     p = n;
% %     n = n + m;
%     while 0
%         iter = iter + 1;
%         [val0,x0,B,status0] = dual_simplex_bb(c,A,b,eps,x0,B);        
%         if status0 == 2 % no feasiable solution
%             x = old_x; 
%             val = old_v; 
%             status = status0;
%             break;
%         end
%         if( x0-floor(x0+eps) < eps)
%             x = x0(1:p);
%             val = val0;
%             status = 1;
%             break;
%         end
%         N = setdiff([1:n],B);
%         D = A(:,B) \ A(:,N);
%         %set br_var = most fractional
%         [~,br_var] = max(abs(x0-round(x0)));
%         br_value = floor(x0(br_var)+eps);
%         %get z= the gomery cutting plane and add the cut to the system
%         z = D(br_var,:);
%         A = [A zeros(m,1);zeros(1,n+1)];
%         m = m+1;
%         n = n+1;
%         
%         A(m,br_var) = 1;
%         A(m,N) = floor(z);
%         b = [b;br_value];
%         c = [c;0];
%         x0 = [x0;b(m)-sum(A(m,1:n-1))*x0];
%         B = [B;n];
%     end
