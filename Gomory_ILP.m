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
    model.N=N;
    x0 = model.b;
    status0 = 1;
    val0 = model.zeta;
        
    old_x = [];
    old_v = [];
    iter = 0;    
    %maxIter = 1000;

    while 1
        iter = iter + 1;                       
        sf = [ '%d. x0 = [' repmat(' %.2f',1,N) '] val0=%.2f '];
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
        
        %  cutting plane / in which row ?
        % add a new constraint
        Abar = model.A;
        [m,n] = size(Abar);
        Anew = - ( Abar(idx,:) - floor(Abar(idx,:)+eps) );
        f0= - ( model.b(idx) - floor(model.b(idx)+eps) );
        model.A = [model.A;Anew];
        model.b = [model.b;f0];
        model.basis = [model.basis,m+n+1];
        
        [x0,val0,status0,model] = dual_simplex(model);   %%remove: get the result of linear programming when first running.
        val0 = -val0;
    end
    if status0 <=0 %|| iter > maxIter
        disp('no feasibale solution');
        x = old_x; 
        val = old_v; 
        status = 0;
        return;
    end
end

