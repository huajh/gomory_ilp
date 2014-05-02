function [x,val,status] = branch_ILP(f,A,b,d,k)
    
    
    % initialization
    % tableau
    %______________________________
    %   basis | Abar  | bbar    
    %         | cbar' | zeta
    %  
    
    % Use exact rational arithmetic 
    format rat;
    eps = 2^-24;
    
    [N,N] = size(A);
    model.A = 1/(2^(k-1))*(ones(N,N) - 2*(-A') );
    model.b = d/2^(k-1)*ones(N,1);
    model.c = 1/2^(k-1)*ones(N,1);
    model.zeta = -(2^k - 1)*d/(2^(k-1));
    model.basis = 1:N;
    model.nonbasis = N+1:N+N;
    model.N = N;
    x0 = model.b;
    val0 = model.zeta;
    bound = inf;
    [x,val,bound,status] = recursion_tree(x0,val0,model,bound);
    val = - val;
end

function [xx,vval,bbound,status] = recursion_tree(x,val,model,bound)
    
    format rat;
    eps = 2^-24;
    [x0,val0,status0,model0] = dual_simplex(model);
    if status0<=0 || -val0 > bound  
        xx=x; vval=val; status=status0; bbound=bound;
        return;
    end
    
    ind = find( abs(x0-round(x0))>eps);  
    if isempty(ind) % find integer solution   
        status=1;        
        if -val0 < bound    % this solution is better than the current solution hence replace            
            xx = round(x0);
            vval=val0;
            bbound=-val0;
        else
            xx=x;  % return the input solution
            vval=val;
            bbound=bound;
        end
        sf = [ ' val0=%.2f  x0 = [' repmat(' %.2f',1,length(x0)) ']'];
        str = sprintf(sf,val,xx);
        disp(str);        
        return
    end
    
    privot_indx = ind(1);
    priovt_x = x0(privot_indx);
    tt = find(model0.basis == privot_indx); % always only one
    
    [m,n] = size(model0.A);
    
    % depth first search (dfs)
    
    % left child :  x_i <= floor(priovt_x) 
    model1 = model0;   
    model1.A = [model0.A;zeros(1,n)];
    model1.A(end,:) = -model1.A(tt,:);
    model1.b = [model0.b;floor(priovt_x) - model0.b(tt)];
    model1.basis = [model1.basis,m+n+1];
    
    [x1,val1,bound1,status1] = recursion_tree(x0,val0,model1,bound);
    status=status1;
    if status1 > 0 && bound1 < bound % if the solution was successfull and gives a better bound
       xx=x1;
       vval=val1;
       bound = bound1;
       bbound = bound1;
    else
        xx=x0;
        vval=val0;
        bbound=bound;
    end
    
    
    %right child :  x_i >= ceil(priovt_x)
    model2 = model0;   
    model2.A = [model0.A;zeros(1,n)];
    model2.A(end,:) = model2.A(tt,:);
    model2.b = [model0.b;-ceil(priovt_x) + model0.b(tt)];
    model2.basis = [model2.basis,m+n+1];
    
    [x2,val2,bound2,status2] = recursion_tree(x0,val0,model2,bound);
    
    if status2 >0 && bound2<bound % if the solution was successfull and gives a better bound
        status = status2;
        xx=x2;
        vval=val2;
        bbound=bound2;
    end 
    
end

