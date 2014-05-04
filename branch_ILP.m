function [x,val,status,Count] = branch_ILP(f,A,b,d,k,Maxrep)
    
%
%  @author: Junhao Hua
%  @email:  huajh7@gmail.com
%  
%  create time: 2014/4/27
%  last update: 2014/5/4
%    
    % initialization
    % tableau
    %______________________________
    %   basis | Abar  | bbar    
    %         | cbar' | zeta
    %
    %
    % status
    % 0: No integer solution find
    % 1: Stop becausefind optimal solution 
    % 2: Stop because exceed max repeat times (default: 5000)   
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
    if nargin < 6
        Maxrep = 5000;
    end
    model.Maxrep = Maxrep; 
    x0 = model.b;
    val0 = model.zeta;
    bound = inf;
    Count.iter = 0;
    Count.left = 0;
    Count.right = 0;
    Count.rep = 0;
    [x,val,bound,status,Count] = bb_tree(x0,val0,model,bound,Count);
    val = - val;
end


function [xx,vval,bbound,status,Count] = bb_tree(x,val,model,bound,Count)
% Branch and bound (BB or B&B) tree
%    
    if  Count.rep > model.Maxrep
        xx = x;
        vval = val;
        bbound = bound;
        status = 2;
        return;
    end
    format rat;
    eps = 2^-24;
    Count.iter = Count.iter + 1;
    %minimize value of f'x
    [x0,val0,status0,model0] = dual_simplex(model);   
    % no optimal solution or optimal value can not be improved (\zeta^*
    % should be a integer)
    
    if status0 == 0 || -val0 > (bound-1+eps)
        xx=x; vval=val; status=status0; bbound=bound;
        Count.rep = Count.rep + 1;
        return;
    end
    
    ind = find( abs(x0-round(x0))>eps);  
    if isempty(ind) % find integer solution   
        status = 1;
        newbound = round(-val0);
        if newbound < bound   % this solution is better than the current solution hence replace            
            xx = round(x0);
            vval = val0;
            bbound = newbound;
            Count.rep = 0;    % reset
        else                  % return the input solution
            Count.rep = Count.rep + 1;
            xx=round(x);  
            vval=val;
            bbound=bound;
        end
       sf = [ 'Iter=%d: left=%d,right=%d,bound=%d,val=%.2f,x=[' repmat(' %d',1,length(xx)) ']'];
       str = sprintf(sf,Count.iter,Count.left,Count.right,bbound,vval,xx);
        disp(str);        
        return
    end
    
    privot_indx = ind(1);
    priovt_x = x0(privot_indx);
    tt = find(model0.basis == privot_indx); % always only one
    
    [m,n] = size(model0.A);
    
    % depth first search (dfs)
    
    % left child :  x_i <= floor(priovt_x) 
    Count.left = Count.left + 1;
    model1 = model0;   
    model1.A = [model0.A;zeros(1,n)];
    model1.A(end,:) = -model1.A(tt,:);
    model1.b = [model0.b;floor(priovt_x) - model0.b(tt)];
    model1.basis = [model1.basis,m+n+1];
   % disp('left');
    [x1,val1,bound1,status1,Count] = bb_tree(x0,val0,model1,bound,Count);
    Count.left = Count.left - 1;
    status=status1;    
    if status1 > 0 && bound1 < bound % if the solution was successfull and gives a better bound
       xx=x1;
       vval=val1;
       bound = bound1;
       bbound = bound1;
    else                              % not update
        xx=x0;
        vval=val0;
        bbound=bound;
    end
    if status == 2               % stop
        return;
    end
    
    %right child :  x_i >= ceil(priovt_x)
    Count.right = Count.right + 1;
    model2 = model0;   
    model2.A = [model0.A;zeros(1,n)];
    model2.A(end,:) = model2.A(tt,:);
    model2.b = [model0.b;-ceil(priovt_x) + model0.b(tt)];
    model2.basis = [model2.basis,m+n+1];
    %disp('right');
    [x2,val2,bound2,status2,Count] = bb_tree(x0,val0,model2,bound,Count);
    Count.right = Count.right - 1;       
    if status2 > 0 && bound2<bound % if the solution was successfull and gives a better bound        
        status = status2;
        xx=x2;
        vval=val2;
        bbound=bound2;
    end
    if status2 == 2
        status = 2;
    end
end

