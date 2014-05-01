
%%%%%%%
% A binary linear [n,k,d] code is a (vector) subspace C \in F^n_2 of 
% dimension k and minimum distance d.
%
clear;clc;

for k=5
    N = 2^k-1;
    d_max = ceil((k-1)/2)*2^(k-1);
    %list_d = [3,4,5,6]; 
    %list_d = 3;%[3,5,6,7,9,11,13,19,20];
    for d=3%1:d_max-1%length(list_d)
        %d = list_d(i);
        % n = ?
        %Griesmer bound 
        Gries_bound = 0;
        for t=0:k-1
            Gries_bound = Gries_bound + ceil(d/2^t);
        end
        
        % solve integer linear programming by Gomory's algorithm       
        f = ones(N,1);
        A = -bin_mat_A(k);
        b = -d*ones(N,1);
        lb = zeros(N,1);
        ub = inf*ones(N,1);
        M = 1:N;
        yidx = true(N,1);
        e=2^-24;
        itermax = 10^5;
        [x v s]= IP1(f,A,b,[],[],lb,ub,M,e);
       %[x,v,s] = Gomory_ILP(f,A,b);
        %options = optimoptions('intlinprog','Display','off');
        %options = optimoptions('intlinprog','MaxTime',2000);    
        %[x, v, s]  = intlinprog(f,M,A,b,[],[],lb,ub,options); 
        n = round(v);
       % if v > Gries_bound
            str0 = sprintf('Gries_bound = %d, ',Gries_bound);
            str = sprintf('[k, d, n] = [%d, %d, %d], ',k,d,n);
            sf = [ 'x = [' repmat(' %d',1,N) ']\n\n'];
            str2 = sprintf(sf,x');
            disp([str0,str, str2]);
       % end
        

    end
end




