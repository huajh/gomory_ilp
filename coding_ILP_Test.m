
%%%%%%%
% A binary linear [n,k,d] code is a (vector) subspace C \in F^n_2 of 
% dimension k and minimum distance d.
%
%  @author: Junhao Hua
%  @email:  huajh7@gmail.com
%  
%  create time: 2014/4/27
%  last update: 2014/5/4
%
clear;clc;
for k=5
    N = 2^k-1;
    d_max = ceil((k-1)/2)*2^(k-1);
    %list_d = [13,21,29]; 
    %  4,6,7,8,12,14
    %  3,5,9,10,11,13,19,20
    %  3,4,5,6,7,8,9,10,11,12,13,14,19,20
    % list_d = [3,4,5,6,7,8,9,10,11,12,13,14,19,20];
    %list_d = [13,14];
    for d= 1:d_max-1  %i=1:length(list_d) %d= 1:d_max-1
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
       tic;
       Maxrep = 80000;
       %[x,v,status,Count] = branch_ILP(f,A,b,d,k,Maxrep);
       [x,v,status,Count] = Gomory_ILP(f,A,b,d,k);
        %options = optimoptions('intlinprog','Display','off');
        %options = optimoptions('intlinprog','MaxTime',2000);    
        %[x, v, s]  = intlinprog(f,M,A,b,[],[],lb,ub,options); 
        tot_time = toc;
        n = round(v);
        if status == 0
            stop_reason = 'No integer solution find';
        elseif status == 1
            stop_reason = 'Find optimal solution';
        elseif status == 2
            stop_reason = sprintf('Stop because exceed Max repeat times,Maxrep=%d',Maxrep);
        end
        larger = 0;
        if Gries_bound < n
            larger = n - Gries_bound;
        end
        disp(stop_reason);            
        sf = ['strictly larger=%d, Total iter=%d,time=%fs\n Gries_bound=%d, [k,d,n]=[%d,%d,%d], x=[' repmat(' %d',1,N) ']\n'];
        str = sprintf(sf,larger,Count.iter,tot_time,Gries_bound,k,d,n,x');                
        disp(str);        
    end
end




