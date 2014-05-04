
%%%%%%%
% A binary linear [n,k,d] code is a (vector) subspace C \in F^n_2 of 
% dimension k and minimum distance d.
%
clear;clc;
%addpath('.\mip');
%addpath('.\Kartik');
for k=6
    N = 2^k-1;
    d_max = ceil((k-1)/2)*2^(k-1);
    %list_d = [13,21,29]; 
    %list_d = 3;%[3,5,9,11,13,19,20];
    for d= 1:d_max-1 %1:length(list_d) 
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
       % [x v s]= IP1(f,A,b,[],[],lb,ub,M,e);
       %[x,v] = miprog(f,A,b,[],[],lb,ub,yidx);
       %[c,A,b,v,x,B,iter] = CuttingPlane(f,A,b,eps,M,itermax)
       Maxrep = 30000;
       [x,v,status,Count] = branch_ILP(f,A,b,d,k,Maxrep);
       %[x,v,s] = Gomory_ILP(f,A,b,d,k);
        %options = optimoptions('intlinprog','Display','off');
        options = optimoptions('intlinprog','MaxTime',2000);    
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
            disp(stop_reason);
            sf = ['Total iter=%d, Gries_bound=%d, [k,d,n]=[%d,%d,%d], x=[' repmat(' %d',1,N) '],time=%fs\n\n'];
            str = sprintf(sf,Count.iter,Gries_bound,k,d,n,x',tot_time);
            disp(str);
       % end
        

    end
end




