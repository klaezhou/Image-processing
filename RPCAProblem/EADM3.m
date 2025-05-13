function [M,Y,X] = EADM3(Psize,Omega,b,opts)

tau   = opts.tau;   %%% parameter in model
mu    = opts.mu;
beta  = opts.beta;  %%% penalty of ADMM
maxit = opts.MaxIt; 

%% 初始值
X  = zeros(Psize);
Y  = zeros(Psize);
Z  = zeros(Psize);
Lbd= zeros(Psize);
M     = zeros(Psize);     M(Omega) = b;
rankr = zeros(1,maxit);
SP    = zeros(1,maxit);

for iter  = 1:maxit    
    %%%% Step 1. X - subproblem (low-rank)
    T = M + Lbd/beta - Y - Z;
    [U,D,V] = svd(T,'econ');
    D = diag(D);
    ind = find(D > 1/beta);
    D = diag(D(ind) - 1/beta);
    X = U(:,ind) * D * V(:,ind)';
        
    %%%%%% Step 2. Y - subproblem (sparse)
    T = M + Lbd/beta - X - Z;
    Y = sign(T).*max(abs(T)-tau/beta, 0);    
    
    %%%%% Step 3. Z - subprolbme 
    Z = M + Lbd/beta - X - Y;
    sfactor = mu*beta/(1+mu*beta);    
    Z(Omega)= Z(Omega)*sfactor;   
   
    %%%%% update Lbd
    Lbd = Lbd - beta * (X + Y + Z - M); 
    
    %%%%%%%%%%%%%%%%%%%% 记录结果
    rankr(iter)= length(ind);
    SP(iter)   = sum(abs(Y(:))>1);
    fprintf('It:%2d, rank:%3d, spa:%7d,\n',iter,rankr(iter),SP(iter));   
end 


