clear;
clc;
tic;
% set parameters
alpha = 1/3;
beta = 0.99;
sigma = 2;
a_low = 0;
rho = 0.5;
delta = 0.025;
sigma_std = 0.2;



% several steps/parts to solution



% 3. Discretize the exogenous state variable using R.. method
m = 5; % a grid with m points
[zgrid, p] = rouwenhorst(rho,sigma_std,m);
zgrid = exp(zgrid);
p_all = p;
% then find invariant p
for i = 1 : 10000
    p_all = p_all * p_all;
end

PI = p_all(1,:);
% then aggregate effective labor supply
N_s = zgrid * PI';

% 4: Discretize the endougenous state variable
a_low = 0;
a_hi_guess = 5;
num_a = 500;
a = linspace(a_low,a_hi_guess,num_a);

% 5:
% a) guess an equilibrium value for agg K, associated with the factor
% prices
err = 1;

K_max = 6;
K_min = 0.5;
K_guess = (K_max+K_min)/2;
step = 0;
data = [];
while err>0.0001

 K_guess = (K_min+K_max)/2;
%K_guess = 4.95;
w = (1-alpha)*(K_guess/N_s)^alpha;
r = alpha * (N_s/K_guess)^(1-alpha) + 1 - delta;



% policy function:

% CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, r*a',a);
    cons = bsxfun(@plus, cons, permute(zgrid*w, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons < 0) = -inf;
    
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(5, num_a);
    
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    iteration = 0;
    while v_tol >.0001
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        v_mat = ret + beta * repmat(permute(p * v_guess,[3 2 1]), [num_a, 1 ,1]);
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        
        [vfn, pol_indx] = max(v_mat, [], 2);
        vfn = shiftdim(vfn,2);
        % The distance between value function and v_guess
        v_tol = max(abs(vfn - v_guess));
        v_guess = vfn;
        iteration = iteration +1;
  
    end
    
    % KEEP DECSISION RULE
    %pol_fn = a(pol_indx);
    pol_indx = permute(pol_indx, [3 1 2]);
    % policy function
    pol_fn = a(pol_indx);
    
    
    % Find steady state distribution
    % set up initial distribution
    Mu = ones(5,num_a)/(5 * num_a);
    
    
         dis= 1;
      iter = 0;
     while dis > 0.00001
     %for i = 1 :50
    % ITERATE OVER DISTRIBUTIONS
    [state, a_num, mass] = find(Mu); % find non-zero indices
    
    iter = iter +1;
    MuNew = zeros(size(Mu));
    for ii = 1:length(state)
        apr_num = pol_indx(state(ii), a_num(ii)); % which a_num to go next period by policy fn?
        MuNew(:, apr_num) = MuNew(:, apr_num) + (p(state(ii), :) * mass(ii))' ;% which mass of households goes to which exogenous state?
        dis = max(abs(Mu(:) - MuNew(:)));
        
    end
    Mu = MuNew;
    end
    % Find stationary distribution
    Mu;
    
    % Check total capital 'market'
    
    K = Mu.*pol_fn;
    agg_K = sum(K(:));
    if agg_K>K_guess
        K_min = (K_min+K_guess)/2;
       
    else
        K_max = (K_max+K_guess)/2;
       
    end
    step = step + 1;
    KK=[ step,K_guess,agg_K]'
    data = [data KK];

    err = abs(K_guess-agg_K)
    
end
    
agg_K;
toc;
figure(1)
plot(data(1,:),data(2,:),'r')
hold on
plot(data(1,:),data(3,:),'bl')
xlabel('step')
title('convergence process')
legend('K-guess','agg-K','location','southeast')

figure(2)
plot(pol_fn(1,:),'g')
hold on
plot(pol_fn(2,:),'b')
hold on
plot(pol_fn(3,:),'bl')
hold on
plot(pol_fn(4,:),'y')
hold on
plot(pol_fn(5,:),'p')



% Wealth and Gini index
wealth_dist = Mu.*a;
plot(wealth_dist(1,:))

