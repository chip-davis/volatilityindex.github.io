function V = GetV(options,T,spot)
%-------------------------------------------------------------------------%
% ESTIMATE OPTION-BASED SPOT VARIANCE
% input:
% options(vector [K O IV]): out-of-the-money option price (O) and Black-Scholes
% implied volatility (IV) sorted by strike price (K) in ascending order
% T(scalar): time-to-maturity in year unit
% spot(scalar): spot/forward price
% output:
% V(scalar): annualized option-implied spot variance
%-------------------------------------------------------------------------%

K = options(:,1);
O = options(:,2);
IV = options(:,3);
dK = diff(K);
F = spot;

% ATM BS implied vol
BSIV = IV((abs(F - K)==min(abs(F - K))));

% set u grid [0.1,ubar] based on ATM BSIV
ubar = sqrt((2/T)*(log(1/0.05)/BSIV^2));
uvec = (0.1:0.1:ubar)';

% L(u) -- Eqaution (3)
L_fun = @(K,dK,O,F,u) (1 - (u^2+1i*u)*(sum(dK.*exp(1i*u*(log(K(1:end-1)) ...
    - log(F))).*(O(1:end-1)./(K(1:end-1).^2)))));

L = zeros(length(uvec),1);
for ui=1:length(uvec)
    u = uvec(ui);
    L(ui,1) = L_fun(K,dK,O,F,u);
end

% uhat1 based on abs(L(u)): 1st time abs(L(u)) <= 0.2
uhat1_ind = find(abs(L)<=0.2);
if ~isempty(uhat1_ind)
    uhat1 = uvec(uhat1_ind(1));
else
    uhat1 = ubar;
end
% uhat2 -- abs(L(u)) attains minimum on [0,ubar]
uhat2 = uvec((abs(L)== min(abs(L))));
% uhat is the minimum of uhat1 and uhat2
uhat = min(uhat1,uhat2);

% V(uhat) - option-implied spot variance Eqaution (2)
V = (-2/(T*uhat^2))*real(log(L_fun(K,dK,O,F,uhat)));






