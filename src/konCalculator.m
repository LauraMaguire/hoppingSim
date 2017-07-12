function kon = konCalculator(ks, cm, L, k, tether_locations)
% Calculates macroscopic kon in units of per molarity per time.

syms x n
int = exp(-0.5.*k.*(x-L/2).^2)./sum(exp(0.5.*k.*(x-tether_locations).^2));


%fun = @(x) exp(-0.5.*k.*(x-L/2).^2)./sum(exp(0.5.*k.*(x-tether_locations).^2));
kon = (ks/(cm*L))*integral(fun,0,L);
end