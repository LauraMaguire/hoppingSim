syms x
test = x - tl;
test2 = sum(x-tl);
%%
x=1;

%%
denom = sum(exp(-0.5.*k.*(x-tl).^2));

%%
int = exp(-0.5.*k.*(x-L/2).^2)./sum(-exp(0.5.*k.*(x-tl).^2));

%%
g = matlabfunction(int);

%%
test = integral(int,0,L);