function [y, dydx, dydA, d_dydx_dA] =CSTairfoil(A,x)


N1 = 0.5;
N2 = 1;

C = ((x.^N1)).*(1-x).^N2;

%% create Bernstein polynomial

n = length(A);

for v = 0:n-1
    Sx(v+1,:) = nchoosek(n-1,v)*x.^v.*(1-x).^(n-1-v);
    dSdx(v+1,:) = x.^v.*nchoosek(n - 1, v).*(1 - x).^(n - v - 2).*(v - n + 1) + v*x.^(v - 1).*nchoosek(n - 1, v).*(1 - x).^(n - v - 1);

end

%%
yb = zeros(1,length(x));

for i = 1:n
    yb(1,:) = yb(1,:) + A(i).*Sx(i,:);
    dydA(i,:) = C.*Sx(i,:);
end

y = C.*yb;

%%

dCdx = N1*x.^(N1 - 1).*(1 - x).^N2 - N2*x.^N1.*(1 - x).^(N2 - 1);

dybdx = zeros(1,length(x));

for i=1:n
    dybdx = dybdx + A(i)*dSdx(i,:);
    d_dydx_dA(i,:) = dCdx.*Sx(i,:) + C.*dSdx(i,:);
end

dydx = dCdx.*yb + C.*dybdx;

dydx(x==0) = inf;
dydx(x==1) = dydx(end-1);

d_dydx_dA(:,end) = d_dydx_dA(:,end-1);
