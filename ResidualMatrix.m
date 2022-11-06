C = ConvectMat2D(p,t,bx1,bx2);
LH = R;
RH = R*(u-up)/k_l+C*u;
Res = LH\RH; 
Res = Res/max(abs(u-mean(u)));
for K = 1:nt 
    beta_K = max(sqrt(bx1(t(1:3,K)).^2+(bx2(t(1:3,K))).^2));
    Res_K = max(abs(Res(t(1:3,K))));
    h = ComputeDiameter(p(1,t(1:3,K)),p(2,t(1:3,K)));
    eps(K) = min(Cvel*h*beta_K,Crv*h^2*Res_K);
end
A = StiffMat2D(p,t,eps);