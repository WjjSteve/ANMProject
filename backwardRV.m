function u = backwardRV(p,e,t,T)

n = size(p,2);
u = zeros(n,1);
un = u;
f = ones(n,1);
for i = 1:n   
    f(i,1) = u_initial2(p(1,i),p(2,i));
end 

CFL = 1.5;
nt = size(t,2);
bx1 = -2*pi*p(2,:);
bx2 = 2*pi*p(1,:);   
nor = max(sqrt(bx1.^2+bx2.^2));
h = sqrt(2*pi/nt);
k = CFL*h/nor;
time = 1;

if time>T
    if time-k<T
        k = time-T; 
    end    
    R = ReactMat2D(p,t);
    C = ConvectMat2D(p,t,bx1,bx2);
    b = LoadAssembler2D(p,t,f);
    LH = R/k-C/2;
    RH = (R/k+C/2)*u+b;
    I = eye(length(p));
    LH(e(1,:),:) = I(e(1 ,:) ,:); 
    RH(e(1 ,:)) = 0;
    un = u;
    u = LH\RH;
    time = time+k;
end

nt = size(t,2);
eps = zeros(nt,1);
Cvel = 0.25;
Crv = 1;

while time > T
    LH = R;
    RH = -R*(un-u)/k-C*u;
    Res = LH\RH; 
    Res = Res/max(abs(u-mean(u)));
    for K = 1:nt 
        beta_K = max(sqrt(bx1(t(1:3,K)).^2+(bx2(t(1:3,K))).^2));
        Res_K = max(abs(Res(t(1:3,K))));
        h = ComputeDiameter(p(1,t(1:3,K)),p(2,t(1:3,K)));
        eps(K) = min(Cvel*h*beta_K,Crv*h^2*Res_K);
    end
    if time-k<T
        k = time-T; 
    end     
    A = StiffMat2D(p,t,eps);
    LH = R/k-C/2+A/2;
    RH = (R/k+C/2-A/2)*u+b;
    LH(e(1,:),:) = I(e(1 ,:) ,:); 
    RH(e(1 ,:)) = 0;
    un = u;
    u = LH\RH;
    time = time-k;
end

end


