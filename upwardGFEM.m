function [up,u,k] = upwardGFEM(p,e,t,T,eps)

n = size(p,2);
u = zeros(n,1);
for i = 1:n
    u(i,1) = u_initial(p(1,i),p(2,i));
end
up = u;

CFL = 0.1;
nt = size(t,2);
bx1 = -2*pi*p(2,:);
bx2 = 2*pi*p(1,:);   
nor = max(sqrt(bx1.^2+bx2.^2));
h = sqrt(2*pi/nt);
k = CFL*h/nor;
R = ReactMat2D(p,t);
C = ConvectMat2D(p,t,bx1,bx2);
A = StiffMat2D(p,t,ones(nt,1));
LH = R/k+C/2+eps*A/2;
I = eye(length(p));
LH(e(1,:),:) = I(e(1 ,:) ,:); 
time = 0;

while time < T
    if time+k>T
        k = T-time;
        LH = R/k+C/2+eps*A/2;
        LH(e(1,:),:) = I(e(1 ,:) ,:); 
    end    
    RH = (R/k-C/2-eps*A/2)*u;
    RH(e(1 ,:)) = 0;
    up = u;
    u = LH\RH;
    time = time+k;
end

end