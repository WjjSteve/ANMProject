function z = backwardGFEM(p,e,t,T,eps)

n = size(p,2);
z = zeros(n,1);
nt = size(t,2);

f = ones(n,1);
for i = 1:n
    f(i,1) = u_initial2(p(1,i),p(2,i));
end 

CFL = 0.1;
bx1 = -2*pi*p(2,:);
bx2 = 2*pi*p(1,:);   
nor = max(sqrt(bx1.^2+bx2.^2));
h = sqrt(2*pi/nt);
k = CFL*h/nor;
R = ReactMat2D(p,t);
C = ConvectMat2D(p,t,bx1,bx2);
A = StiffMat2D(p,t,ones(nt,1));
b = LoadAssembler2D(p,t,f);
LH = R/k-C/2+eps*A/2;
I = eye(length(p));
LH(e(1,:),:) = I(e(1 ,:) ,:); 
time = 1;

while time > T
    if time-k<T
        k = time-T;
        LH = R/k-C/2+eps*A/2;
        LH(e(1,:),:) = I(e(1 ,:) ,:); 
    end    
    RH = (R/k+C/2-eps*A/2)*z+b;
    RH(e(1 ,:)) = 0;
    z = LH\RH;
    time = time-k;
end

end