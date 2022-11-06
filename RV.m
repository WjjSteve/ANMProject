% This is the code that implement the RV
close all;
clear;

g = Rectg(-2,-2.5,2,1.5);
hmax = 1/8;
[p,e,t] = initmesh(g,'hmax',hmax);
% get the trangular elements of the domain

n = size(p,2);
u = zeros(n,1);
for i = 1:n
    u(i,1) = u_initial(p(1,i),p(2,i));
end

CFL = 0.5;
T = 1;

R = ReactMat2D(p,t);

bx1 = cos(u);
bx2 = -sin(u);
% vector field f'(u)
nor = max(sqrt(bx1.^2+bx2.^2));
k = CFL*hmax;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = ConvectMat2D(p,t,bx1,bx2);
LH = R;
RH = -C;
A_r = LH\RH;
w1 = A_r*u;
bx1 = cos(u+k/2*w1);  bx2 = -sin(u+k/2*w1);
C = ConvectMat2D(p,t,bx1,bx2);
RH = -C;
A_r = LH\RH;
w2 = A_r*(u+k/2*w1);
bx1 = cos(u+k/2*w2);  bx2 = -sin(u+k/2*w2);
C = ConvectMat2D(p,t,bx1,bx2);
RH = -C;
A_r = LH\RH;
w3 = A_r*(u+k/2*w2); 
bx1 = cos(u+k*w3);  bx2 = -sin(u+k*w3);
C = ConvectMat2D(p,t,bx1,bx2);
RH = -C;
A_r = LH\RH;
w4 = A_r*(u+k*w3);
w = w1+2*w2+2*w3+w4;
w(e(1,:)) = 0;
up = u;
u = u + k/6*w;
   
time = k;
k_l = k;

nt = size(t,2);
eps = zeros(nt,1);
Cvel = 0.25;
Crv = 1.0;
Res = zeros(n,1);

while time<T
    bx1 = cos(u);
    bx2 = -sin(u);
    % vector field f'(u)
    nor = max(sqrt(bx1.^2+bx2.^2));
    k = CFL*hmax/nor;
    % Time steps
    if time+k>T
        k=T-time;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LH = R;
    ResidualMatrix;
    RH = -(C+A);
    A_r = LH\RH;
    w1 = A_r*u;
    
    bx1 = cos(u+k/2*w1);  bx2 = -sin(u+k/2*w1);
    C = ConvectMat2D(p,t,bx1,bx2);
    RH = -(C+A);
    A_r = LH\RH;
    w2 = A_r*(u+k/2*w1);

    bx1 = cos(u+k/2*w2);  bx2 = -sin(u+k/2*w2);
    C = ConvectMat2D(p,t,bx1,bx2);
    RH = -(C+A);
    A_r = LH\RH;
    w3 = A_r*(u+k/2*w2); 

    bx1 = cos(u+k*w3);  bx2 = -sin(u+k*w3);
    C = ConvectMat2D(p,t,bx1,bx2);
    RH = -(C+A);
    A_r = LH\RH;
    w4 = A_r*(u+k*w3);

    w = w1+2*w2+2*w3+w4;
    w(e(1,:)) = 0;
    up = u;
    u = u + k/6*w;
    k_l = k;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% RK4 %%%%%%%%%%%%%%
    time = time+k;
end

pdeplot(p,e,t,"XYData",u); 
pbaspect([1 1 1]);
caxis([0.5 1.3]);
colormap turbo;
saveas(gcf,"T22RVh8.png");