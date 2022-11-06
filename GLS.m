% This is the code that implement the GLS
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
u0 = u;

CFL = 0.5;
T = 1;

R = ReactMat2D(p,t);
time = 0;

while time<T
    bx1 = cos(u);
    bx2 = -sin(u);
    % vector field f'(u)
    nor = max(sqrt(bx1.^2+bx2.^2));
    k = CFL*hmax;
    tau = ((2/k)^2+(2*nor/hmax)^2)^(-1/2);
    % Time steps
    if time+k>T
        k=T-time;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C = ConvectMat2D(p,t,bx1,bx2);
    S = SDAssembler2D(p,t,bx1,bx2); 
    LH = R+tau*C';
    RH = -(C+tau*S);
    A =  LH\RH;
    w1 = A*u;
     
    bx1 = cos(u+k/2*w1);    bx2 = -sin(u+k/2*w1);
    C = ConvectMat2D(p,t,bx1,bx2);
    S = SDAssembler2D(p,t,bx1,bx2); 
    LH = R+tau*C';
    RH = -(C+tau*S);
    A =  LH\RH;
    w2 = A*(u+k/2*w1);

    bx1 = cos(u+k/2*w2);    bx2 = -sin(u+k/2*w2);
    C = ConvectMat2D(p,t,bx1,bx2);
    S = SDAssembler2D(p,t,bx1,bx2); 
    LH = R+tau*C';
    RH = -(C+tau*S);
    A =  LH\RH;
    w3 = A*(u+k/2*w2);

    bx1 = cos(u+k*w3);      bx2 = -sin(u+k*w3);
    C = ConvectMat2D(p,t,bx1,bx2);
    S = SDAssembler2D(p,t,bx1,bx2); 
    LH = R+tau*C';
    RH = -(C+tau*S);
    A =  LH\RH;
    w4 = A*(u+k*w3);

    w = w1+2*w2+2*w3+w4;
    w(e(1,:)) = 0;
    u = u + k/6*w;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% RK4 %%%%%%%%%%%%%%
    time = time+k;
end

pdeplot(p,e,t,"XYData",u); 
pbaspect([1 1 1]);
caxis([0.5 1.3]);
colormap turbo;
saveas(gcf,"T22GLSh8.png");



