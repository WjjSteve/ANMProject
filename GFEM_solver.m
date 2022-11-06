% This is the code that implement the GFEM
close all;
clear;

g = Rectg(-2,-2.5,2,1.5);
hmax = 1/16;
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
time = 0;
while time<T
    bx1 = cos(u);
    bx2 = -sin(u);
    % vector field f'(u)

    C = ConvectMat2D(p,t,bx1,bx2);
    nor = max(sqrt(bx1.^2+bx2.^2));
    k = CFL*hmax/nor;

    % Time steps
    if time+k>T
        k=T-time;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    LH = R/k+C/2;
    RH = (R/k-C/2)*u;
    I = eye(length(p));
    LH(e(1,:),:) = I(e(1 ,:) ,:); 
    RH(e(1 ,:)) = pi/4;
    u = LH\RH;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Crank-Nicolson %%%%%%%%
    time = time+k;
end

pdeplot(p,e,t,"XYData",u); 
pbaspect([1 1 1]);
caxis([0.5 1.3]);
colormap turbo;
saveas(gcf,"T21GFEMh16.png");
