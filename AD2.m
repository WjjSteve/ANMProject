close all;
clear;

g = @circleg;
[p,e,t] = initmesh(g,"Hmax",1/4);
% initialize the mesh

T = 1;
time = 0;
k = 0.01;
% time interval

eps = 0;
TOL = 5e-3;

RES = [];
NT = [];
while time < T

    %[up,u,kl] = upwardGFEM(p,e,t,time+k,eps);
    %z = backwardGFEM(p,e,t,time+k,eps);
    [up,u,kl] = upwardRV(p,e,t,time+k);
    z = backwardRV(p,e,t,time+k);
    n = size(p,2);
    nt = size(t,2);
    bx1 = -2*pi*p(2,:);
    bx2 = 2*pi*p(1,:); 
    R = ReactMat2D(p,t);
    C = ConvectMat2D(p,t,bx1,bx2);
    LH = R;
    RH = R*(u-up)/kl+C*u;
    Res = LH\RH;
    C = ConvectMat2D(p,t,ones(1,n),ones(1,n));
    RH = C*z;
    Dz = LH\RH;
   
    etas = zeros(nt,1);
    for K = 1:nt
       loc2glb = t(1:3,K); 
       x = p(1,loc2glb); 
       y = p(2,loc2glb); 
       area = polyarea(x,y); 
       MK = [2 1 1;
             1 2 1;
             1 1 2]/12*area; 
       h = ComputeDiameter(x,y)/2;
       etas(K) = h*sqrt(Res(loc2glb)'*MK*Res(loc2glb))*sqrt(Dz(loc2glb)'*MK*Dz(loc2glb));
    end

    RES = [RES,sum(etas)];
    NT = [NT,nt];

    if T*sum(etas)<TOL
        break;
    end

    figure
    pdemesh(p,e,t);
    [maxetas,elements] = maxk(etas,round(0.1*nt));
    [p,e,t] = refinemesh(g,p,e,t,elements,'longest');
    time = time+k;
end

pdeplot(model,"XYData",u,"ZData",u)
pbaspect([1 1 1]);
caxis([-15 0]);
colormap turbo;
plot(NT,RES);
pdemesh(p,e,t);
saveas(gcf,"p3f2e0me.png");
pdemesh(p,e,t)
saveas(gcf,"p3f1e01me.png");

loglog(NT1,RES1,'g-s');%,'color','#5ed1fd');
hold on
loglog(NT2,RES2,'b-o');%,'color');%'#8a5efd');
hold on
loglog(NT1,NT1.^-1,'--');%,'color');%,'#808080');
xlabel('log_{10}#nodes');
ylabel('log_{10}\eta');
title('Error indicators \eta');
legend('Traget1','Traget2','N^{-1}');
saveas(gcf,"p3e0.png");