%close all;
%clear;

m = 401;
a = 1;
C = 2;
epsillon = 1e-1;
x_l = -1;
x_r = 1;
l = x_r-x_l;
h = l/(m-1);
x = x_l:h:x_r;
c1 = 1/2;
cr = 1.5;
T = 0.4;
%CFL = 7;
%CFL = 5.2;
CFL = 1.9;
k = CFL*h*h;
u0 = C-a*tanh(a*x/(2*epsillon));
v=u0';
vp1=v;
vp2=vp1;

eps = zeros(m,1);
alpha = 2/3;
SBP6_Variable;
L = [e_1';e_m'];
P = eye(m)-HI*L'*(L*HI*L')^(-1)*L;
B = HI*L'*(L*HI*L')^(-1);
t = 0;

while t<T
    Res = abs((v-4/3*vp1+1/3*vp2)/(2/3*k)+1/2*DD_3*diag(v)*v);
    Res = [Res(1);Res;Res(end)];
    if t+k>T
        k=T-t;
    end
    norm_para = norm(v-mean(v));
    wavespeed = v;
    wavespeed = [wavespeed(1);wavespeed;wavespeed(end)];
    adv = [v(1);v;v(end)];
    nv = zeros(m,1);
    for i=1:m
        nv(i) = abs(max([adv(i,1),adv(i+1,1),adv(i+2,1)])-min([adv(i,1),adv(i+1,1),adv(i+2,1)])-norm_para);
    end
    nv = [nv(1);nv;nv(end)];
    for i=1:m
        local_res = max([Res(i)/nv(i),Res(i+1)/nv(i+1),Res(i+2)/nv(i+2)]);
        local_ws = max([abs(wavespeed(i)),abs(wavespeed(i+1)),abs(wavespeed(i+2))]);
        eps(i)=min([c1*h*local_ws,cr*h^2*local_res]);
    end
    g = [gl(t,a,C,epsillon);gr(t,a,C,epsillon)];
    gt = [gpl(t,a,C,epsillon);gpr(t,a,C,epsillon)];
    w1=B*gt-alpha/2*P*D1*diag(P*v+B*g)*(P*v+B*g)-(1-alpha)*P*diag(P*v+B*g)*D1*(P*v+B*g)+epsillon*P*D2*(P*v+B*g)+eps.*P*D2*(P*v+B*g);
    g = [gl(t+k/2,a,C,epsillon);gr(t+k/2,a,C,epsillon)];
    gt = [gpl(t+k/2,a,C,epsillon);gpr(t+k/2,a,C,epsillon)];
    w2=B*gt-alpha/2*P*D1*diag(P*(v+k/2*w1)+B*g)*(P*(v+k/2*w1)+B*g)-(1-alpha)*P*diag(P*(v+k/2*w1)+B*g)*D1*(P*(v+k/2*w1)+B*g)+epsillon*P*D2*(P*(v+k/2*w1)+B*g)+eps.*P*D2*(P*(v+k/2*w1)+B*g);
    w3=B*gt-alpha/2*P*D1*diag(P*(v+k/2*w2)+B*g)*(P*(v+k/2*w2)+B*g)-(1-alpha)*P*diag(P*(v+k/2*w2)+B*g)*D1*(P*(v+k/2*w2)+B*g)+epsillon*P*D2*(P*(v+k/2*w2)+B*g)+eps.*P*D2*(P*(v+k/2*w2)+B*g);
    g = [gl(t+k,a,C,epsillon);gr(t+k,a,C,epsillon)];
    gt = [gpl(t+k,a,C,epsillon);gpr(t+k,a,C,epsillon)];
    w4=B*gt-alpha/2*P*D1*diag(P*(v+k*w3)+B*g)*(P*(v+k*w3)+B*g)-(1-alpha)*P*diag(P*(v+k*w3)+B*g)*D1*(P*(v+k*w3)+B*g)+epsillon*P*D2*(P*(v+k*w3)+B*g)+eps.*P*D2*(P*(v+k*w3)+B*g);
    vp2 = vp1;
    vp1 = v;
    v=v+k/6*(w1+2*w2+2*w3+w4);
    t=t+k;
end

u = C-a*tanh(a*(x-C*T)/(2*epsillon));
figure(2)
plot(x,u,"*");
hold on
plot(x,v);
legend("Analytic solution","FD solution",'Location','best');
%saveas(gcf,"T8ep000001.png");

e = u-v';
no = sqrt(h)*norm(e)

%0.0051 8.2865e-04 1.3854e-04

