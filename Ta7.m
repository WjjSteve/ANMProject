close all;
clear;

m = 200;
a = 1;
C = 2;
epsillon = 1e-2;
x_l = -1;
x_r = 1;
l = x_r-x_l;
h = l/(m-1);
x = x_l:h:x_r;
gamma = 1;

T = 0.4;
CFL = 7;
%CFL = 5.2;
%CFL = 1.9;
k = CFL*h*h;
u0 = C-a*tanh(a*x/(2*epsillon));
v=u0';
figure(1)
plot(x,v);

alpha = 2/3;
SBP2_Variable;
L = [e_1';e_m'];
P = eye(m)-HI*L'*(L*HI*L')^(-1)*L;
B = HI*L'*(L*HI*L')^(-1);
t = 0;
while t<T
    if t+k>T
        k=T-t;
    end
    g = [gl(t,a,C,epsillon);gr(t,a,C,epsillon)];
    gt = [gpl(t,a,C,epsillon);gpr(t,a,C,epsillon)];
    w1=B*gt-alpha/2*P*D1*diag(P*v+B*g)*(P*v+B*g)-(1-alpha)*P*diag(P*v+B*g)*D1*(P*v+B*g)+epsillon*P*D2*(P*v+B*g)-gamma*P*HI*DD_1'*DD_1*(P*v+B*g);
    g = [gl(t+k/2,a,C,epsillon);gr(t+k/2,a,C,epsillon)];
    gt = [gpl(t+k/2,a,C,epsillon);gpr(t+k/2,a,C,epsillon)];
    w2=B*gt-alpha/2*P*D1*diag(P*(v+k/2*w1)+B*g)*(P*(v+k/2*w1)+B*g)-(1-alpha)*P*diag(P*(v+k/2*w1)+B*g)*D1*(P*(v+k/2*w1)+B*g)+epsillon*P*D2*(P*(v+k/2*w1)+B*g)-gamma*P*HI*DD_1'*DD_1*(P*(v+k/2*w1)+B*g);
    w3=B*gt-alpha/2*P*D1*diag(P*(v+k/2*w2)+B*g)*(P*(v+k/2*w2)+B*g)-(1-alpha)*P*diag(P*(v+k/2*w2)+B*g)*D1*(P*(v+k/2*w2)+B*g)+epsillon*P*D2*(P*(v+k/2*w2)+B*g)-gamma*P*HI*DD_1'*DD_1*(P*(v+k/2*w2)+B*g);
    g = [gl(t+k,a,C,epsillon);gr(t+k,a,C,epsillon)];
    gt = [gpl(t+k,a,C,epsillon);gpr(t+k,a,C,epsillon)];
    w4=B*gt-alpha/2*P*D1*diag(P*(v+k*w3)+B*g)*(P*(v+k*w3)+B*g)-(1-alpha)*P*diag(P*(v+k*w3)+B*g)*D1*(P*(v+k*w3)+B*g)+epsillon*P*D2*(P*(v+k*w3)+B*g)-gamma*P*HI*DD_1'*DD_1*(P*(v+k*w3)+B*g);
    v=v+k/6*(w1+2*w2+2*w3+w4);
    t=t+k;
end

u = C-a*tanh(a*(x-C*T)/(2*epsillon));
figure(2)
plot(x,u,"*");
hold on
plot(x,v);
legend("Analytic solution","FD solution",'Location','best');
%saveas(gcf,"T72Ae0.png");

e = u-v';
no = sqrt(h)*norm(e)

%7.9532e-05 1.2952e-06 2.5190e-08
%1.7255e-04 1.2470e-05 1.1056e-06
%0.1729 0.1241 0.0878
%0.0698 0.0493 0.0349
%0.0807 0.0572 0.0404
% 0.1300 0.0754
% 0.0501 0.0146
% 0.0437 0.0091