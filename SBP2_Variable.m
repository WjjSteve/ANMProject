%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2nbd order SBP Finite differens         %%%
%%% operatorers by Ken Mattsson             %%%
%%% Here with variable coeff. second        %%%
%%% derivative (c(x)*u_x)_x, and also       %%%
%%% including the building block for AD     %%%
%%%                                         %%%
%%% Variable c(x) needed as a vector with   %%%
%%% m komponents to build D2                %%% 
%%%                                         %%%
%%% Theses OP will be useful in first part  %%%
%%% of project 3, in course ANM             %%%
%%%                                         %%%
%%% H      The norm                         %%%
%%% D1=   (First deivative)                 %%%
%%% D2=   (Variable second derivative)      %%%
%%%                                         %%%
%%% DD_1 undivided first derivative FD OP.  %%%
%%%                                         %%%
%%% Add DI_2= -gamma_2*HI*DD_1'*DD_1 to RHS %%%
%%% for some suitable value on gamma_2 > 0  %%%
%%%                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%m=11; %problemstorlek
%h=1/(m-1);
%n=(m-1)/2;
%x=-n*h:h:n*h;x=x';

%x0=x.^0/fac(1); 
%x1=x.^1/fac(1);
%x2=x.^2/fac(2);
%x3=x.^3/fac(3);
%x4=x.^4/fac(4);

c=ones(m,1);

D1=((.5*diag(ones(m-1,1),1)-.5*diag(ones(m-1,1),-1)));
D1(1,1)=-1;D1(1,2)=1;D1(m,m-1)=-1;D1(m,m)=1;
D1(m,m-1)=-1;D1(m,m)=1;
D1=D1/h;

e_1=zeros(m,1);e_1(1)=1;
e_m=zeros(m,1);e_m(m)=1;

S_U=[-3/2, 2, -1/2]/h;
d_1=zeros(1,m);
d_1(1:3)=S_U;
d_m=zeros(1,m);
d_m(m-2:m)=fliplr(-S_U);

H=(eye(m,m));H(1,1)=0.5;H(m,m)=0.5;
H=h*H;
HI=inv(H);

M=zeros(m,m);

for i=2:m-1
M(i,i-1:i+1)=[-c(i-1) / 0.2e1 - c(i) / 0.2e1 c(i-1) / 0.2e1 + c(i) + c(i+1) / 0.2e1 -c(i) / 0.2e1 - c(i+1) / 0.2e1;];
end

M(1:2,1:2)=[c(1) / 0.2e1 + c(2) / 0.2e1 -c(1) / 0.2e1 - c(2) / 0.2e1; -c(1) / 0.2e1 - c(2) / 0.2e1 c(1) / 0.2e1 + c(2) + c(3) / 0.2e1;];
M(m-1:m,m-1:m)=[c(m-2) / 0.2e1 + c(m-1) + c(m) / 0.2e1 -c(m-1) / 0.2e1 - c(m) / 0.2e1; -c(m-1) / 0.2e1 - c(m) / 0.2e1 c(m-1) / 0.2e1 + c(m) / 0.2e1;];
M=M/h;

D2=HI*(-M-diag(c)*e_1*d_1+diag(c)*e_m*d_m);

%AD
DD_1=(diag(ones(m-1,1),+1)-diag(ones(m,1),0));
DD_1(m,m-1:m)=[-1 1];



