% d is the dimension of the problem, 1 <=p<=d is the number of LEs and
% Steklov averages being approximated
d=3; p=3;
% H is the window-length of the Steklov averarges;
H = 0.005;
% eps is the tolerance used in ode45
eps=1e-8;
% U0 is the initial conditions for the initial value problem
U0=zeros(d,1);
U0(1,1)=8.15;
U0(2,1)=50;
U0(3,1)=0;
% Q0 is the initial condition for the Q0, can use identity matrix or some
% random matrix
Q0=zeros(d,p);
QQ0=zeros(d,d);
QQ0(1,1)=2; QQ0(1,2)=-2; QQ0(1,3)=1;
QQ0(2,1)=1; QQ0(2,2)=2; QQ0(2,3)=2;
QQ0(3,1)=2; QQ0(3,2)=1; QQ0(3,3)=-2;
QQ0=QQ0/3;
Q0(1:d,1:p)=QQ0(1:d,1:p);

[Q0, junk]=mgs(Q0);
X0(1:d,1)=U0(1:d,1);
X0(d+1:d+d*p,1)=reshape(Q0,d*p,1);
% Tfinal is the final time
% Tfinal=22.24;
Tfinal=.1;
% Ttransient is the time of the transient behavior; used in the limsup
% calculation
Ttransient=Tfinal/2;
prob=5;
work(1)=1.0;   %.064;
work(2)=0.0;
work(2)=.375;
work(3)=.1148;
Tspan=[0 Tfinal];
options = odeset('RelTol',eps,'AbsTol',eps);
% This line integrates for the u' = f(u,t) variables
tic
sol= ode15s(@(t,x) frhs(t, x, d,prob,work), Tspan, U0,options);
toc
tic
% this line computes the discrete QR iteration
[Q,R,Rdiag,FMS] = disQRfun(sol,Q0,d,prob,work);
T=sol.x;
% this line computes the upper and lower LEs and the approximate LEs along
% the length of Tspan
[appules,applles,apples] = applesfun_disQR(T,R,Ttransient,d,p);
toc
T=sol.x;
Xd=(sol.y)';
Q0=reshape(Q0,d,p);

solR=sol;
solR.y=(Rdiag)';
% this line computes the Steklov averages
[Tstek,stek]=stekfun_disQR(sol,H,p,100); 
 work(1)=1.0; %.064;
 work(2)=0.0;
 work(2)=.375;
 work(3)=0.1158;

 options = odeset('RelTol',eps,'AbsTol',eps);
tic
sol2= ode15s(@(t,x) frhs(t, x, d,prob,work), Tspan, U0,options);
toc
tic
T2=sol2.x;
[Q2,R2,Rdiag2] = disQRfun(sol2,Q0,d,prob,work);
[appules2,applles2,apples2] = applesfun_disQR(T2,R2,Ttransient,d,p);
toc

Xd2=(sol2.y)';
Q0=reshape(Q0,d,p);

solR=sol;
solR.y=(Rdiag2)';
[Tstek2,stek2]=stekfun_disQR(sol2,H,p,100); 


solR=sol;
solR.y=(Rdiag)';
[Tstek,stek]=stekfun_disQR(sol,H,p,100); 
 work(1)=1.0; %.064;
 work(2)=0.0;
 work(2)=.375;
 work(3)=.2;

 options = odeset('RelTol',eps,'AbsTol',eps);
tic
sol3= ode15s(@(t,x) frhs(t, x, d,prob,work), Tspan, U0,options);
toc
tic
T3=sol3.x;
[Q3,R3,Rdiag3] = disQRfun(sol3,Q0,d,prob,work);
[appules3,applles3,apples3] = applesfun_disQR(T3,R3,Ttransient,d,p);
toc

Xd3=(sol3.y)';
Q0=reshape(Q0,d,p);

solR=sol;
solR.y=(Rdiag3)';
[Tstek3,stek3]=stekfun_disQR(sol3,H,p,100); 
