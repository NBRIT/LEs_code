% d is the dimension of the problem, 1 <=p<=d is the number of LEs and
% Steklov averages being approximated
d=2; p=2;
% H is the window-length of the Steklov averarges;
H = 1;
% eps is the tolerance used in ode45
eps=1e-6;
% U0 is the initial conditions for the initial value problem
U0=zeros(d,1);
U0(1,1)=0.5;
U0(2,1)=.5;
% Q0 is the initial condition for the Q0, can use identity matrix or some
% random matrix
%Q0=zeros(d,p);
Q0=zeros(d,p);
theta=pi/4;
Q0(1,1)=cos(theta); Q0(1,2)=sin(theta);
Q0(2,1)=-sin(theta); Q0(2,2)=cos(theta);
Q00=Q0;
Q0=eye(d,p);
Q0=reshape(Q0,d*p,1);
X0 = [U0 ; Q0 ] ;
% Tfinal is the final time
Tfinal=120;
% Ttransient is the time of the transient behavior; used in the limsup
% calculation
Ttransient=Tfinal/2;
prob=7;
work(1)=0.5;
work(2)=0.376;
work(3) = 1;
Tspan=[0 Tfinal];
options = odeset('RelTol',eps,'AbsTol',eps);
% This line integrates for the Q variables and the u' = f(u,t) variables
sol = ode45(@(T,X) fullrhs(T, X, d,p,prob,work), Tspan, X0, options);
T=sol.x;
X=(sol.y)';
% applesfun approximates the upper and lower Lyapunov exponents
[appules , applles] = applesfun(T,X,Ttransient ,prob,work,d,p);
% steklov approximates the Steklov averages
inc=1;
[Tstek,stek] = stekfun(T,X ,prob,work,H,d,p,inc);

work(2)=.377;
sol2 = ode45(@(T,X) fullrhs(T, X, d,p,prob,work), Tspan, X0, options);
T2=sol2.x;
X2=(sol2.y)';
% applesfun approximates the upper and lower Lyapunov exponents
[appules2 , applles2] = applesfun(T2,X2,Ttransient ,prob,work,d,p);
% steklov approximates the Steklov averages
inc=1;
[Tstek2,stek2] = stekfun(T2,X2 ,prob,work,H,d,p,inc);


work(2)=.378;
sol3 = ode45(@(T,X) fullrhs(T, X, d,p,prob,work), Tspan, X0, options);
T3=sol3.x;
X3=(sol3.y)';
% applesfun approximates the upper and lower Lyapunov exponents
[appules3 , applles3] = applesfun(T2,X2,Ttransient ,prob,work,d,p);
% steklov approximates the Steklov averages
inc=1;
[Tstek3,stek3] = stekfun(T3,X3 ,prob,work,H,d,p,inc);
Q3=zeros(length(T3),d,p); fpp3=zeros(length(T3),d); qse3=zeros(length(T3),3);

% These function calls compute the angle between the columns of Q for
% between the (1st and 2nd ) and the (2nd and 3rd ) paramter sets
[Tangle,Qangle] = Qanglefun( sol,sol2,d,p,100);
[Tangle2,Qangle2] = Qanglefun( sol2,sol3,d,p,100);
