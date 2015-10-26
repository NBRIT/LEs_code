% d is the dimension of the problem, 1 <=p<=d is the number of LEs and
% Steklov averages being approximated
d=1; p=1;
% H is the window-length of the Steklov averarges;
H = 5e-1;
% eps is the tolerance used in ode45
eps=1e-8;
% U0 is the initial conditions for the initial value problem
U0=zeros(d,1);
U0(1,1)=1;
% Q0 is the initial condition for the Q0, can use identity matrix or some
% random matrix
Q0=zeros(d,p);
for j=1:p
    Q0(j,j)=1;
end
Q0=reshape(Q0,d*p,1);
X0 = [U0 ; Q0 ] ;
% Tfinal is the final time
Tfinal=25;
% Ttransient is the time of the transient behavior; used in the limsup
% calculation
Ttransient=10;
prob=2;
work(1)=0.5;
work(2)=0.3;
work(3)=0.5;
Tspan=[0 Tfinal];
options = odeset('RelTol',eps,'AbsTol',eps);
% This line integrates for the Q variables and the u' = f(u,t) variables
[T,X] = ode45(@(T,X) fullrhs(T, X, d,p,prob,work), Tspan, X0, options);

% applesfun approximates the upper and lower Lyapunov exponents
[appules , applles] = applesfun(T,X,Ttransient ,prob,work,d,p);
% steklov approximates the Steklov averages
stek = stekfun(T,X ,prob,work,H,d,p);
 
 
 
 