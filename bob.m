function [ T,X,tip ] = bob(work)

%This function file works the same as driv_les_mrc, but also uses the
%event_function need to check for tipping in the Compost Bomb problem

%[ applles,appules,stek,T,X,tip ] = bob(work)
% d is the dimension of the problem, 1 <=p<=d is the number of LEs and
% Steklov averages being approximated
d=2; p=2;
% H is the window-length of the Steklov averarges;
H = 5e-1;
% eps is the tolerance used in ode45
eps=1e-8;
% U0 is the initial conditions for the initial value problem
U0=zeros(d,1);
U0(1,1)=8.15;
U0(2,1)=50;
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
Ttransient=1;
prob=5; %code prob=3 is BiLin Example, prob=4 is UnLin Example
%work(1)=0.064; %for prob=5 (Compost Bomb) work(1) is epsilon %delta=0.5 for prob2-4
%work(2)=0.07; %for prob=5 (Compost Bomb) work(2) is v, the rate parameter 
                %a=rate parameter
%work(3)=0.5;
%Other parameters for Compost Bomb
CapitalPi=1.055;
r=0.01;
alpha=log(2.5)/10;
lambda=5.049e6;
BigA=3.9e7;
Tspan=[0 Tfinal];
options = odeset('RelTol',eps,'AbsTol',eps,'MaxStep',1e-2,'Events',@event_function);
% This line integrates for the Q variables and the u' = f(u,t) variables
TE=0;
tip=0;
[T,X,TE,XE,IE] = ode45(@(T,X) fullrhs(T, X, d,p,prob,work), Tspan, X0, options);
if TE>0
    tip=1;
end    
% applesfun approximates the upper and lower Lyapunov exponents
%[appules , applles] = applesfun(T,X,Ttransient ,prob,work,d,p);
% steklov approximates the Steklov averages
%stek = stekfun(T,X ,prob,work,H,d,p);

%This would make log scale pictures of solutions.
%figure
%     semilogy(T,work(2)*T+BigA*CapitalPi/lambda,'--')
%     s=num2str(work(2));
%     y=num2str(work(1));
%     title(strcat('epsilon=',y,',','v=',s))
%     hold on
%     semilogy(T,X(:,1))
%     hold off

end

