% Function getA forms the Jacobian along a solution trajectory at a given
% point (t,u)
%
% Inputs: t,u,d,prob,work
% t - A scalar that represents the time
% u - A d x 1 vector that represents the value the solution u(t) at time t
% d - The dimension of the phase space
% prob - A scalar number that corresponds to the ODE being solved (used to change between different
% problems
% work - A vector of parameters of arbitrary size
%
% Outputs: A
% A - The d x d Jacobian matrix A(t)=Df(t,u(t))
%
% Problems: 
% prob==1 : work(4) is the rate at which the equilibrium changes
function A=getA(t,u,up,d,prob,work)  
    A=zeros(d,d);
    if prob==1 
         % Lorenz '63 
         sigma = 10;
         beta=8/3;
         rho=28;
         J=zeros(d,d);
         J(1,1)=-sigma;
         J(1,2)=sigma;
         J(2,1)=rho-u(3);
         J(2,2)=-1;
         J(2,3)=-u(1);
         J(3,1)=u(2);
         J(3,2)=u(1);
         J(3,3)=-beta;
         A=J;
    elseif prob==2 % Bistable Logistic
         del=work(1);
         b0=1e-6;
         b1=(1-b0)/b0;
         a=work(2);
         b=exp(a*t)/(b1+exp(a*t));      
         A(1,1)=-((u(1)-b)*(u(1)-b-work(3)*del)+(u(1)-b)*(u(1)-b+del)+(u(1)-b-work(3)*del)*(u(1)-b+del));
    elseif prob==3 %% this is the Unique Logistic Code
         del=work(1);
         b0=1e-6;
         b1=(1-b0)/b0;
         a=work(2);
         b=exp(a*t)/(b1+exp(a*t));
         A(1,1)=-2*(u(1)-b)-del;
    elseif prob ==4 %% this is the Bistable Linear Code
        del=work(1);
        b0=0;
        a=work(2);
        b=a*t+b0;
        A(1,1)=-(u(1)-b)*(u(1)-b-del)-(u(1)-b)*(u(1)-b+del)-(u(1)-b-del)*(u(1)-b+del);
    elseif prob==5 %% this is the Unique Linear Code
         del=work(1);
         b0=0;
         a=work(2);
         b=a*t+b0;
         A(1,1)=-2*(u(1)-b)-del;
    elseif prob == 6
        a=work(2);
        b0=0;
        b=a*t+b0;
        epsilon=work(1);
        CapitalPi=1.055;
        r=0.01;
        alpha=log(2.5)/10;
        lambda=5.049e6;
        BigA=3.9e7;
        J=zeros(d,d);
        J(1,1)=(1/epsilon)*(alpha*u(2)*r*exp(alpha*u(1))-lambda/BigA);
        J(1,2)=(1/epsilon)*(r*exp(alpha*u(1)));
        J(2,1)=-alpha*u(2)*r*exp(alpha*u(1));
        J(2,2)=-r*exp(alpha*u(1));
        A=J;
        f(1)= (1/epsilon)*(u(2)*r*exp(alpha*u(1))-(lambda/BigA)*(u(1)-b));
        f(2)=CapitalPi-u(2)*r*exp(alpha*u(1));
    elseif prob==7
         del=work(1);
         b0=1e-6;
         b1=(1-b0)/b0;
         a=work(2);
         b=exp(a*t)/(b1+exp(a*t));      
         A(1,1)=-((u(1)-b)*(u(1)-b-work(3)*del)+(u(1)-b)*(u(1)-b+del)+(u(1)-b-work(3)*del)*(u(1)-b+del));
         A(2,1)= u(1);
         A(1,2)=0;
         A(2,2)=-1;
    end
end
