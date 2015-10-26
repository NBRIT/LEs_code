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
    elseif prob==2
         del=work(1);
         b0=1e-6;
         b1=(1-b0)/b0;
         a=work(2);
         b=exp(a*t)/(b1+exp(a*t));      
         A(1,1)=-((u(1)-b)*(u(1)-b-work(3)*del)+(u(1)-b)*(u(1)-b+del)+(u(1)-b-work(3)*del)*(u(1)-b+del));
    end
  
end