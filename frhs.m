% Function frhs forms the right hand side of du/dt = f(t,u(t))
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
% f- A d x 1 vector that represents the right hand side of du/dt = f(t,u(t))
%
% Problems: 
% prob==1 : Lorenz '63 with sigma=10; beta = 8/3; rho =28;
% prob==2 : NBRIT example; work(1) is the distance the QSE are apart and work(2) is the rate at which the equilibrium changes
% , and work(3) is a parameter to break the symetry of the distance the QSE
% are apart
function f=frhs(t,u,d,prob,work)
    f=zeros(d,1);

    if prob==1
     % Lorenz '63    
        sigma = 10;
        beta=8/3;
        rho=28;
        f(1)=sigma*(u(2)-u(1));
        f(2)=u(1)*(rho-u(3))-u(2);
        f(3)=u(1)*u(2)-beta*u(3);
    elseif prob==2
        del=work(1);
         b0=1e-6;
         b1=(1-b0)/b0;
         a=work(2);
         b=exp(a*t)/(b1+exp(a*t));      
         f(1)=-(u(1)-b)*(u(1)-b-del)*(u(1)-b+work(3)*del);
    end
 
end