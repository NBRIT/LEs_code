% Function rhs returns a vector that corresponds to the full equation that
% computes the solution of u'=f(t,u), Q'=Q*S, u(0)=u_0, Q(0)=Q_0
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
% rhs - A (d+d*d) x 1 vector with du/dt = f(t,u)=rhs(1:d) and dQ/dt = Q*S = reshape(Qdot,d,d)
%
% Problems: 
% prob==1 : work(4) is the rate at which the equilibrium changes
function rhs=fullrhs(t,x,d,p,prob,work)  
    rhs=zeros(d,1); 
    rhs(1:d,1)=frhs(t,x(1:d),d,prob,work);
    A=getA(t,x(1:d),rhs(1:d,1),d,prob,work);
    Q=reshape(x(d+1:d+d*p),d,p);
    [Q,junk]=mgs(Q);
    Qdot=getQdot(Q,A,d,p);
    rhs(d+1:d+d*p,1)=reshape(Qdot,d*p,1);
end
