% getQdot returns the matrix right hand side of the equation dQ/dt=Q*S
%
% Inputs: Q,A,d
% Q - The d x d matrix Q(t)
% A - The d x d matrix A(t) = Df(t,u(t))
% d - The dimension of the problem
% p - the number of LEs being approximated
%
% Outputs: Qdot
% Qdot - The d x d matrix dQ/dt = Q*S
function Qdot = getQdot(Q,A,d,p)
  S = zeros(p,p);
  QtAQ = Q' * A * Q ;
  % This double loop computes S(Q)
  for i = 1:p
    for j = 1:p
        if i > j
            S(i,j) = QtAQ(i,j);
        elseif i < j
            S(i,j) = -QtAQ(j,i);
        else % i = j
            S(i,j) = 0;
        end       
    end
   end       
   % Qdot is the RHS of the Q'(t) = ... equation
   if d==p
       Qdot = Q*S; 
   else
       Qdot=A*Q-Q*(QtAQ-S);
   end
 %  B=QtAQ-Q'*Qdot;
end
