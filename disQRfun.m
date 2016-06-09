% Function disQRfun computes the Q and R factors and related quantities of a discrete QR iteration
% on a solution of u'=f(u,t)
%
% Inputs: sol,Q0,d,prob,work
%
% sol - A solution from a Matlab ode solver such as ode45 or ode15s
% Q0 - An initial orthogonal d x d matrix
% d - the dimension of the problem
% work - The vector of parameters that for u'=f(u,t)
%
% Outputs: Q, R, Rdiag, FMS
%
% Q - An N x d x d array of orthogonal factors from the discrete QR
% iteration where N is the number of time-steps of sol.x
% R - An N x d x d array of upper trianuglar factors from the discrete QR
% iteration 
% Rdiag - An N x d array of the diagonal elements of the upper triangular
% factors from the discrete QR iteration
% FMS - N x d x d array of the solution of X'=A(t)X where X(t_n)=I along
% the sequence of time-step of the solver (see discrete QR algorithm)
%
function [Q,R,Rdiag,FMS] = disQRfun(sol,Q0,d,prob,work)
    T=sol.x;
    X=(sol.y)';
    N=length(T);
   %  4th order dormand prince method for first 4 steps
 Ark=[0 0 0 0 0 0 0 ; 
          1/5 0 0 0 0 0 0; 
          3/40 9/40 0 0 0 0 0 ; 
          44/45 -56/15 32/9 0 0 0 0 ; 
          19372/6561 -25360/2187 64448/6561 -212/729 0 0 0 ;
          9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0 ;
          35/384 0 500/1113 125/192 -2187/6784 11/84 0];
    % b1 is 4th order accurate
          b1=[5179/57600 ; 0 ; 7571/16695 ; 393/640 ; -92097/339200 ; 187/2100 ; 1/40];
    % b2 is 5th order accurate
          b2=[ 35/384 ; 0 ; 500/1113 ; 125/192 ; -2187/6784; 11/84; 0 ];
          c=[0 ; 1/5; 3/10; 4/5; 8/9; 1 ; 1];
          s=7;
    Qtemp=Q0; R=zeros(N-1,d,d); Rdiag=zeros(N-1,d); Q=zeros(N-1,d,d);
    FMS=zeros(N-1,d,d);
    for n=1:N-1
        h=T(n+1)-T(n); Xn=zeros(d,d); 
        for j=1:d
            e=zeros(d,1); e(j)=1;
        % xnew1 is the 4th order update, xnew2 is the 5th order update
            [xnew1,xnew2]=rk_step(sol,T(n),X(n,:)',e,d,prob,h,Ark,b1,b2,c,s,work);
          %   Xn(1:d,j)=
             Xn(1:d,j)=xnew1(1:d,1);
        end
        T(n)
     %   getA(T(n),X(n,:)',zeros(d,1),d,prob,work)
        [Qtemp,Rtemp]=qr(Xn*Qtemp);
        for j=1:d
            if Rtemp(j,j) < 0
                Rtemp(j,:)=-1*Rtemp(j,:);
                P=eye(d);
                P(j,j)=-1;
                Qtemp=Qtemp*P;
            end
        end
        R(n,1:d,1:d)=Rtemp(1:d,1:d);
        for j=1:d
            Rdiag(n,j)=Rtemp(j,j);
        end
        Q(n,1:d,1:d)=Qtemp(1:d,1:d);
        FMS(n,1:d,1:d)=Xn(1:d,1:d);
    end
   
end