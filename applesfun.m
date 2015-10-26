% getQdot returns the matrix right hand side of the equation dQ/dt=Q*S
%
% Inputs: T, X, Ttransient, prob, work, d , p
% T - vector of times
% X - vector of solutions of u'=f(t,u) and Q'=S*Q
% Ttransient - time at which the limsup approximation beings
% prob - A scalar number that corresponds to the ODE being solved (used to change between different
% problems
% work - A vector of parameters of arbitrary size
% d - The dimension of the problem
% p - the number of LEs being approximated
%
% Outputs: appules, applles
% appules - p x 1 vector of the approximate upper Lyapunov exponents
% applles - p x 1 vector of the apporximate lower Lypaunov exponents
function [appules,applles] = applesfun(T,X,Ttransient ,prob,work,d,p)
    Tn=length(T); 
    u=zeros(d,1);
    Bint=zeros(p,1); appules=zeros(p,1); applles=zeros(p,1);
    [t_delta, t_ind] = min(abs(T-Ttransient));
    t_ind=t_ind+1;
    g=zeros(Tn-t_ind-1,d);
    for n=1:Tn-1
        t=T(n); h=T(n+1)-T(n);
        u(1:d)=X(n,1:d);
        up=frhs(t,u,d,prob,work);
        A=getA(t,u,up,d,prob,work) ;
        Q=reshape(X(n,d+1:d+d*p),d,p);
        QtAQ=Q'*A*Q;
        S=zeros(p,p);
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
        B=triu(QtAQ-S);
        for j=1:p
            Bint(j)=Bint(j)+h*B(j,j);
            if n > t_ind             
                g(n-t_ind,j)=Bint(j)/t;
            end
        end
  
        
    end       
    for j=1:p
        appules(j)=max(g(:,j)    );
        applles(j)=min(g(:,j));
    end
end

