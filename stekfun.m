% stekfun returns
%
% Inputs: T, X, Ttransient, prob, work, d , p
% T - vector of times
% X - vector of solutions of u'=f(t,u) and Q'=S*Q
% prob - A scalar number that corresponds to the ODE being solved (used to change between different
% problems
% work - A vector of parameters of arbitrary size
% H - the window length of the Steklov averages
% d - The dimension of the problem
% p - the number of LEs being approximated
% increment - The number of time steps in T which are between Steklov
% average calculations
%
% Outputs: stek
% stek - p x 1 vector of Steklov averages
function [Tstek,stek] =stekfun(T,X,prob,work,H,d,p,increment) 
    u=zeros(d,1);
    mfin=find(T < T(end)-H,1,'last');
    stek=zeros(mfin,p);
    Tstek=zeros(mfin,1);
    for counter=1:mfin
        tprev=T(counter);
        m=find(T < tprev+H,1,'last');
        tnew=T(m);
        for n=counter:increment:m    
            Bint=zeros(p,1);
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
            end
        end
        for j=1:p
            stek(counter,j)=Bint(j)/H;
            Tstek(counter,j)=tprev;
        end
    end
end
