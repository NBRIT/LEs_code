% applesfun_disQR returns the approximate upper and lower LEs
%
% Inputs: T,R,Ttransient, d, p
% T - vector of times
% R - vector of upper triangular factors from the discrete QR iteration
% Ttransient - time at which the limsup approximation beings
% d - The dimension of the problem
% p - the number of LEs being approximated
%
% Outputs: appules, applles
% appules - p x 1 vector of the approximate upper Lyapunov exponents
% applles - p x 1 vector of the apporximate lower Lypaunov exponents
% apples - N x p vector of the value of the approximate LE along the
% interval of the length of T
function [appules,applles,apples] = applesfun_disQR(T,R,Ttransient,d,p)
    apples=zeros(length(T),p); appules=zeros(p,1); applles=zeros(p,1); 
    for n=1:length(T)-1
        for j=1:p
          apples(n+1,j)=(T(n)*apples(n,j)+log(R(n,j,j)))/T(n+1);
        end
    end
     [c index] = min(abs(T-Ttransient));
     for j=1:p
         appules(j)=max(apples(index:end,j));
         applles(j)=min(apples(index:end,j));
     end
end
