% function Qanglefun returns the
%
% inputs:  sol1,sol2,d,p,number
%
% sol1,sol2 - the solution structure returned by a Matlab ODE solver
% d - the dimension of the problem
% p - the number of times at which the angle is compute
%
% outputs: Tangle, Qangle
%
% Tangle - Vector of time steps at which the angles are compute
% Qangle - Angle between the columns of Q computed with two different
% paramter sets
%
function [Tangle,Qangle] = Qanglefun( sol1,sol2,d,p,number)
    Tangle=zeros(number,1); Qangle=zeros(number,p);
    T=sol1.x;
    Tangle(1)=T(1); dt=(T(end)-T(1))/(number);
    for n=1:number-1
        Tangle(n+1)=Tangle(n)+dt;
        X1=deval(sol1,Tangle(n+1));
        X2=deval(sol2,Tangle(n+1));
        Q1=reshape(X1(d+1:end),d,p);
        [Q1,junk]=mgs(Q1);
        Q2=reshape(X2(d+1:end),d,p);
        [Q2,junk]=mgs(Q2);
        for j=1:p
            Qangle(n+1,j)= subspace(Q1(:,j),Q2(:,j));
        end
    end
end

