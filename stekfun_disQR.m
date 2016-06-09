% stekfun_disQR returns the Steklov averages from a discrete QR iteration
%
% Inputs: solR,H,p,number
%
% solR - A Matlab structure where solR.x is the time-sequence of the
% solution of u'=f(u,t) and solR.R is the diagonal of the upper triangular
% factors of the discrete QR algorithm
% H - The window length of the Steklov averages
% p - the dimension of the upper triangular factors R (usuallly equal to d)
% number - the number of points at which the Steklov averages are computed
%
% Outputs: Tstek,stek
% Tstek - vector of the times at which the Steklov averages are computed
% stek - p x 1 vector of Steklov averages
function [Tstek,stek]=stekfun_disQR(solR,H,p,number) 
    mesh=20;
    dh=H/mesh;
    T=solR.x;
    t0=T(1);
    dt=(T(end)-H-T(1))/(number-1);
    Bdiag=zeros(mesh,p);
    Tstek=zeros(number,1);
    stek=zeros(number,p);
    for n=1:number
        th=t0;
        Ttemp=zeros(mesh,1);
        for m=1:mesh            
            Rdiag=deval(th,solR);
            for j=1:p
                Bdiag(m,j)=log(Rdiag(j));
            end
            Ttemp(m,1)=th;
            th=th+dh;
        end
        for j=1:p
            stek(n,j) = trapz(Ttemp,Bdiag(:,j))/H;
        end
        Tstek(n)=t0;
        t0=t0+dt;
    end
end
    
    
    