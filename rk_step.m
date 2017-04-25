function [xnew1,xnew2]= rk_step(t,x,d,prob,h,Ark,b1,b2,c,s,par)
    g=zeros(s,d);
    for i=1:s
        if i ==1
            g(i,:)=frhs(t,x,d,prob,par);
        else
            g(i,:)=frhs(t+c(i)*h,x+h*Ark(i,1:i-1)*g(1:i-1,1),d,prob,par)'; 
        end
    end
    xnew1=x+h*g(:,:)'*b1;
    xnew2=x+h*g(:,:)'*b2;
end
