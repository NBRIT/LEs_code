function xnew= rk_step(t,x,d,prob,h,Ark,b,c,s,par)
    g=zeros(s,d);
    for i=1:s
        if i ==1
            g(i,:)=frhs(t,x,d,prob,par);
        else
            g(i,:)=frhs(t+c(i)*h,x+h*Ark(i,1:i-1)*g(1:i-1,1),d,prob,par)'; 
        end
    end
    xnew=x+h*g(:,:)'*b;
end
