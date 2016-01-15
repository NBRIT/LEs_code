%checking for a critical v value for different epsilons in the compost
%bomb. This file runs the for loops, then calculates the minimum value of v
%for which there is tipping (within a certain time limit), and then creates
%a plot of epsilon v. vcrit

%file uses the bob.m function file (which uses the event_function.m file)

close all
vvalues=0.07:0.0005:0.12;
epvalues=logspace(-4,-1,30);
tiparray=zeros(length(epvalues),length(vvalues));
for i=1:length(epvalues)
    work=zeros(1,2);
    work(1)=epvalues(i);
for n=1:length(vvalues)
    work(2)=vvalues(n);
    %[applles,appules,stek,T,X,tip]=bob(work);
    [T,X,tip]=bob(work);
    tiparray(i,n)=tip;
    [i,n]
end
end
tiparray
%%
tippy=tiparray;
vcrit=zeros(length(epvalues),1);
for i=1:length(epvalues)
    k=find(tippy(i,:));
   if isempty(k)
       vcrit(i)=1
   else vcrit(i)=vvalues(min(k));
   end
end
vcrit
%picture:
semilogx(epvalues(1:length(epvalues)-1),vcrit(1:length(epvalues)-1))
