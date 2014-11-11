clear all
close all
clc
B=[1 2;2 1];
C=[0 1;1 0];
Cu=B;

x=[.01;.01];
RHS=[1;2];


itermax=1600;
iter=1;
while(iter<=1)
    
%     [df,jf]=jaco(Cu,x,RHS);
%     [u,s,v]=svd(df);
%     ds=s;
%     for i=1:length(s);
%         if(s(i,i)~=0)
%         ds(i,i)=1/s(i,i);
%         end
%     end
    [df,jf]=jaco(Cu,x,RHS);
    x=x-inv(df)*panelEQ(Cu,x,RHS);
    if(panelEQ(Cu,x,RHS)<=(10^-16))
        break;
    end
    iter=iter+1;
end