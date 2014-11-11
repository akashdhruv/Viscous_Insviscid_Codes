function [df,Jf]=jaco(A,var,Bt,sig,RHS,xj_d,Hi,cf_2,Naw,dc_Hstar,Hstari)
m = sym(zeros(size(var)));
for k = 1:numel(var)
    m(k) = sym(sprintf('m%d', k));
end
f=panelEQ(A,m,Bt,sig,RHS,xj_d,Hi,cf_2,Naw,dc_Hstar,Hstari);
Jf=jacobian(f,m);
df=subs(Jf,m,var);
end