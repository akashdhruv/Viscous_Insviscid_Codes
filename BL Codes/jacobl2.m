function [df,Jf]=jacobl2(theta,Hi,Uej,cf_2,xj_d,dc_Hstar,Hstari,Naw)
s = sym(zeros(size(Uej)));
for k = 1:numel(Uej)
    s(k) = sym(sprintf('s%d', k));
end
f=blEQ2(theta,Hi,s,cf_2,xj_d,dc_Hstar,Hstari,Naw);
Jf=jacobian(f,s);
df=subs(Jf,s,Uej);
end