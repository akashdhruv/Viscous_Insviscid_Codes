function [df,Jf]=jacobl1(theta,Hi,Uej,cf_2,xj_d,Naw)
s = sym(zeros(size(theta)));
for k = 1:numel(theta)
    s(k) = sym(sprintf('s%d', k));
end
f=blEQ1(s,Hi,Uej,cf_2,xj_d,Naw);
Jf=jacobian(f,s);
df=subs(Jf,s,theta);
end