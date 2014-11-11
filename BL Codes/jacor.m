function df=jacor(B,sig,RHS)
s = sym(zeros(size(RHS)));
for k = 1:numel(RHS)
    s(k) = sym(sprintf('s%d', k));
end
f=rhsEQ(B,sig,s);
Jf=jacobian(f,s);
df=subs(Jf,s,RHS);
end