function f=fH1(H)
if(H<=1.6)
    f=3.3 + 0.8234*sign(H-1.1)*((abs(H-1.1))^(-1.287));
else
    f=3.3 + 1.5501*sign(H-0.6678)*((abs((H-0.6678)))^(-3.064));
end
end