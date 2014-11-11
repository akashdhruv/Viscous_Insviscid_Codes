function f=panelEQ(A,var,Bt,sig,RHS,xj_d,Hi,cf_2,Naw,dc_Hstar,Hstari)
mu=var(1:Naw);
Uej=var(Naw+1:Naw+Naw+1);
theta=var(2*(Naw+1):2*(Naw+1)+Naw);

q=A*mu+Bt*sig-RHS;

for i=1:Naw
    r1(i,1)=((theta(i+1)-theta(i))/0.5*(theta(i+1)+theta(i)))+((2+0.5*(Hi(i+1)-Hi(i)))*((Uej(i+1)-Uej(i))/0.5*(Uej(i+1)+Uej(i))))-...
            (0.5*(cf_2(i+1)+cf_2(i))*(xj_d(i+1)-xj_d(i))/(0.5*(theta(i+1)+theta(i))));
end

for i=1:Naw
    r2(i,1)=((Hstari(i+1)-Hstari(i))/0.5*(Hstari(i+1)+Hstari(i)))+((1-0.5*(Hi(i+1)-Hi(i)))*((Uej(i+1)-Uej(i))/0.5*(Uej(i+1)+Uej(i))))+...
            ((0.5*(cf_2(i+1)+cf_2(i))-0.5*(dc_Hstar(i+1)+dc_Hstar(i)))*(xj_d(i+1)-xj_d(i))/(0.5*(theta(i+1)+theta(i))));
end
f=[q;r1;r2];
end