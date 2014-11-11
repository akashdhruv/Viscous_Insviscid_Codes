function r=blEQ1(theta,Hi,Uej,cf_2,xj_d,Naw)
    %r=(delTHT./avgTHT)+((avgH+2).*(delUEI./avgUEI))-(avgCF2.*(delX./avgTHT));
for i=1:Naw
    r(i,1)=((theta(i+1)-theta(i))/0.5*(theta(i+1)+theta(i)))+((2+0.5*(Hi(i+1)-Hi(i)))*((Uej(i+1)-Uej(i))/0.5*(Uej(i+1)+Uej(i))))-...
            (0.5*(cf_2(i+1)+cf_2(i))*(xj_d(i+1)-xj_d(i))/(0.5*(theta(i+1)+theta(i))));
end
end