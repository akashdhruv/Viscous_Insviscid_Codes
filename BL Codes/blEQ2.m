function r=blEQ2(theta,Hi,Uej,cf_2,xj_d,dc_Hstar,Hstari,Naw)
%     r=(delHS./avgHS)+((1-avgH).*(delUEI./avgUEI))+((avgCF2-avgDCH).*(delX./avgTHT));
for i=1:Naw
    r(i,1)=((Hstari(i+1)-Hstari(i))/0.5*(Hstari(i+1)+Hstari(i)))+((1-0.5*(Hi(i+1)-Hi(i)))*((Uej(i+1)-Uej(i))/0.5*(Uej(i+1)+Uej(i))))+...
            ((0.5*(cf_2(i+1)+cf_2(i))-0.5*(dc_Hstar(i+1)+dc_Hstar(i)))*(xj_d(i+1)-xj_d(i))/(0.5*(theta(i+1)+theta(i))));
end
end