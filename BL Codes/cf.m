function c=cf(theta,H,Ve,nu)
c=0.246*10^(-0.678*H)*sign(theta*Ve/nu)*(abs((theta*Ve/nu))^(-0.268));
end