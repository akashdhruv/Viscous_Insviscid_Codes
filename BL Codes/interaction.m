function f=interaction(var)
load('constants.mat','A','C','a','c','RHS','Na','nu','Naw','xj_d','Uinf','Vinf','t_com')
mu=var(:,1);
mdef=var(:,2);
theta=var(:,3);


f(:,1)=A*mdef+C*mu-RHS;

Uei=([Uinf Vinf]*t_com)'+a*mdef+c*mu;

dstar=mdef./Uei;

H=dstar./theta;

Hstar=calc_hstar(H,Na);

Re_THT = (Uei.*theta)/nu;
cf_2 = calc_cf(H,Na);
cf_2 = (1./Re_THT).*cf_2;

dcHS = calc_DCHstar(H,Na);
dc_Hstar = (1./Re_THT).*dcHS;


for i=1:Naw
    if i==1
    f(i,2)=((theta(i))/0.5*(theta(i)))+((2+0.5*(H(i)))*((Uei(i))/0.5*(Uei(i))))-...
            (0.5*(cf_2(i))*(xj_d(i))/(0.5*(theta(i))));
    else
    f(i,2)=((theta(i)-theta(i-1))/0.5*(theta(i)+theta(i-1)))+((2+0.5*(H(i)-H(i-1)))*((Uei(i)-Uei(i-1))/0.5*(Uei(i)+Uei(i-1))))-...
            (0.5*(cf_2(i)+cf_2(i-1))*(xj_d(i)-xj_d(i-1))/(0.5*(theta(i)+theta(i-1))));
    end
end

for i=1:Naw
    if i==1
    f(i,3)=((Hstar(i))/0.5*(Hstar(i)))+((1-0.5*(H(i)))*((Uei(i))/0.5*(Uei(i))))+...
            ((0.5*(cf_2(i))-0.5*(dc_Hstar(i)))*(xj_d(i))/(0.5*(theta(i))));
    else
    f(i,3)=((Hstar(i)-Hstar(i-1))/0.5*(Hstar(i)+Hstar(i-1)))+((1-0.5*(H(i)-H(i-1)))*((Uei(i)-Uei(i-1))/0.5*(Uei(i)+Uei(i-1))))+...
            ((0.5*(cf_2(i)+cf_2(i-1))-0.5*(dc_Hstar(i)+dc_Hstar(i-1)))*(xj_d(i)-xj_d(i-1))/(0.5*(theta(i)+theta(i-1))));
    end
end

end