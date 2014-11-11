function f=vi_int(x)
load('aero_variables.mat','Uinf','Vinf','nu');
load('wing_geom.mat','xi_d','zi_d','si_d','xj_d','zj_d','sj_d','xi_s','zi_s','si_s','xj_s','zj_s','sj_s','n_com','t_com','Na','Nw','Naw');
load('influences.mat','B','C','sigma');

mu=(x(1:end-1,1)+x(2:end,1))/2;
mdef=(x(:,2));
theta=(x(:,3));

    for i=1:Naw
    sig(i,1)=(mdef(i+1)-mdef(i))/(xj_s(i)-xj_s(i+1));
    end
       
    for i = 1:Naw
        Vt=0;
        Vn=0;
    for j=1:Naw
        [ud,wd]=doublet_2Dc(mu(j),xi_d(i),zi_d(i),xj_d(j),zj_d(j),xj_d(j+1),zj_d(j+1));
        [us,ws]=source_2Dc(sig(j),xi_d(i),zi_d(i),xj_d(j),zj_d(j),xj_d(j+1),zj_d(j+1));
        Vt=Vt-ud-us;
        Vn=Vn-wd-ws;
    end
    Vt=Vt+Uinf;
    Vn=Vn+Vinf;
    Uei(i,1)=[Vt Vn]*t_com(:,i);
    end
    
    dif = sqrt((xi_d(1)-xj_d(1))^2+(zi_d(1)-zj_d(1))^2);
    Uej = ppval(pchip(si_d+dif,Uei),sj_d);
    Uej=Uej';
   
    
    dstar = mdef./Uej;
    Hi = dstar./theta;
    Hstari = calc_hstar(Hi,Na);

    Re_THT = (Uej.*theta)/nu;
    cf_2 = calc_cf(Hi,Na);
    cf_2 = (1./Re_THT).*cf_2;

    dcHS = calc_DCHstar(Hi,Na);
    dc_Hstar = (1./Re_THT).*dcHS;


for i=1:Naw
    f(i,1)=((theta(i+1)-theta(i))/0.5*(theta(i+1)+theta(i)))+((2+0.5*(Hi(i+1)-Hi(i)))*((Uej(i+1)-Uej(i))/0.5*(Uej(i+1)+Uej(i))))-...
            (0.5*(cf_2(i+1)+cf_2(i))*(xj_d(i+1)-xj_d(i))/(0.5*(theta(i+1)+theta(i))));
end

for i=1:Naw
    f(i,2)=((Hstari(i+1)-Hstari(i))/0.5*(Hstari(i+1)+Hstari(i)))+((1-0.5*(Hi(i+1)-Hi(i)))*((Uej(i+1)-Uej(i))/0.5*(Uej(i+1)+Uej(i))))+...
            ((0.5*(cf_2(i+1)+cf_2(i))-0.5*(dc_Hstar(i+1)+dc_Hstar(i)))*(xj_d(i+1)-xj_d(i))/(0.5*(theta(i+1)+theta(i))));
end

f(:,3)=C*mu+B*sig+B*sigma;
end