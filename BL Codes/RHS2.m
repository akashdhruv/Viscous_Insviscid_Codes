function f=RHS2(H1,theta,H,Ve,nu,dVeds)
f=-(H1/theta)*RHS1(theta,H,Ve,nu,dVeds) - (H1/Ve)*dVeds + F(H1)/theta;
end