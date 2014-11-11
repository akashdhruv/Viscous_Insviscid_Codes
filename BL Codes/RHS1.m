function f=RHS1(theta,H,Ve,nu,dVeds)
f=0.5*cf(theta,H,Ve,nu) - theta/Ve*(2+H)*dVeds;
end