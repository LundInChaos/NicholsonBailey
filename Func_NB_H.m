function H = Func_NB_H(v1,v2,v3,v4,v5,v6)

Ho = v1;
Hotot = v2;
Po = v3;
lam = v4;
K = v5;
a = v6;

if(Ho < 1e-14)
    H = 0;
else
    H = Ho*exp(log(lam)*(1 - Ho/K) - a*Po);
end
end