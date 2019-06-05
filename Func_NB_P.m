function P = Func_NB_P(v1,v2,v3,v4)

Ho = v1;
Po = v2;
a = v3;
c = v4;

%if(Po < 1e-14)
%    P = 0;
%else
    P = c*Ho*(1 - exp(-a*Po));
%end
end