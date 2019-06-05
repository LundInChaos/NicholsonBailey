function M = Func_NB_Mix(v1,v2,v3)

HPAA = v1;
HPAa = v2;
HPtot = v3;

if(HPtot < 1e-14)
    HPtot = 1e-14;
end

M = (HPAA + HPAa/2)/HPtot;

end