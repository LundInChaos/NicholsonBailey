function M = Func_NB_Mix(v1,v2,v3)
% A function that is used to calculate how the genotypes of host/parasite
% gets distributed after mixing. When calculating for hosts only host
% populations are used, for parasites only parasite populations.

HPAA = v1;      % Host/Parasite of genotype AA/BB
HPAa = v2;      % Host/Parasite of genotype aA/bB
HPtot = v3;     % Total population of host/Parasite

if(HPtot < 1e-14)   % Avoid divission by zero
    HPtot = 1e-14;
end

M = (HPAA + HPAa/2)/HPtot;

end