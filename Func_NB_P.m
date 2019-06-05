function P = Func_NB_P(v1,v2,v3,v4)
% A function that calculates the parasite population (before mixing)

Ho = v1;    % Total hosts of relevant genotype (bB -> aA etc.)
Po = v2;    % Total parasite population
a = v3;     % parasite searching efficiency
c = v4;     % fecundity, "parasite infection success"

%if(Po < 1e-14)     % If pop too low, it is extinct
%    P = 0;
%else
    P = c*Ho*(1 - exp(-a*Po));
%end
end