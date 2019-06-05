function H = Func_NB_Hptot(v1,v2,v3,v4,v5,v6)
% A function that calculates the host population (before mixing)
% More advanced version considering genotypes

Ho = v1;        % Host population
Hotot = v2;     % Total host population
Po = v3;        % Relevant parasite population (aA -> bB etc.)
lam = v4;       % Base host growth factor
K = v5;         % Host max carrying capacity
a = v6;         % parasite searching efficiency

%if(Ho < 1e-14)     % If pop too low it goes extinct
%    H = 0;
%else
    H = Ho*lam*exp(-Hotot*log(lam)/K)*exp(-a*Po);
%end
end