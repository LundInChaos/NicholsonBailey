% A script to iterate populations of hosts and parasites using the Nicholson-
% -Bailey model. It Includes three populations (genotypes) of hosts 
% (aa, aA and AA) and three populations of parasites 
% (bb, bB and BB) where  bB only infects aA etc. 
% The mixing facors then determine how much genetic mixing there is in host
% and parasite populations, normally 1 (complete mixing) for hosts and 
% somewhere in between 0 and 1 for parasites.
% The functions used are in the file "HostPar_funcs.py" and come from an artice
% by Flatt, J. theor. Biol. (2001) 212, 345}354 (doi:10.1006/jtbi.2001.2380).

K = 10;     % Host carrying capacity
a = 0.45;   % paracite searching efficiency
c = 1;      % fecundity, "parasite infection success"

% Lambda is the base host growth factor. Here we let it vary between 1 and 50
% like in the artice by Flatt et al.
% Range of lambdas
L1 = 1;
L2 = 50;
Ls = 1/2;

% Mixing factors
SFH = 1;
SFP = [0.01 0.5 1];

% Number time steps
N = 1000;

H = zeros(3,N);     % Hosts
P = zeros(3,N);     % Parasites

% Starting densities
H(1,1) = 4;        % Pop AA
H(2,1) = 3;        % Pop Aa
H(3,1) = 3;        % Pop aa

P(1,1) = 0.3;         % Pop BB
P(2,1) = 0.36;        % Pop Bb
P(3,1) = 0.34;        % Pop bb

% Total pop
Htot = zeros(1,N);     % Hosts
Ptot = zeros(1,N);     % Parasites

% Starting densities
Htot(1) = H(1,1) + H(2,1) + H(3,1);        % Pop H
Ptot(1) = P(1,1) + P(2,1) + P(3,1);        % Pop P

% The mixing variable is calculated with the population densities:
% Initial mixing variable:
Hp = Func_NB_Mix(H(1,1),H(2,1),Htot(1));
Pq = Func_NB_Mix(P(1,1),P(2,1),Ptot(1));

% The following vectors are used to store the last population densities of 
% each population. They are used to plot the bifurcation diagrams:
% End values
EN = 40;
E = zeros(3,EN,(L2-1)/Ls+1);
F = zeros(3,EN,(L2-1)/Ls+1);

% These vectors store the populations for a certain lambda (Thelambda)
% It is usefull for plotting and analysis purposes.
Thelambda = 15;

Hti1 = zeros(3,N);     % Hosts
Hti2 = zeros(3,N);     % Hosts
Hti3 = zeros(3,N);     % Hosts

Pti1 = zeros(3,N);     % Parasites
Pti2 = zeros(3,N);     % Parasites
Pti3 = zeros(3,N);     % Parasites

Hti1(:,1) = H(1,1);
Hti2(:,1) = H(2,1);
Hti3(:,1) = H(3,1);

Pti1(:,1) = P(1,1);
Pti2(:,1) = P(2,1);
Pti3(:,1) = P(3,1);

% Just to have the initial conditions clearly written:

disp(K)
disp(a)
disp(c)
disp(Htot(1))
disp(Ptot(1))
disp(Hp)
disp(Pq)
disp(SFH)
disp(SFP)

% For each parasite mixing factor
for s = 1:3
    % For each lambda
for l = 1:(L2-1)/Ls+1
    lambda = 1 + ((l-1)*Ls);   % Chose a new lambda
    
    % Calculate H and P:
    for n = 2:N
        for nn = 1:3
            % First calculation of populations:
            H(nn,n) = Func_NB_Hptot(H(nn,n-1),Htot(n-1),P(nn,n-1),lambda,K,a);
            P(nn,n) = Func_NB_P(H(nn,n-1),P(nn,n-1),a,c);
        end
        % Total pops:
        Htot(n) = H(1,n) + H(2,n) + H(3,n);
        Ptot(n) = P(1,n) + P(2,n) + P(3,n);
        
        % Mixing factors:
        Hp = Func_NB_Mix(H(1,n),H(2,n),Htot(n));
        Pq = Func_NB_Mix(P(1,n),P(2,n),Ptot(n));
        
        % New pops:
        H(1,n) = (1-SFH)*H(1,n) +    SFH*    Hp^2        *Htot(n);
        H(2,n) = (1-SFH)*H(2,n) +    SFH*    2*Hp*(1-Hp) *Htot(n);
        H(3,n) = (1-SFH)*H(3,n) +    SFH*    (1-Hp)^2    *Htot(n);
        
        P(1,n) = (1-SFP(s))*P(1,n) + SFP(s)* Pq^2        *Ptot(n);
        P(2,n) = (1-SFP(s))*P(2,n) + SFP(s)* 2*Pq*(1-Pq) *Ptot(n);
        P(3,n) = (1-SFP(s))*P(3,n) + SFP(s)* (1-Pq)^2    *Ptot(n);
        
        % Store the populations of a certain lambda
        if(lambda == Thelambda)
            Hti1(s,n) = H(1,n);
            Hti2(s,n) = H(2,n);
            Hti3(s,n) = H(3,n);
            Pti1(s,n) = P(1,n);
            Pti2(s,n) = P(2,n);
            Pti3(s,n) = P(3,n);
        end
                
        % New total pops:
        Htot(n) = H(1,n) + H(2,n) + H(3,n);
        Ptot(n) = P(1,n) + P(2,n) + P(3,n);
    end    
    for g = 1:EN
        % Store the last value:
        E(s,g,l) = Htot(N-(g-1));
        F(s,g,l) = Ptot(N-(g-1));
    end
end
end

% More plots than you want:

figure(1)
plot(Xvals, Hti1(1,:), 'r');
hold on
plot(Xvals, Hti2(1,:), 'g');
hold on
plot(Xvals, Hti3(1,:), 'b');
ylabel('Pop. dens.')
xlabel('time')
title('Populations, lambda = 15, s = 0')
legend('Host1','Host2','Host3')

figure(2)
plot(Xvals, Hti1(2,:), 'r');
hold on
plot(Xvals, Hti2(2,:), 'g');
hold on
plot(Xvals, Hti3(2,:), 'b');
ylabel('Pop. dens.')
xlabel('time')
title('Populations, lambda = 15, s = 0.5')
legend('Host1','Host2','Host3')

figure(3)
plot(Xvals, Hti1(3,:), 'r');
hold on
plot(Xvals, Hti2(3,:), 'g');
hold on
plot(Xvals, Hti3(3,:), 'b');
ylabel('Pop. dens.')
xlabel('time')
title('Populations, lambda = 15, s = 1')
legend('Host1','Host2','Host3')

figure(4)
plot(Xvals, Pti1(1,:), 'r');
hold on
plot(Xvals, Pti2(1,:), 'g');
hold on
plot(Xvals, Pti3(1,:), 'b');
ylabel('Pop. dens.')
xlabel('time')
title('Populations, lambda = 15, s = 0')
legend('Parasite1','Parasite2','Parasite3')

figure(5)
plot(Xvals, Pti1(2,:), 'r');
hold on
plot(Xvals, Pti2(2,:), 'g');
hold on
plot(Xvals, Pti3(2,:), 'b');
ylabel('Pop. dens.')
xlabel('time')
title('Populations, lambda = 15, s = 0.5')
legend('Parasite1','Parasite2','Parasite3')

figure(6)
plot(Xvals, Pti1(3,:), 'r');
hold on
plot(Xvals, Pti2(3,:), 'g');
hold on
plot(Xvals, Pti3(3,:), 'b');
ylabel('Pop. dens.')
xlabel('time')
title('Populations, lambda = 15, s = 1')
legend('Parasite1','Parasite2','Parasite3')

lvals = L1:Ls:L2;

figure(7)
for g = 1:EN
    scatter(lvals, E(1,g,:), '.','k')
    hold on
    scatter(lvals, F(1,g,:), '.','r')
    hold on
end
title('Parasite mixing factor = 0.01')
ylabel('Pop. dens.')
xlabel('lambda')
legend('Host','Parasite')

figure(8)
for g = 1:EN
    scatter(lvals, E(2,g,:), '.','k')
    hold on
    scatter(lvals, F(2,g,:), '.','r')
    hold on
end
title('Parasite mixing factor = 0.5')
ylabel('Pop. dens.')
xlabel('lambda')
legend('Host','Parasite')

figure(9)
for g = 1:EN
    scatter(lvals, E(3,g,:), '.','k')
    hold on
    scatter(lvals, F(3,g,:), '.','r')
    hold on
end
title('Parasite mixing factor = 1')
ylabel('Pop. dens.')
xlabel('lambda')
legend('Host','Parasite')

