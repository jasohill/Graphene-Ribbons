% get_subbands_3NN() computes the subbands of a GNR using
% the 3NN tight binding model (TBM) with optical parameters
% and the first order perturbation theory result for interactions;
% also computes effective mass for charge calculations.
%
% inputs:
% N should be preset in define_design.m
%
% Author: Jason E. Hill, Ph.D., Post-doc research fellow at Texas Tech
% Updated: 21 MAY 2017
% CREDIT: based on previous version known simply as get_subbands()
% Sources:
% (1) Paolo Michetti and Giuseppe Iannaccone, "Analytical Model of One-Dimensional Carbon-Based
% Schottky-Barrier Transistors," IEEE TRANSACTIONS ON ELECTRON DEVICES, VOL. 57, NO. 7, JULY 2010
% (2) Jason Edward Hill, “One-dimensional Electron Systems on Graphene Edges,” 
% Ph.D. Thesis, University of Texas Press, (December 22, 2007).
% (3) S. Reich, J. Maultzsch, C. Thomsen, and P. Ordejon, "Tight-binding description of graphene"
% PHYSICAL REVIEW B 66, 035412 ~(2002).


%% Constant definitions
%-------------------------

% physical constants
q_e = 1.60217646e-19;  % fundamental electron charge magnitude [C]
epsilon_0 = 1/(4*pi*10e-6*physconst('lightspeed')^2); 
                       % vacuum permitivity [F/cm]
epsilon_ox = 3.45e-13;
h = 4.13566766225e-15; % Planck's constant [eV*s]
h_bar = h/(2*pi);      % Reduced Planck's constant [eV*s]
k_B = physconst('Boltzmann')/q_e; % Boltzmann's constant [eV/K]
T = 300; % room temperature [K] operation (~295 + ~5 due to electonic heating)
kT = k_B*T;
beta = 1/(kT); % represents beta in 1/ev

% graphene crystal parameters:
a = sqrt(3)*1.42e-10; % [m]

% 3NN TBM parameters (optimized to optical band):
t = 2.7;           % [eV] (NN TBM hopping parameter)
energy_2p = -2.03; % [eV]
gamma_0   = -2.79; % [eV]
gamma_1   = -0.68; % [eV]
gamma_2   = -0.30; % [eV]
s_0       = 0.300; 
s_1       = 0.046; 
s_2       = 0.039; 

% first order perturbation theory correction factor:
d = 0.12; 

%% Perform calculations
alpha = 1:N; 

A = cos(pi*alpha/(N+1)); % A sub alphas
B = cos(2*pi*alpha/(N+1)); 
C = cos(3*pi*alpha/(N+1)); 

f0 = 1 + 4.*A + 4.*A.^2;
f2 = 1 + 4.*B + 4.*B.^2;
u0 = 4.*A + 2.*B;
g0 = 2*u0 + 4.*C + 2;

E_0 = (energy_2p + gamma_1*u0).*(1 + s_1*u0);
E_1 = 2*s_0*gamma_0*f0 + (s_0*gamma_2 - s_2*gamma_0).*g0 + 2*s_2.*gamma_2.*f2;
E_2 = (energy_2p + gamma_1*u0).^2 - (gamma_0^2)*f0 - gamma_0*gamma_2*g0 - (gamma_2^2).*f2;
E_3 = (1 + s_1*u0).^2 - (s_0^2).*f0 - s_0.*s_2.*g0 - (s_2^2).*f2;

FOPT_correction = (2.*(A >-0.5)-1).*4.*d.*t./(N+1).*(sin(alpha.*pi./(N+1))).^2;

%E_zero = (0+t*sqrt(f0))./(1+0.0);
E_zero = (2*E_0-E_1+sqrt((2*E_0-E_1).^2-4*E_2.*E_3))./(2*E_3);


E_alpha = E_zero + FOPT_correction;

E_final = sort(E_alpha); 

% Finding the old index for first n_sbbd subbands

sbbd = zeros(1,n_sbbd);
mass = zeros(1,n_sbbd);
for s=1:n_sbbd
    n = 0;
    sbbd(s) = E_final(s);
    
    for i = 1:N
        if(E_alpha(i) == sbbd(s))
            n = i;
        end
    end
    mass(s) = q_e*abs(2/3*h_bar^2*sbbd(s)/(a^2*t^2*A(n))); % effective mass [kg]
end

%clear f0 f2 u0 g0 a d E_0 E_1 E_2 E_3 s_1 s_2 s_3;
%clear A E_final E_alpha E_zero alpha FOPT_correction i n;
