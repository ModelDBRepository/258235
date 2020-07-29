% A mathematical model for a bursting neuron with Nav1.6-type Na+ currents
% The model reproduces bursting activity in brainstem proprioceptive
% Mesencephalic V neurons.
% Author: Sharmila Venugopal, Ph.D.
% Date: 10/02/2018
% Email: vsharmila@g.ucla.edu
% Department of Integrative Biology & Physiology, UCLA
% If you benefit from this model/code, please cite us.
%% 

clc;
clear all;

% % Initial values for membrane potential and channel gating variables
V=-60;   % Initial membrane voltage
hp=0.3; % Initial inactivation gate for persistent sodium
n=1/(1+exp(-(V-(-43))/3.9));  % Initial activation gate for fast sodium
h=(1/(1+exp((V+55)/7.1)));   % Initial inactivation gate for fast sodium
br=0.9;  % Initial activation gate for resurgent sodium
hr=0;  % Initial inactivation gate for resurgent sodium

y0=[V; hp; n; h; br; hr];

% Simulation time
tmax=5000;
tstep=0.01;
tspan=0:tstep:tmax;

[time,sol]=ode45(@BurstingNeuron,tspan,y0);     % Solve using ode45

hp=sol(:,2);
h=sol(:,4);
br=sol(:,5);
hr=sol(:,6);
V = sol(:,1);

%%% 
figure
plot(time,V,'k');
xlabel('Time (ms)');
ylabel('Membrane Voltage, V (mV)');
set(gcf, 'Position', [100, 100, 600, 500]);

%% 

% BurstingNeuron(t,y): This function models a bursting neuron with a novel
% formulation for resurgent sodium current.
% Burst patterns resemble observed in vitro bursting in Mesencephalic V
% neurons.
function dydt = BurstingNeuron(t,y)

    % Setting input array to corresponding model variables
    V=y(1);
    hp=y(2);
    n=y(3);
    h=y(4);
    br=y(5);
    hr=y(6);

    % Membrane Constants
    C=1;            % Membrane capacitance
    % Leak current parameters (Ileak)
    Eleak=-62;      % Leak reversal potential
    gleak=2;        % Leak conductance

    % Slower potassium current parameters (IK)
    EK=-80;         % Potassium reversal potential
    tauV=4;         
    KV12=-43;
    Kk=3.9;
    gK=7;           % K+ conductance

    % Parameters of Sodium currents
    NaV12=-50;
    Nak=6.4;
    tauphVbs=100;
    tauphVst=10000;
    tauh=2;         

    Ah=3;           % Default value provided. Values ranged from 1 to 3 for various simulations presented to fit experimental data
    Ak=9;           % Default value provided. Values ranged from 9 to 11 for various simulations presented to fit experimental data
    Bh=2.5;         % Default value provided. Values ranged from 1 to 2.5 for various simulations presented to fit experimental data
    Ab=0.05;        % Default value provided. Values ranged from 0.05 to 0.1 for various simulations presented to fit experimental data 
    kb=0.3;         % Default value provided. Values ranged from 0.3 to 1 for various simulations presented to fit experimental data 
%     
    % % Maximal conductances
    gNat=22;        % Transient Na+ conductance
    gNar=5;         % Resurgent Na+ conductance. Values ranged from 3 to 12 nS to reproduce dynamic-clamp data
    gNap=0.45;       % Persistent Na+ conductance Values ranged from 0.45 to 0.9 nS to reproduce dynamic-clamp data
    ENa=55;         % Sodium reversal potential
    
    % Injected current
    Iapp = 4;

    % Voltage-dependent channel gating functions
    function mt = mtinf(arg1) 
        mt = 1/(1+exp(-(arg1+35)/4.3));
    end

    function h = hinf(arg1) 
        h = 1/(1+exp((arg1+55)/7.1));
    end

    function mp = mpinf(arg1) 
        mp = 1/(1+exp(-(arg1-NaV12)/Nak));
    end

    function hp = hpinf(arg1) 
        hp = 1/(1+exp((arg1+52)/14));   
    end

    function hhr = hrinf(arg1) 
        hhr = 1/(1+exp((arg1+40)/20));        
    end

    function alphah = alphahr(arg1)
        alphah = (Ah*(1/(1+exp(-(arg1+45)/Ak))));
    end

    function betah = betahr(arg1)
        betah = (0.5*(1/(1+exp(-(arg1+40)/15))));
    end

    function bbr = brinf(arg1) 
        bbr = 1/(1+exp((arg1+40)/12));        
    end

    function beta = betabr(arg1)
       beta = (2*(1/(1+exp(-(arg1-40)/8))));
    end

    function nn = ninf(arg1) 
        nn = 1/(1+exp(-(arg1-KV12)/Kk));
    end

    function taup = tauphV(arg1) 
        taup = (tauphVbs+tauphVst/(1+exp((arg1+60)/10)));        
    end

    ivK = gK*n*(V-EK);      % TEA sensitive Kv1.2 type K+ current
    ivL = gleak*(V-Eleak);  % Leak current
    ivNa = ((gNap*mpinf(V)*hp) + (gNat*mtinf(V)*h) + (gNar*(1-br)^3*hr^5))*(V-ENa); % Nav1.6 type Na+ current with transient, resurgent and persistent components

    % Membrane potential differential equation
    dydt = [((Iapp)-ivL-(ivNa)-(ivK))/C; (hpinf(V)-hp)/tauphV(V); (ninf(V)-n)/tauV; (hinf(V)-h)/tauh; (Ab*(1-br)*brinf(V))-kb*(betabr(V)*br); (alphahr(V)*hrinf(V))-Bh*betahr(V)*hr];
end