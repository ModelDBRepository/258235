% A mathematical model for Nav1.6-type Na+ currents with novel resurgent
% current formulation
% The model simulates a voltage-clamp protocol. The total Ina = Inat + Inap
% + Inar.
% Author: Sharmila Venugopal, Ph.D.
% Date: 10/02/2018
% Email: vsharmila@g.ucla.edu
% Department of Integrative Biology & Physiology, UCLA
% If you benefit from this model/code, please cite us.
%% 
% Voltage clamp simulations of transient, resurgent and persistent sodium
% currents
clear all;
clc;

% Global model parameters
global ENa gNat gNar gNap;

% Vclamp protocol
vhold=-90;  % holding potential
vdep=30;    % unblocking potential
vtest=-70:10:-10;  % Test potentials to reveal resurgent current

dhold=50;   % duration of 1st pulse at holding potential
ddep=5;     % duration of 2nd pulse, at test potential
dtest=100;  % duration of 3rd pulse at unblocking potential

dt=0.1;     % integration interval

setParams();
y0 = setInit(vhold);

% Simulating Vclamp experiment
for i=1:length(vtest)
    [t1,y1]=ode23s(@rates1,0:dt:dhold,y0,[],vhold);     % Solve using ode23s
    [t2,y2]=ode23s(@rates1,0:dt:ddep,y1(end,:),[],vdep);
    [t3,y3]=ode23s(@rates1,0:dt:dtest,y2(end,:),[],vtest(i));
    [t4,y4]=ode23s(@rates1,0:dt:dhold,y3(end,:),[],vhold);
    
    V(:,i)=[(vhold.*ones(length(t1),1))' (vdep.*ones(length(t2),1))' (vtest(i).*ones(length(t3),1))' (vhold.*ones(length(t4),1))'];
    Inat(:,i)=gNat*[((mtinf(vhold)*y1(:,2))*(vhold-ENa))' ((mtinf(vdep)*y2(:,2))*(vdep-ENa))' ((mtinf(vtest(i))*y3(:,2))*(vtest(i)-ENa))' ((mtinf(vhold)*y4(:,2))*(vhold-ENa))'];
    Inar(:,i)=gNar*[((1-y1(:,3)).^3.*(y1(:,4).^5)*(vhold-ENa))' ((1-y2(:,3)).^3.*(y2(:,4).^5)*(vdep-ENa))' ((1-y3(:,3)).^3.*(y3(:,4).^5)*(vtest(i)-ENa))' ((1-y4(:,3)).^3.*(y4(:,4).^5)*(vhold-ENa))'];
    Inap(:,i)=gNap*[((mpinf(vhold)*y1(:,1))*(vhold-ENa))' ((mpinf(vdep)*y2(:,1))*(vdep-ENa))' ((mpinf(vtest(i))*y3(:,1))*(vtest(i)-ENa))' ((mpinf(vhold)*y4(:,1))*(vhold-ENa))'];
    Ina(:,i)=Inat(:,i) + Inar(:,i) + Inap(:,i);
end

t=[t1' (t2+dhold)' (t3+dhold+ddep)' (t4+dhold+ddep+dtest)'];
makeplots(t,V,Ina,Inat,Inap,Inar);
%% 
% setParams(): Sets values for global parameters
function setParams()
    % Global model parameters
    global ENa gNat gNar gNap Ab tauh;
    global NaV12 Nak; 
    global tauphVbs tauphVst;
    global Ah Ak Bh kb;

    % Parameters of Sodium currents
    NaV12=-50;
    Nak=6.4;
    tauphVbs=100;
    tauphVst=10000;
    tauh=1.5;         % Range: 1.5 to 2

    Ah=3; 
    Ak=9;
    Bh=2.5;         % Default value 0.8. Range tested: 0.5 to 1
    Ab=0.05;        % Default value 0.08. Range tested: 0.08 to 0.1
    kb=0.3;           % Default value 1. Ramge tested: 0.8 to 1.2
%   
    % % Maximal conductances
    gNat=22;        % Transient Na+ conductance
    gNar=3;         % Resurgent Na+ conductance: values ranged from 5+/- 2nS for Vclamp simulations
    gNap=0.45;       % Persistent Na+ conductance
    ENa=55;         % Sodium reversal potential

 end
%% 
% Voltage-dependent functions of sodium currents
function mt = mtinf(arg1) 
    mt = 1/(1+exp(-(arg1+35)/4.3));
end

function h = hinf(arg1) 
  h = 1/(1+exp((arg1+55)/7.1));
end

function mp = mpinf(arg1)
    global NaV12 Nak; 
    mp = 1/(1+exp(-(arg1-NaV12)/Nak));
end

function hp = hpinf(arg1) 
    hp = 1/(1+exp((arg1+52)/14));   
end

function hhr = hrinf(arg1) 
    hhr = 1/(1+exp((arg1+40)/20));        
end

function alphah = alphahr(arg1)
    global Ah Ak;
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

function taup = tauphV(arg1) 
    global tauphVbs tauphVst;
    taup = (tauphVbs+tauphVst/(1+exp((arg1+60)/10)));        
end
%% 
% setInit(): sets and returns an initial value vector for the gating
% variables
function y0 = setInit(V)
    hp=1/(1+exp((V+52)/14)); % Initial inactivation gate for persistent sodium
    h=(1/(1+exp((V+55)/7.1))); % Initial inactivation gate for fast sodium
    br=1/(1+exp((V+40)/12));  % Initial activation gate for resurgent sodium
    hr=1/(1+exp((V+40)/20));  % Initial inactivation gate for resurgent sodium

    y0=[hp; h; br; hr]; % Initial conditions vector
end
%% 
% rates1(t,y,v) 
function dydt = rates1(t,y,V)
    global  Ab tauh kb Bh;
    % Setting input array to corresponding model variables
    hp=y(1);
    h=y(2);
    br=y(3);
    hr=y(4);
    
    % Model Nav1.6 type Na+ current gating equations
    dydt = [(hpinf(V)-hp)/tauphV(V); (hinf(V)-h)/tauh; (Ab*(1-br)*brinf(V))-kb*(betabr(V)*br); (alphahr(V)*hrinf(V))-Bh*betahr(V)*hr];
end
%% 
% makeplots(): Graphs the result plots
function makeplots(t,V,Ina,Inat,Inap,Inar)
    XMIN = 30;
    XMAX = 185;
    YMIN = -600;
    YMAX = 50;
    
    figure
    subplot(2,1,1);
    plot(t,V,'LineWidth',1.5);
    axis([XMIN XMAX -100 YMAX]);
    xlabel('Time (ms)');
    ylabel('V_{cmd}');
    
    subplot(2,1,2);
    plot(t,Ina,'LineWidth',1.5);
    axis([XMIN XMAX YMIN YMAX]);
    xlabel('Time (ms)');
    ylabel('I_{Na}');
    set(gcf, 'Position', [100, 100, 300, 700]);
    
    figure
    subplot(2,2,1);
    plot(t,Ina(:,4),'k','LineWidth',1.5);
    axis([XMIN XMAX YMIN YMAX]);
    xlabel('Time (ms)');
    ylabel('I_{Na}');

    subplot(2,2,2);
    plot(t,Inat(:,4),'LineWidth',1.5);
    axis([XMIN XMAX YMIN YMAX]);
    xlabel('Time (ms)');
    ylabel('I_{transient}');
    
    subplot(2,2,3);
    plot(t,Inap(:,4),'m','LineWidth',1.5);
    axis([XMIN XMAX -100 YMAX]);
    xlabel('Time (ms)');
    ylabel('I_{persistent}');
    
    subplot(2,2,4);
    plot(t,Inar(:,4),'r','LineWidth',1.5);
    axis([XMIN XMAX -300 YMAX]);
    xlabel('Time (ms)');
    ylabel('I_{resurgent}');
    set(gcf, 'Position', [600, 100, 500, 700]);
    
end