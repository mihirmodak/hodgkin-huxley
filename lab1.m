clear
clc
close all

%% Define variables

T = 310.15; % K
R = 8.3144598; % Joule/mol*K

F = 96485.332/1000; 
% Every time we use Faraday's constant, e.g. Nernst Equation or GHK
% Equation, we need to convert from V to mV. So, we need to multiply by
% 1000. This is the same as dividing by 1/1000. So we divide F by 1000
% since F is in the denominator of these equations

z_Na = 1; z_K = 1; z_Cl = -1;
C_Na_in = 10   ; C_Na_out = 145; %mM
C_K_in = 140; C_K_out = 5; %mM
C_Cl_in = 10; C_Cl_out = 110; %mM

%% 3.1 Calculate Membrane Potentials

Em_Na = (R*T)/(z_Na * F) .* log(C_Na_out/C_Na_in)
Em_K = (R*T)/(z_K * F) .* log(C_K_out/C_K_in)
Em_Cl = (R*T)/(z_Cl * F) .* log(C_Cl_out/C_Cl_in)

P_Na = 0.02; P_K = 1; P_Cl = 0.5;
Vm = (R*T/F) * log((P_Na*C_Na_out + P_K*C_K_out + P_Cl*C_Cl_in)/...
    (P_Na*C_Na_in + P_K*C_K_in + P_Cl*C_Cl_out) )

%% 3.2 Calculate Resting Potential

g_K_max = 36; %mS/cm2
g_Na_max = 120; %mS/cm2
n_initial = 0.337; m_initial = 0.061; h_initial = 0.552;

gk = @ (n) g_K_max .* (n^4);
gNa = @ (m,h) g_Na_max .* (m^3) .* h;

gL = 0.3; % mS/cm2
E_L = -55; % mV

Vr = (gk(n_initial)*Em_K + gNa(m_initial, h_initial)*Em_Na + ...
    gL*E_L)/(gk(n_initial) + gNa(m_initial,h_initial) + gL)


%% 3.3 Develop MATLAB Code

V_initial = -75; % mV
n_initial = 0.337; m_initial = 0.061; h_initial = 0.552;

stim_params = {[0], [0 0]}; % see section 3.4 for explanation
% stim_params = {300, [1.0 1.1 40 40.1]}; % [magnitude of stimulus(uA), 
                                          % start time, end time]
                                          
options = odeset('MaxStep', 0.01);

[t,y] = ode45(@ (t,y) hh_model(t,y, Vr, Em_K, Em_Na, E_L, g_K_max, ...
    g_Na_max, stim_params), [0 20], [V_initial; n_initial; m_initial;...
    h_initial], options);

Vm = y(:,1); n=y(:,2); m=y(:,3);h=y(:,4);

figure;
subplot(2,1,1)
plot(t,Vm);
title("Hodgkin Huxley Model at rest")
ylim([-100 100])
ylabel("V_m [mV]")

subplot(2,1,2)
plot(t, n)
hold on
plot(t, m)
plot(t,h)
ylim([0 1])
title("Gating Parameters")
xlabel("time [ms]")

legend('n', 'm', 'h')


%% 3.4 Generating an Action Potential

V_initial = -75; % mV
n_initial = 0.337; m_initial = 0.061; h_initial = 0.552;
options = odeset('MaxStep', 0.01);

% Hodgkin Huxley Model with stimulus

stim_params = {3000, [1.0 1.1]}; % [magnitude of stimulus(uA), start time,
                                % end time]
% I_stim_params is a vecotr containing 2 elements: (1) vector of the 
% magnitudes of successive stimuli, and (2) a vector containing the start 
% and stop times of successive stimuli e.g. [1.0 1.1 40.0 40.1]


[t2,y2] = ode45(@ (t2,y2) hh_model(t2,y2, Vr, Em_K, Em_Na, E_L, g_K_max,...
    g_Na_max, stim_params), [0 20], [V_initial; n_initial; m_initial;...
    h_initial], options);
Vm2 = y2(:,1); n2=y2(:,2); m2=y2(:,3);h2=y2(:,4);

figure;
subplot(2,1,1)
plot(t2,Vm2)
title(['Hodgkin Huxley Model with ',num2str(cell2mat(stim_params(1))) ,...
    ' nA stimulus'])
ylabel("V_m [mV]");
ylim([-100 100]);

subplot(2,1,2)
plot(t2, n2)
hold on
plot(t2, m2)
plot(t2,h2)
ylim([0 1])
title("Gating Parameters")
xlabel("time [ms]")

legend('n', 'm', 'h')

%% 3.4 Continued

% Hodgkin Huxley Model with stimulus

V_initial = -75; % mV
n_initial = 0.337; m_initial = 0.061; h_initial = 0.552;
options = odeset('MaxStep', 0.01);

% Below threshold
% stim_params = {78, [1.0 1.1]}; % [magnitude of stimulus(uA), start time, 
                                 % end time]

% Above threshold
stim_params = {[79], [1.0 1.1]}; % [magnitude of stimulus(uA), 
                                   % start time, end time]

[t3,y3] = ode45(@ (t3,y3) hh_model(t3,y3, Vr, Em_K, Em_Na, E_L, g_K_max,...
    g_Na_max, stim_params), [0 20], [V_initial; n_initial; m_initial; ...
    h_initial], options);
Vm3 = y3(:,1); n3=y3(:,2); m3=y3(:,3);h3=y3(:,4);

figure;
subplot(2,1,1)
plot(t3,Vm3)
title(['Hodgkin Huxley Model with ',num2str(cell2mat(stim_params(1))) ,...
    ' nA stimulus'])
ylabel("V_m [mV]"); 
ylim([-100 100]);

subplot(2,1,2)
plot(t3, n3)
hold on
plot(t3, m3)
plot(t3, h3)
ylim([0 1])
title("Gating Parameters")
xlabel("time [ms]")

legend('n', 'm', 'h')

%% 3.5 Hyperpolarized Action Potential

% Hodgkin Huxley Model with stimulus

V_initial = -75; % mV
n_initial = 0.337; m_initial = 0.061; h_initial = 0.552;
options = odeset('MaxStep', 0.01);

% Below threshold
stim_params = {-70, [1.0 1.1]}; % [magnitude of stimulus(uA), start time,
                                  % end time]

% Above threshold
% stim_params = {-72, [1.0 1.1]}; % [magnitude of stimulus(uA), start time,
                                % end time]

[t4,y4] = ode45(@ (t4,y4) hh_model(t4,y4, Vr, Em_K, Em_Na, E_L, g_K_max,...
    g_Na_max, stim_params), [0 20], [V_initial; n_initial; m_initial; ...
    h_initial], options);
Vm4 = y4(:,1); n4=y4(:,2); m4=y4(:,3);h4=y4(:,4);

figure;
subplot(2,1,1)
plot(t4,Vm4)
title(['Hodgkin Huxley Model with ',num2str(cell2mat(stim_params(1))) ,...
    ' nA stimulus'])
ylabel("V_m [mV]");
ylim([-100 100]);

subplot(2,1,2)
plot(t4, n4)
hold on
plot(t4, m4)
plot(t4, h4)
ylim([0 1])
title("Gating Parameters")
xlabel("time [ms]")

legend('n', 'm', 'h')

%% 3.6 Refractory Periods

% Hodgkin Huxley Model with 2 stimuli in close proximity

V_initial = -75; % mV
n_initial = 0.337; m_initial = 0.061; h_initial = 0.552;
options = odeset('MaxStep', 0.01);

% Below threshold
stim_params = {[480,388], [1.0 1.1 8.0 8.1]}; % [magnitude of stimulus
                                                 % (nA), start time, end 
                                                 % time]


[t5,y5] = ode45(@ (t5,y5) hh_model(t5,y5, Vr, Em_K, Em_Na, E_L, g_K_max,...
    g_Na_max, stim_params),[0 35], [V_initial; n_initial; m_initial; ...
    h_initial], options);
Vm5 = y5(:,1); n5=y5(:,2); m5=y5(:,3);h5=y5(:,4);

figure;
subplot(2,1,1)
plot(t5,Vm5)

title('Hodgkin Huxley Model with multiple stimuli')
ylabel("V_m [mV]"); 
xlim([0 35]);ylim([-100 100]);

subplot(2,1,2)
plot(t5, n5)
hold on
plot(t5, m5)
plot(t5, h5)
xlim([0 40]);ylim([0 1])
title("Gating Parameters")
xlabel("time [ms]")
legend('n', 'm', 'h')

%% 3.6 Refractory Period Plot

% Import the data from the manually created excel file
opts=detectImportOptions('thresholds.csv');
opts.VariableNamesLine = 2;
opts.Delimiter =',';
thresholds = readtable('thresholds.csv', opts, 'ReadVariableNames', true);

thresholds(1,:) = []; % remove the row containing the initial stimulus 
                      % at 1.0 - 1.1

relative_refractory = thresholds(2:end, 1:end-1)
absolute_refractory = table2array( thresholds(1,3) )

diff = table2array( relative_refractory(:,3) );
thresh = table2array( relative_refractory(:,4) );

figure;
scatter(diff, thresh)
xline(absolute_refractory, 'LineWidth', 1)
text(0, -50, 'Absolute'); 
text(15, -50, 'Relative'); 
ylabel('Current (nA)'); xlabel('Time difference (ms)')
title('Refractory Period')

%% 3.7 Experimental Analysis Part 1

% Hodgkin Huxley Model with a more negative Vr, single stimulus

V_initial = -75; % mV
n_initial = 0.337; m_initial = 0.061; h_initial = 0.552;
options = odeset('MaxStep', 0.01);

% Set the stimulus
stim_params = {300, [1.0 1.1]}; 

Vr_new = - 90;

% Compute the Vm using hh_model

[t6,y6] = ode45(@ (t6,y6) hh_model(t6,y6, Vr_new, Em_K, Em_Na, E_L, g_K_max,...
    g_Na_max, stim_params),[0 100], [V_initial; n_initial; m_initial; ...
    h_initial], options);
Vm6 = y6(:,1); n6=y6(:,2); m6=y6(:,3);h6=y6(:,4);


figure;
subplot(2,1,1)
plot(t6,Vm6)
title("Hodgkin Huxley Model")


subplot(2,1,2)
plot(t6, n6)
hold on
plot(t6, m6)
plot(t6, h6)
ylim([0 1])
title("Gating Parameters")
xlabel("time [ms]")
legend('n', 'm', 'h')

%% 3.7 Experimental Analysis Part 2

% Hodgkin Huxley Model with a more positive Vr

V_initial = -75; % mV
n_initial = 0.337; m_initial = 0.061; h_initial = 0.552;
options = odeset('MaxStep', 0.01);

% Set the stimulus
stim_params = {750, [1.0 1.1]}; 

Vr_new = - 50;

% Compute the Vm using hh_model

[t7,y7] = ode45(@ (t7,y7) hh_model(t7,y7, Vr_new, Em_K, Em_Na, E_L, g_K_max,...
    g_Na_max, stim_params),[0 40], [V_initial; n_initial; m_initial; ...
    h_initial], options);
Vm7 = y7(:,1); n7=y7(:,2); m7=y7(:,3);h7=y7(:,4);


figure;
subplot(2,1,1)
plot(t7,Vm7)
title("Hodgkin Huxley Model")


subplot(2,1,2)
plot(t7, n7)
hold on
plot(t7, m7)
plot(t7, h7)
xlim([0 40]);ylim([0 1])
title("Gating Parameters")
xlabel("time [ms]")
legend('n', 'm', 'h')

%% Creating a gif from a plot 
% with modified conditions of Experimental Analysis Part 1

% Hodgkin Huxley Model with a more negative Vr, single stimulus

V_initial = -75; % mV
n_initial = 0.337; m_initial = 0.061; h_initial = 0.552;
options = odeset('MaxStep', 0.01);

% Set the stimulus
stim_params = {300, [1.0 1.1]}; 

Vr_new = - 90;

% Compute the Vm using hh_model

% [t_part,y_part] = ode45( ...
%     @ (t_part,y_part) hh_model(t_part,y_part, Vr_new, Em_K, Em_Na, E_L, g_K_max, g_Na_max, stim_params), ...
%     [0 100], ...
%     [V_initial; n_initial; m_initial; h_initial], ...
%     options ...
% );
% Vm = y_part(:,1); n=y_part(:,2); m=y_part(:,3);h=y_part(:,4);

% % The figure as it was before
% figure;
% subplot(2,1,1)
% plot(t_part, Vm)
% title("Hodgkin Huxley Model")
% 
% subplot(2,1,2);
% plot(t_part, n, 'Color',[0, 0.4470, 0.7410])
% hold on
% plot(t_part, m, 'Color', [0.8500, 0.3250, 0.0980])
% plot(t_part, h, 'Color', [0.9290, 0.6940, 0.1250])
% ylim([0 1])
% title("Gating Parameters")
% xlabel("time [ms]")
% legend('n', 'm', 'h')

% Trying to make an animated gif
figure;
for t = 1:length(t6)

    [t_part,y_part] = ode45( ...
        @ (t_part,y_part) hh_model(t_part,y_part, Vr_new, Em_K, Em_Na, E_L, g_K_max, g_Na_max, stim_params), ...
        [0 t], ...
        [V_initial; n_initial; m_initial; h_initial], ...
        options ...
    );
    Vm = y_part(:,1); n=y_part(:,2); m=y_part(:,3);h=y_part(:,4);

    subplot(2,1,1);
    plot(t_part,Vm, 'Color', [0, 0.4470, 0.7410])
    title("Hodgkin Huxley Model")
%     xlim([0 200])
    
    subplot(2,1,2);
    plot(t_part, n, 'Color',[0, 0.4470, 0.7410])
    hold on
    plot(t_part, m, 'Color', [0.8500, 0.3250, 0.0980])
    plot(t_part, h, 'Color', [0.9290, 0.6940, 0.1250])
    hold off
%     xlim([0 200])
    ylim([0 1])
    title("Gating Parameters")
    xlabel("time [ms]")
    legend('n', 'm', 'h')
    
    % gif utilities
    set(gcf,'color','w'); % set figure background to white
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    outfile = 'hh_viz.gif';
    
    % On the first loop, create the file. In subsequent loops, append.
    if t==1
        imwrite(imind,cm,outfile,'gif','DelayTime',0,'loopcount',inf);
    else
        imwrite(imind,cm,outfile,'gif','DelayTime',0,'writemode','append');
    end

end
