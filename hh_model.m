function hh = hh_model(t,y, Vr, Em_K, Em_Na, E_L, g_K_max, g_Na_max, I_stim_params)
    
    % parse y into its individual components
    Vm = y(1); n = y(2); m = y(3); h = y(4);
    
    % parse I_stim_params into its components
    stim_vals = cell2mat(I_stim_params(1));
    t_bounds = cell2mat(I_stim_params(2));
    if length(stim_vals) == 1
        stim_vals = ones(1,length(t_bounds)/2) * stim_vals;
    end
    
    % calculate conductances based on gating parameters
    gk = g_K_max .* (n^4);
    gNa = g_Na_max .* (m^3) .* h;
    gL = 0.3; % mS/cm2

    Cm = 1; % uF/cm2, Given value
    
    
    [toStim, current_stim_val] = RepetitiveStimulator(t_bounds, t, stim_vals);
    I_stim = (toStim) .* current_stim_val ;
    % Similar to a_n and a_m, create a conditional value by multiplying the
    % default by a conditional false.

    % Hodgkin Huxley Model
    dVmdt =  (-gNa*(Vm - Em_Na) - gk*(Vm-Em_K)...
        - gL*(Vm-E_L) + I_stim)/Cm;

    % Gating Parameters
    % Using the n, m, and h from y
    vm = Vm - Vr;

    %n
    a_n_default_val = ( 0.01*(10-vm) )/( exp((10-vm)/(10)) - 1);
    a_n = @ (vm) ( (vm ~= 10) .* (a_n_default_val-0.1) ) + 0.1;
    % If vm =/= 10 or 25, then the first set of parenthesis evaluates to 0. So
    % the output of the function is 0.1 because of the +0.1 term at the end. If
    % vm == 10 or 25, the parenthesis evaluates to 1 * the default val. The
    % -0.1 and +0.1 cancel out.
    B_n = 0.125 * exp(-vm/80);
    dndt = a_n(vm) * (1-n) - B_n * n;
    
    % m
    a_m_default_val = ( 0.1*(25-vm) )/( exp((25-vm)/(10)) - 1);
    a_m = @(vm) ( (vm ~= 25) .* (a_m_default_val-1) ) + 1;
    % Similar to the a_n function, this acheives a conditional output in an
    % anonymous function by multiplying the default value by a conditional
    % false.
    B_m = 4 * exp(-vm/18);
    dmdt = a_m(vm) .* (1-m) - B_m * m;
    
    % h
    a_h = 0.07 * exp(-vm/20);
    B_h = (exp( (30-vm)/10 ) + 1)^-1;
    dhdt = a_h * (1-h) - B_h * h;
    
    hh = [dVmdt; dndt; dmdt; dhdt];
    
end