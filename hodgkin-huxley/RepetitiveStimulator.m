function [toStimulate_overall, current_stim_val] = RepetitiveStimulator(t_bounds, t, stim_vals)
    toStimulate = [];
    
    for i=1:2:length(t_bounds)
        stim_decision = (t > t_bounds(i) & t < t_bounds(i+1));
        toStimulate = [toStimulate, stim_decision];
    end
    toStimulate_overall = any(toStimulate);
    
    current_stim_val = stim_vals(logical(toStimulate));
    
    if isempty(current_stim_val)
       current_stim_val = 0;
    end
end