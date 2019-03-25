function [el] = new_electrode(start_point, end_point, radius, zi)
    % Creates an Electrode struct
    
    % TODO argument (error) checking
    if nargin == 0
        el.start_point = [0.0, 0.0, 0.0];
        el.middle_point = [0.0, 0.0, 0.0];
        el.end_point = [0.0, 0.0, 0.0];
        el.length = 0.0;
        el.radius = 0.0;
    else
        el.start_point = start_point;
        el.middle_point = (start_point + end_point)/2;
        el.end_point = end_point;
        el.length = norm(start_point - end_point);
        el.radius = radius;
    end
    if nargin < 4
        el.zi = 0.0j; % internal impedance
    else
        el.zi = zi;
    end
end