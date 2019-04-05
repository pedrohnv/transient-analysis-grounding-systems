% Splits `electrode` in `n` segments of equal length.
% 
% Parameters
% ----------
%   electrode : the electrode struct to be split
%   num_segments : number of segments into which electrode will be split
% 
% Returns
% -------
%   segments : vector of electrode struct
%   nodes : vector of unique nodes to which the segments are connected
function [segments, nodes] = segment_electrode(electrode, num_segments)
    nn = num_segments + 1;
    nodes = zeros(nn, 3);
    startp = [0., 0., 0.];
    endp = [0., 0., 0.];
    for k = 1:3
        startp(k) = electrode.start_point(k);
        endp(k) = electrode.end_point(k);
    end
    increment = (endp - startp)/num_segments;
    for k = 0:num_segments
        nodes(k+1,:) = startp + k*increment;
    end
    segments = repmat(new_electrode(), num_segments, 1); %pre-allocate
    for k = 1:num_segments
        segments(k) = new_electrode(nodes(k,:), nodes(k+1,:), electrode.radius, electrode.zi);
    end
end