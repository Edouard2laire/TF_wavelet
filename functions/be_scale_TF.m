function power = be_scale_TF(power, ref)
%be_scale_TF Scale the TF representation power based on the power in ref. 
% Ref is a vector containing 1 value per sensor used to normalized the TF
% If Ref is not provided, compoute the median power accross time for each
% sensor.

    if nargin < 2
        ref         = median(sqrt(sum(power.^2, 2)),3);
    end

    assert( size(ref,1) == size(power,1) && size(ref,2) == 1);
    power = power ./ ref;

end

