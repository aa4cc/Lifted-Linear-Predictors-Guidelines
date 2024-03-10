function out = bound_val(in,bound_min, bound_max)

out = max(min(   in     , bound_max),bound_min);

end

