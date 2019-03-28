% Type of integration to be done.
%   NONE: exp(-gamma*rbar)/rbar * ls * lr
%   DOUBLE: integral2(exp(-gamma*r)/r, dls, dlr)
%   EXP_LOGNF: integral(exp(-gamma*rbar) * log(Nf), dlr)
%   LOGNF: exp(-gamma*rbar) * integral(log(Nf), dlr)
classdef Integration_type < int16
   enumeration
      NONE (1)
      DOUBLE (2)
      EXP_LOGNF (3)
      LOGNF (4)
   end
end