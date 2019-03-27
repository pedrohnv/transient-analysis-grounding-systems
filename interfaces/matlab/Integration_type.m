% Type of integration to be done.
%   NONE: exp(-gamma*r)/r * integral2(1, dls, dlr)
%   DOUBLE: integral2(exp(-gamma*r)/r, dls, dlr)
%   EXP_LOGNF: integral(exp(-gamma*r) * log(Nf), dlr)
%   LOGNF: exp(-gamma*r) * integral(log(Nf), dlr)
classdef Integration_type < int16
   enumeration
      INTG_NONE (1)
      INTG_DOUBLE (2)
      INTG_EXP_LOGNF (3)
      INTG_LOGNF (4)
   end
end