% The norm that is used to measure the error and determine convergence properties
% during integration.
% -   `ERROR_L1`, `ERROR_L2`, `ERROR_LINF`: the absolute error is measured
%     as |e| and the relative error as |e|/|v|, where |...| is the
%     [L?](https://en.wikipedia.org/wiki/Taxicab_geometry),
%     [L?](https://en.wikipedia.org/wiki/Euclidean_distance), or
%     [L?](https://en.wikipedia.org/wiki/Uniform_norm)
%     [norm](https://en.wikipedia.org/wiki/Norm_(mathematics)), respectively. (|x| in the
%     L? norm is the sum of the absolute values of the components, in
%     the L? norm is the root mean square of the components, and in the
%     L? norm is the maximum absolute value of the components)
% 
% -   `ERROR_INDIVIDUAL`: Convergence is achieved only when each integrand
%     (each component of v and e) individually satisfies the requested
%     error tolerances.
% 
% -   `ERROR_PAIRED`: Like `ERROR_INDIVIDUAL`, except that the integrands
%     are grouped into consecutive pairs, with the error tolerance applied
%     in an L? sense to each pair. This option is mainly useful for
%     integrating vectors of complex numbers, where each consecutive pair
%     of real integrands is the real and imaginary parts of a single complex
%     integrand, and you only care about the error in the complex plane rather
%     than the error in the real and imaginary parts separately.
classdef Error_norm < int16
   enumeration
      INDIVIDUAL (0)
      PAIRED (1)
      L2 (2)
      L1 (3)
      LINF (4)
   end
end