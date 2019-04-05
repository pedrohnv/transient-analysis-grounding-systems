% Compile the HP_HEM library interface for use from MATLAB
% By default, only the impedance matrices calculation is compiled.
% This requires a C Standard 99 compliant compiler*
% see:
%   https://mingw-w64.org
%   http://tdm-gcc.tdragon.net
%   https://www.gnu.org/
%   https://www.mathworks.com/help/matlab/matlab_external/install-mingw-support-package.html
%
% *That means that Visual Studio C Compiler won't work due to its non-compliance
% to C Standard 99.
% 
% [not implemented] If you also wish to use the linear algebra parts of the library,
% run 'setup(true)'. Beware though that this requires that you have
% Intel's MKL installed.
%
% Author: Pedro H. N. Vieira
% https://github.com/pedrohnv/HP_HEM
function setup(linalg, interleaved)
    if nargin < 1
        linalg = false;
    end
    if verLessThan('matlab', '9.4') % running on a release < R2018a ?
        if ispc % windows?
            mex -R2017b Mcalculate_impedances.c interface_matlab.c ..\\..\\src\\electrode.c ..\\..\\cubature\\hcubature.c ..\\..\\src\\auxiliary.c -I. -I..\\..\\src -I..\\..\\cubature
            mex -R2017b Mimpedances_images.c interface_matlab.c ..\\..\\src\\electrode.c ..\\..\\cubature\\hcubature.c ..\\..\\src\\auxiliary.c -I. -I..\\..\\src -I..\\..\\cubature
            if linalg
                % TODO
            end
        else
            mex -R2017b Mcalculate_impedances.c interface_matlab.c ../../src/electrode.c ../../cubature/hcubature.c ../../src/auxiliary.c -I. -I../../src -I../../cubature
            mex -R2017b Mimpedances_images.c interface_matlab.c ../../src/electrode.c ../../cubature/hcubature.c ../../src/auxiliary.c -I. -I../../src -I../../cubature
            if linalg
                % TODO
            end
        end
    else
        % R2018a onwards has interleaved complex API (better performance),
        % but needs a flag to use it
        if ispc % windows?
            mex -R2018a Mcalculate_impedances.c interface_matlab.c ..\\..\\src\\electrode.c ..\\..\\cubature\\hcubature.c ..\\..\\src\\auxiliary.c -I. -I..\\..\\src -I..\\..\\cubature
            mex -R2018a Mimpedances_images.c interface_matlab.c ..\\..\\src\\electrode.c ..\\..\\cubature\\hcubature.c ..\\..\\src\\auxiliary.c -I. -I..\\..\\src -I..\\..\\cubature
            if linalg
                % TODO
            end
        else
            mex -R2018a Mcalculate_impedances.c interface_matlab.c ../../src/electrode.c ../../cubature/hcubature.c ../../src/auxiliary.c -I. -I../../src -I../../cubature
            mex -R2018a Mimpedances_images.c interface_matlab.c ../../src/electrode.c ../../cubature/hcubature.c ../../src/auxiliary.c -I. -I../../src -I../../cubature
            if linalg
                % TODO
            end
        end
    end
end