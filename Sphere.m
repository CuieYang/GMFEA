function obj = Sphere(var,M,opt)
%Sphere function
%   - var: design variable vector
%   - opt: shift vector
    var = (M*(var-opt)')';
    obj=var*var';
end