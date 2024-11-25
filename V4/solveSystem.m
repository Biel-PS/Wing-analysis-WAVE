function [u,r] = solveSystem(data,K,f,up,vp)
    vf = setdiff((1:data.ndof)',vp);
    u = zeros(data.ndof,1);
    u(vp) = up;
    u(vf) = inv(K(vf,vf)) * (f(vf) - K(vf,vp)*u(vp)); % free DOFs
    r = K(vp,:)*u-f(vp); % reaction loads
end