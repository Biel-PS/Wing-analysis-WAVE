function [up,vp] = applyBC(data,p)
    data.np = size(p,1);
    vp = zeros(data.np,1);
    up = zeros(data.np,1);
    for i = 1:data.np
        vp(i) = nod2dof(data.ni,p(i,1),p(i,2)); % DOF index
        up(i) = p(i,3);                         % DOF value
    end
end