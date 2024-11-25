function [sigVM,posmax] = VonMises(open,h1,h2,d,nnodes,data,m,Sel,Mbel,Mtel)

[x, Tn, Tm] = node_pos(h1,h2,d,nnodes,open);

data.nnod = size(x,1); % Number of nodes 
data.nd = size(x,2);   % Problem dimension
data.ndof = data.nnod*data.ni;  % Total number of degrees of freedom
data.nel = size(Tn,1); % Number of elements 
data.nne = size(Tn,2); % Number of nodes in a bar

[Xo,Yo,Xs,Ys,Atot,Ixx,Iyy,Ixy,J,Ain] = SectionProperties(data,x,Tn,m,Tm,open);

sigVM = zeros(1,size(Sel,2));
posmax = zeros(2,size(Sel,2));
for bnod = 1:size(Sel,2)
    Mx_p = -Mbel(1,bnod);
    My_p = 0;
    [sigma,~] = normalStress(data,x,Tn,Xo,Yo,Ixx,Iyy,Ixy,Mx_p,My_p);
    
    Sx_p = 0; Sy_p = Sel(1,bnod);
    [tau_s,~] = TangentialShear(data,x,Tn,m,Tm,Xo,Yo,Xs,Ys,Ain,Ixx,Iyy,Ixy,Sx_p,Sy_p,open);
    
    Mz_p = Mtel(1,bnod);
    [tau_t,~] = TangentialTorsion(data,x,Tn,m,Tm,Mz_p,J,Ain,open);
    
    
        for e = 1:data.nel
            control = sqrt(sigma(1,e)^2+3*(tau_s(1,e)+tau_t(1,e))^2);
            if control > sigVM(1,bnod)
                sigVM(1,bnod) = control;
                posmax(:,bnod) = [x(e,1); x(e,2)];
            end
        end
end

end
