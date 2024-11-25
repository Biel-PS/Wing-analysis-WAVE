function [tau_s,s] = TangentialShear(data,x,Tn,m,Tm,Xo,Yo,Xs,Ys,Ain,Ixx,Iyy,Ixy,Sx_p,Sy_p,open)

s = zeros(2,data.nel);
tau_s = zeros(2,data.nel);

    if open == true
        q = 0;
    else
        q = (Sy_p*(Xo-Xs)-Sx_p*(Yo-Ys))/(2*Ain);
    end

    Ix = Ixx - Ixy^2/Iyy;
    Iy = Iyy - Ixy^2/Ixx;
    Sx = Sx_p - Sy_p*Ixy/Ixx;
    Sy = Sy_p - Sx_p*Ixy/Iyy;

    for e = 1:1:data.nel

        t = m(Tm(e),1);

        delta_x = x(Tn(e,2),1) - x(Tn(e,1),1);
        delta_y = x(Tn(e,2),2) - x(Tn(e,1),2);
        l = sqrt(delta_x^2+delta_y^2);

        if e>1
            s(1,e) = s(2,e-1);
        end
        s(2,e) = s(1,e) + l;

        tau_s(1,e) = q/t; % shear stress at first node
        q = q - Sx*t*l*(delta_x/2 + x(Tn(e,1),1)-Xo)/Iy - Sy*t*l*(delta_y/2 + x(Tn(e,1),2)-Yo)/Ix;
        tau_s(2,e) = q/t; % shear stress at second node

    end


end

