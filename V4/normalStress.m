function [sigma,s] = normalStress(data,x,Tn,Xo,Yo,Ixx,Iyy,Ixy,Mx_p,My_p)

s = zeros(2,data.nel);
sigma = zeros(2,data.nel);

Ix = Ixx - Ixy^2/Iyy;
Iy = Iyy - Ixy^2/Ixx;
Mx = Mx_p + My_p*Ixy/Iyy;
My = My_p + Mx_p*Ixy/Ixx;

    for e = 1:1:data.nel

        delta_x = x(Tn(e,2),1) - x(Tn(e,1),1);
        delta_y = x(Tn(e,2),2) - x(Tn(e,1),2);
        l = sqrt(delta_x^2+delta_y^2);

        if e>1
            s(1,e) = s(2,e-1);
        end
        s(2,e) = s(1,e) + l;

        sigma(1,e) = (x(Tn(e,1),2) - Yo)*Mx/Ix - (x(Tn(e,1),1)-Xo)*My/Iy;
        sigma(2,e) = (x(Tn(e,2),2) - Yo)*Mx/Ix - (x(Tn(e,2),1)-Xo)*My/Iy;
    end

end

