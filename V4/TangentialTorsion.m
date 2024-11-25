function [tau_t,s] = TangentialTorsion(data,x,Tn,m,Tm,Mz_p,J,Ain,open)

s = zeros(2,data.nel);
tau_t = zeros(2,data.nel);

    for e = 1:1:data.nel

        t = m(Tm(e),1);

        delta_x = x(Tn(e,2),1) - x(Tn(e,1),1);
        delta_y = x(Tn(e,2),2) - x(Tn(e,1),2);
        l = sqrt(delta_x^2+delta_y^2);

        if e>1
            s(1,e) = s(2,e-1);
        end
        s(2,e) = s(1,e) + l;

        if open == true
            tau_t(1,e) = Mz_p*t/J;
            tau_t(2,e) = Mz_p*t/J;
        else
            tau_t(1,e) = Mz_p/(2*Ain*t);
            tau_t(2,e) = Mz_p/(2*Ain*t);
        end

    end

end