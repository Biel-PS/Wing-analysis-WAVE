function sig = stressFunction(data,x,Tn,m,Tm,Td,u)
    sig = zeros (data.nel,1);
    
    for e = 1:data.nel
        uel = zeros(data.nne * data.ni,1);
        for i = 1:(data.nne*data.ni)
            uel(i) = u(Td(e,i));
        end
                
        x_el = x(Tn(e,:),:);

        d_x = x_el(2,1) - x_el(1,1);
        d_y = x_el(2,2) - x_el(1,2);
        d_z = x_el(2,3) - x_el(1,3);

        Le = sqrt(d_x^2 + d_y^2 + d_z^2);
       

        E = m(Tm(e),1);
        A = m(Tm(e),2);

        R = 1/Le*[d_x d_y d_z 0 0 0; 0 0 0 d_x d_y d_z];

        E = m(Tm(e),1); 

        epsi = (1/Le)*[-1 1]*R*uel; % element strain
        sig(e) = E*epsi; % element stress
    end
   
end