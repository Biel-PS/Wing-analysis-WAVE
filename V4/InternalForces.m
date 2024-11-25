function [xel,Sel,Mbel,Mtel] = InternalForces(data,x,Tn,Td,Kel,u)
% Compute the loads diagrams along the beam

xel = zeros(2,data.nel);
Sel = zeros(2,data.nel);
Mbel = zeros(2,data.nel);
Mtel = zeros(2,data.nel);

    for e = 1:1:data.nel

        xel(:,e) = x(Tn(e,:));
        uel = zeros(data.nne*data.ni,1);

        for i = 1:1:data.nne*data.ni
            uel(i) = u(Td(e,i));
        end

        fel_int = Kel(:,:,e)*uel;

        Sel(1,e) = -fel_int(1);
        Sel(2,e) = fel_int(4);
        Mbel(1,e) = -fel_int(2);
        Mbel(2,e) = fel_int(5);
        Mtel(1,e) = -fel_int(3);
        Mtel(2,e) = fel_int(6);

    end

end

