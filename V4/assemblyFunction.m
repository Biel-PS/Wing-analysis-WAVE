function [K_global,F_global] = assemblyFunction(data,Td,Kel,fel)

F_global = zeros(data.ndof,1);
K_global = zeros(data.ndof,data.ndof);

    for e = 1:data.nel
        for i = 1:(data.nne*data.ni)
            F_global(Td(e,i)) = F_global(Td(e,i)) + fel(i,e);
            for j = 1:(data.nne*data.ni)
                K_global(Td(e,i),Td(e,j)) = K_global(Td(e,i),Td(e,j)) + Kel(i,j,e);
            end
        end
    end 
end