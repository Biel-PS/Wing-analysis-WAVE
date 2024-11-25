function f = pointLoads(data,Td,f,F)
    for i = 1:size(F,1)
        I = nod2dof(data.ni,F(i,1),F(i,2));
        f(I) = F(i,3) +f(I); % add point loads to the global forces vector
    end

end