function Td = connectDOF(data,Tn)
Td = zeros(data.nel,data.ni*2);

    if data.ni == 1
        Td = Tn;
    elseif data.ni == 2
        for ii = 1:data.nel
            for jj = 1:data.ni
                Td(ii,jj*2-1) = nod2dof(data.ni,Tn(ii,jj),1);
                Td(ii,jj*2) = nod2dof(data.ni,Tn(ii,jj),2);
            end
        end
    elseif data.ni == 3
        for ii = 1:data.nel
            for jj = 1:(data.ni-1)
                Td(ii,jj*3-2) = nod2dof(data.ni,Tn(ii,jj),1);
                Td(ii,jj*3-1) = nod2dof(data.ni,Tn(ii,jj),2);
                Td(ii,jj*3) = nod2dof(data.ni,Tn(ii,jj),3);
            end
        end
    end


end