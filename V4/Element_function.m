function [fe,me,x_nod] = Element_function (be,b,ze,za,zm,data,v_inf,rho,Cl,nnodes,W,We,xi_S)
    
    delta_x = b/(nnodes-1);
    x_nod = zeros(nnodes,3);
    x = 0:delta_x:b;
    x_nod(:,1) = x;
    

    l = 0.5 * rho * v_inf^2 * data.c * Cl * sqrt(1-(x./b).^2);
    fn = zeros (3,nnodes);
    fe = zeros (3,nnodes-1);
    me = zeros (1,nnodes-1);



    for i = 1:nnodes
        fn(2,i) = -W + l(i);
    end

    
    for i = 1:nnodes-1
        fe(:,i) = (fn(:,i)+fn(:,i+1))/2;
        me(i) = (-2*W*(xi_S-zm)+(l(i)+l(i+1))*(xi_S-za))/2; 
    end
    % 
    % index_we = find(x >= (be-delta_x/2) & x <= (be+delta_x/2));
    % fn(2,index_we(1)) = fn(2,index_we(1)) - We;
    % 
    % fe(:,index_we-1) = (fn(:,index_we-1)+fn(:,index_we))/2;
    % fe(:,index_we) = (fn(:,index_we)+fn(:,index_we+1))/2;
    % 
    % me(index_we-1) = me(index_we-1) - We*(xi_S-ze);
    % me(index_we) = me(index_we) - We*(xi_S-ze); 

end