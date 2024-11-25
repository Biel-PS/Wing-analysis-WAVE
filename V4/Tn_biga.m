function [Tn_b,Tm_b] = Tn_biga (nnodes)
    Tn_b= zeros(nnodes-1,2);
    Tm_b = zeros(nnodes-1,1);
    for i = 1:nnodes-1
        Tn_b(i,:) = [i,i+1];
        Tm_b(i) = 1;
    end
end