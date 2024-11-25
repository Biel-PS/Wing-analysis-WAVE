function I = nod2dof(ni,i,j)
% Insert formula that, for a problem with 'ni' degrees of freedom per node,
% returns the global degree of freedom index 'I' corresponding to the
% degree of freedom 'j' (1<j<ni) of node 'i'.

    I = ni*i-(ni-j);

end