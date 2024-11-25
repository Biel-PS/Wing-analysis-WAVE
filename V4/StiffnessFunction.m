function Kel = StiffnessFunction(data_b,x,Tn,m,Tm,J,I)

Kel = zeros(data_b.nne*data_b.ni,data_b.nne*data_b.ni,data_b.nel);

    for e = 1:1:data_b.nel
    
        l = abs(x(Tn(e,2),1)-x(Tn(e,1),1));

        E = m(Tm(e),1);
        G = m(Tm(e),3);


        Kb = (E*I/l^3)*[ 12 6*l 0 -12 6*l 0;
            6*l 4*l^2 0 -6*l 2*l^2 0;
            0 0 0 0 0 0;
            -12 -6*l 0 12 -6*l 0;
            6*l 2*l^2 0 -6*l 4*l^2 0;
            0 0 0 0 0 0;
        ];

        Kt = (G*J/l)*[ 0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 1 0 0 -1;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 -1 0 0 1;
        ];

        Kel(:,:,e) = Kb + Kt;

    end

end

