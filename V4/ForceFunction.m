function fel = ForceFunction(data,x,Tn,felv,mel)

fel = zeros(data.nne*data.ni,data.nel);

    for e = 1:1:data.nel
    
        l = abs(x(Tn(e,2),1)-x(Tn(e,1),1));

        fb = felv(2,e)*l*[1/2; l/12; 0; 1/2; -l/12; 0];
        ft = mel(e)*l*[0; 0; 1/2; 0; 0; 1/2];

        fel(:,e) = fb + ft;

    end

end

