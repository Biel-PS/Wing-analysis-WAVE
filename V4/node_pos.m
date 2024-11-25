function [x,Tn,Tm] = node_pos (h1,h2,d,nod,open)

	x = zeros(nod , 2);
    if open == 0
        Tn = zeros (nod,2);
        Tm = zeros (nod,1);
    else
        Tn = zeros (nod-1,2);
        Tm = zeros (nod-1,1);
    end

    p = 0;
	done = false;
    j = 2;
	increment_d = 4*d/nod;
	increment_h1 = 4*h1/nod;
	increment_h2 = 4*h2/nod;
	
	x(1,2) = -h2/2;

	for i = 2:nod
        if ((i-1) <= nod/4)
		    x(i,1) = x(i-1,1)+increment_d;
		    x(i,2) = -(h1-h2)/(2*d)*x(i,1) - h2/2;
            Tm(j-1,1) = 3;
            
				    
        elseif (i-1) > nod/4 && (i-1) <= 2*nod/4
		    x(i,1) = x(i-1,1);
		    x(i,2) = x(i-1,2) + increment_h1;
            Tm(j-1,1) = 1;
            if x(i,2)>= 0 && x(i-1,2)<= 0 && done == false
               p = i;
            end
				    
        elseif (i-1) > 2*nod/4 && (i-1) <= 3*nod/4
			x(i,1) = x(i-1,1)-increment_d;
			x(i,2) = (h1-h2)/(2*d)*x(i,1) + h2/2;
            Tm(j-1,1) = 3;
			    
        elseif (i-1) > 3*nod/4 && (i-1) <= 4*nod/4
			x(i,1) = x(i-1,1);
			x(i,2) = x(i-1,2) - increment_h2;
            Tm(j-1,1) = 2;
			    
        end
         Tn(j-1,1) = i-1;
         Tn(j-1,2) = i;

        if open && p == i
            j=j-1;
            done = true;
        end
           
            j = j+1;
    end
    if done == false
        Tn(nod,1) = nod;
        Tn(nod,2) = 1;
        Tm(nod,1) = 2;
    else 
        Tn(nod-1,1) = nod;
        Tn(nod-1,2) = 1;
        Tm(nod-1,1) = 2;
    end
end