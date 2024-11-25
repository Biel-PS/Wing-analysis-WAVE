function [Xo,Yo,Xs,Ys,Atot,Ixx,Iyy,Ixy,J,Ain] = SectionProperties(data,x,Tn,m,Tm,open)

l = zeros(data.nel,1);
A = zeros(data.nel,1);
x_c = zeros(data.nel,2);

Atot = 0;
Xo = 0;
Yo = 0;

if open == 0
    tram0 = round(2* data.nel/4 + 1* data.nel/4)/2;
    tram2 = 4* data.nel/4;
end

if open == 1
    tram0 = round(((data.nel/4 +1)+(2* data.nel/4 -1))/2);
    tram2 = round(4* data.nel/4);

end
    for e = 1:1:data.nel

        t = m(Tm(e),1);

        delta_x = x(Tn(e,2),1) - x(Tn(e,1),1);
        delta_y = x(Tn(e,2),2) - x(Tn(e,1),2);
        l(e) = sqrt(delta_x^2+delta_y^2);
        A(e) = t*l(e);
        x_c(e,:) = (x(Tn(e,2),:)+x(Tn(e,1),:))/2;

        Atot = Atot + A(e);
        Xo = Xo + A(e)*x_c(e,1);
        Yo = Yo + A(e)*x_c(e,2);

    end

    Xo = Xo/Atot;
    Yo = Yo/Atot;

    Ixx = 0;
    Iyy = 0;
    Ixy = 0;
    J = 0;
    Ain = 0; % only for closed

    for e = 1:1:data.nel

        t = m(Tm(e),1);

        delta_x = x(Tn(e,2),1) - x(Tn(e,1),1);
        delta_y = x(Tn(e,2),2) - x(Tn(e,1),2);
        Ixx = Ixx + A(e)*delta_y^2/12 + A(e)*(x_c(e,2)-Yo)^2;
        Iyy = Iyy + A(e)*delta_x^2/12 + A(e)*(x_c(e,1)-Xo)^2;
        Ixy = Ixy + A(e)*delta_x*delta_y/12 + A(e)*(x_c(e,1)-Xo)*(x_c(e,2)-Yo);

        if open == 1
            J = J + l(e)*t^3/3;
            Ain = -1;
        else 
            Ain = Ain + norm(cross([x(Tn(e,1),1) - Xo, x(Tn(e,1),2) - Yo, 0],[delta_x, delta_y, 0]))/2;
            J = J + l(e)/t;
        end
    end 

    if open == 0
        J = 4*Ain^2/J;
    end

    q1 = 0; q2 = 0;
    Xs = Xo; Ys = Yo;
    

    for e = tram0:1:tram2
        
        t = m(Tm(e),1);

        delta_x = x(Tn(e,2),1) - x(Tn(e,1),1);
        delta_y = x(Tn(e,2),2) - x(Tn(e,1),2);
        
        A1 = ( Iyy*delta_y/2 - Ixy*delta_x/2 ) / ( Ixx*Iyy-Ixy^2 );

        B1 = ( Iyy*(x(Tn(e,1),2)-Yo) - Ixy*(x(Tn(e,1),1)-Xo) ) / (Ixx*Iyy - Ixy^2);

        A2 = ( Ixx*delta_x/2 - Ixy*delta_y/2) / (Ixx*Iyy - Ixy^2);

        B2 = ( Ixx*(x(Tn(e,1),1)-Xo) - Ixy*(x(Tn(e,1),2) - Yo)) / (Ixx*Iyy - Ixy^2);

        C = ( x(Tn(e,1),1)-Xo ) * delta_y - ( x(Tn(e,1),2)-Yo ) * delta_x;



        Xs = Xs + C*(q1 - t*l(e)*(A1/3+B1/2));
        Ys = Ys + C*( q2 + t*l(e)*(A2/3+B2/2));

        q1 = q1 - t*l(e)*(A1 + B1);
        q2 = q2 + t*l(e)*(A2 + B2);

        
    end
    for e = 1:1:(tram0-1)
            
            t = m(Tm(e),1);
    
            delta_x = x(Tn(e,2),1) - x(Tn(e,1),1);
            delta_y = x(Tn(e,2),2) - x(Tn(e,1),2);
            
            A1 = ( Iyy*delta_y/2 - Ixy*delta_x/2 ) / ( Ixx*Iyy-Ixy^2 );
    
            B1 = ( Iyy*(x(Tn(e,1),2)-Yo) - Ixy*(x(Tn(e,1),1)-Xo) ) / (Ixx*Iyy - Ixy^2);
    
            A2 = ( Ixx*delta_x/2 - Ixy*delta_y/2) / (Ixx*Iyy - Ixy^2);
    
            B2 = ( Ixx*(x(Tn(e,1),1)-Xo) - Ixy*(x(Tn(e,1),2) - Yo)) / (Ixx*Iyy - Ixy^2);
    
            C = ( x(Tn(e,1),1)-Xo ) * delta_y - ( x(Tn(e,1),2)-Yo ) * delta_x;
    
    
    
            Xs = Xs + C*(q1 - t*l(e)*(A1/3+B1/2));
            Ys = Ys + C*( q2 + t*l(e)*(A2/3+B2/2));
    
            q1 = q1 - t*l(e)*(A1 + B1);
            q2 = q2 + t*l(e)*(A2 + B2);
    end
    if open == 0
        Xs = Xo;
        Ys = Yo,
    
    end

end

