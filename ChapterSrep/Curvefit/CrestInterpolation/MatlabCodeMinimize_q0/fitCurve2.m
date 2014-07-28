

function ddr = fitCurve2(r0, rend, drend, dr0, origin)


    hold on;
    plot3([origin(1),origin(1)+r0(1)],[origin(2),origin(2)+r0(2)],[origin(3),origin(3)+r0(3)],'-r');    
    plot3([origin(1),origin(1)+rend(1)],[origin(2),origin(2)+rend(2)],[origin(3),origin(3)+rend(3)],'--b');
    hold off;

    thetamax = acos(dot(r0, rend) / (norm(r0)* norm(rend)));
    
    nrend = norm(rend);    
    
    nr0 = norm(r0);    
    r0n = r0/nr0; 
    
    
    [q0,resnorm,residual,exitflag,output] = fitCurveFindq0(thetamax, drend,dr0, nrend, nr0);   
    
    q2 = 6/thetamax^4*(2*drend*thetamax+q0*thetamax^2-6*nrend+6*nr0+4*dr0*thetamax); 
    q1 = -6/thetamax^3*(-nrend+nr0+dr0*thetamax+(q0*thetamax^2)/2+(q2*thetamax^4)/12);
    
    %q2 = 6/thetamax^4*(2*drend*thetamax+q0*thetamax^2-6*nrend+6*nr0);        
    %q1 = (nrend-nr0-(q0*thetamax^2)/2-(q2*thetamax^4)/12)*6/thetamax^3;
    
    theta = 0:thetamax/20:(thetamax);
    
    ddr = q0 + q1*theta + q2*theta.^2;
    
    graphplot = 0;
    
    if(graphplot == 0)

        hold on;
        
        rotvect = cross(rend,r0);
        rotvect = rotvect/norm(rotvect);
        
        sizetheta = size(theta, 2);
        
        for i=1:sizetheta

            r = norm(r0) + dr0*theta(i) + (q0*theta(i)^2)/2 + (q1*theta(i)^3)/6 + (q2*theta(i)^4)/12;
            
            rotmat = rotationMatrix(rotvect, theta(1,i));
            newvect = (r0n*rotmat);
            newvect = newvect*r;
            newvect = origin + newvect;
            plot3(newvect(1), newvect(2), newvect(3),'og');


        end

        hold off;
    end
end


function rotmat = rotationMatrix(rotvect, theta)

ux = rotvect(1,1);
uy = rotvect(1,2);
uz = rotvect(1,3);
costheta = cos(theta);
sintheta = sin(theta);
rotmat = [costheta + (ux^2)*(1-costheta), ux*uy*(1-costheta)-uz*sintheta, ux*uz*(1-costheta)+uy*sintheta;
          uy*ux*(1-costheta)+uz*sintheta, costheta+(uy^2)*(1-costheta), uy*uz*(1-costheta)-ux*sintheta;
          uz*ux*(1-costheta)-uy*sintheta, uz*uy*(1-costheta)+ux*sintheta, costheta+(uz^2)*(1-costheta)];

end