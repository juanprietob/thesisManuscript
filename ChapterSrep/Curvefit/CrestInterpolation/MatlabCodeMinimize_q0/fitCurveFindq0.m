
function [x,resnorm,residual,exitflag,output] = fitCurveFindq0(thetamax, drend, dr0, nrend, nr0)


%[x,resnorm,residual,exitflag,output] = lsqnonlin(@findbestq0, -100,-200,0);
[x,resnorm,residual,exitflag,output] = lsqnonlin(@findbestq0, -1);

    function F = findbestq0(q0)
        
        q2 = 6/thetamax^4*(2*drend*thetamax+q0*thetamax^2-6*nrend+6*nr0+4*dr0*thetamax); 
        q1 = -6/thetamax^3*(-nrend+nr0+dr0*thetamax+(q0*thetamax^2)/2+(q2*thetamax^4)/12);
        
        %q2 = 6/thetamax^4*(2*drend*thetamax+q0*thetamax^2-6*nrend+6*nr0);        
        %q1 = (nrend-nr0-(q0*thetamax^2)/2-(q2*thetamax^4)/12)*6/thetamax^3;
        
        dtheta = thetamax/20;
        theta = 0:dtheta:thetamax;
        F = q0+q1*theta + q2*theta.^2;
    end    
end 