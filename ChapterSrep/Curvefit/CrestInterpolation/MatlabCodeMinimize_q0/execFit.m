
%S1origin = [0.679702676853675 0.655098003973841 0.162611735194631];%rand(1,3);
S1origin = [0 0 0];
%r1 = rand(1,3);
S1r1 = [0.679702676853675 0.655098003973841 0.162611735194631];
%r2 = rand(1,3);
S1r2 = [0.118997681558377 0.498364051982143 0.959743958516081];
%r0 = rand(1,3);
S1r0 = (S1r1+S1r2)/(norm(S1r1)+norm(S1r2));

S1r1 = S1r1/norm(S1r1);
S1r0 = S1r0/norm(S1r0);

lenght = 1;

drend0 = -1;
dr0 = 0;

S1r0 = S1r0*lenght;
    
%ddr2 = fitCurve(S1r0, S1r1, drend0, drend1, dr0, dr1, S1origin);
ddr2 = fitCurve2(S1r0, S1r1, drend0, dr0, S1origin);
