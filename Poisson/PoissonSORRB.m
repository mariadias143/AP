function [w,ITER]=PoissonSORRB(N, TOL)
% Finds the steady-state solution for the temperature distribution on a 
% square plate, for particular boundary conditions
% USES SUCCESSIVE OVERRELAXATION WITH RELAXATION PARAMETER p
% N is the number of grid points on each direction (including boundary
% points)
% TOL is the tolerance for the stopping citeria
% IN THE GRID, USES i FOR ROW (ASCENDENT), j FOR COLUMN (LEFT-RIGHT)
% THIS USES THE RED-BLACK ORDERING

% SET BOUNDARY VALUES 
% Temperature is zero at the top and 100 on the other boundaries
w(1,1:N)=100;  % lower boundary 
w(1:N,1)=100;  % left boundary 
w(1:N,N)=100;  % right boundary 
w(N,1:N)=0;    % boundary above
w(2:N-1,2:N-1)=50; % initial values for interior points

%RELAXATION PARAMETER 
p=2/(1+sin(pi/(N-1)));


% COMPUTE STEADY STATE SOLUTION
DIFF=TOL+1;
ITER=0;
while DIFF>TOL
    u=w; % a copy of the last iteration
    % UPDATES RED POINTS
    for i=2:N-1
        for j=2+rem(i,2):2:N-1
            w(i,j)=(1-p)*w(i,j)+p*(w(i-1,j)+w(i,j-1)+w(i,j+1)+w(i+1,j))/4;
        end
    end
     % UPDATES BLACK POINTS
    for i=2:N-1 // i=1 -> N-2 
        for j=2+rem(i+1,2):2:N-1
            w(i,j)=(1-p)*w(i,j)+p*(w(i-1,j)+w(i,j-1)+w(i,j+1)+w(i+1,j))/4;
        end
    end
    ITER=ITER+1;
    DIFF=max(max(abs(w-u)));
end

% DISPLAYS A COLOURFUL MAP OF DISTRIBUTIONS
w=flipud(w); image(w,'CDataMapping','scaled')
colorbar






