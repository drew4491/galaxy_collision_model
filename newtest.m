function [t, r, v, cp1, cp2] = newtest(tmax, level, r1, r2, v1, v2, m1, m2)
  
%initializing values to be used
    nt = 2^level + 1;
%3d arrays
    r = zeros(2, 3, nt);
    v = zeros(2, 3, nt);
    t = linspace(0.0, tmax, nt);
    dt = t(2) - t(1);
%setting up poistion and velocity vectors.    
    r(1,:,1) = r1;
    r(2,:,1) = r2;
    v(1,:,1) = v1;
    v(2,:,1) = v2;
    
%initializing position after second time step 
    r(1,:,2) = r1 + dt*v1 + ((1/2)*dt^2*(m2/norm(r1-r2)^3)*(r2-r1));
    r(2,:,2) = r2 + dt*v2 + ((1/2)*dt^2*(m1/norm(r2-r1)^3)*(r1-r2));
    
%initializing cores positions
    cp1 = [[r1(1),r1(2)]; zeros(nt-2,2)];
    cp2 = [[r2(1),r2(2)]; zeros(nt-2,2)];
    
    for n = 2 : nt - 1
        d1 = norm(r(1,:,n) - r(2,:,n));
        d2 = norm(r(2,:,n) - r(1,:,n));
       
        pos1 = r(2,:,n) - r(1,:,n);
        pos2 = r(1,:,n) - r(2,:,n);
        
        r(1,:,n+1) = 2 * r(1,:,n) - r(1,:,n-1) + (dt^2)*(m2/d1^3)*pos1;
        r(2,:,n+1) = 2 * r(2,:,n) - r(2,:,n-1) + (dt^2)*(m1/d2^3)*pos2;
        
        v(1,:,n) = (r(1,:,n+1) - r(1,:,n-1)) / (2*dt);
        v(2,:,n) = (r(2,:,n+1) - r(2,:,n-1)) / (2*dt);
        
        cp1(n,1) = r(1,1,n);
        cp1(n,2) = r(1,2,n);
        cp2(n,1) = r(2,1,n);
        cp2(n,2) = r(2,2,n);
        
    end
   
    v(1,:,nt) = 2 * v(1,:,nt-1) - v(1,:,nt-2);
    v(2,:,nt) = 2 * v(2,:,nt-1) - v(2,:,nt-2);


plot(cp1(:,1), cp1(:,2), 'r');
hold on
plot(cp2(:,1), cp2(:,2), 'g');
pbaspect([1 1 1])
hold off;

end