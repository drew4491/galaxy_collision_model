function [t, r, v, xcore] = test(tmax, level, r1, r2, v1, v2, m1, m2)
  
%initializing values to be used

    nt = (2^level) + 1;
%3d arrays
    xcore = zeros(nt,1);
    r = zeros(2, 3, nt);
    v = zeros(2, 3, nt);
    t = linspace(0, tmax, nt);
    dt = t(2) - t(1);
%setting up poistion and velocity vectors.    
    r(1,:,1) = r1;
    r(2,:,1) = r2;
    v(1,:,1) = v1;
    v(2,:,1) = v2;
    
%initializing position after second time step 
    r(1,:,2) = r1 + dt*v1 + ((1/2)*(dt^2)*(m2/norm(r1-r2)^3)*(r2-r1));
    r(2,:,2) = r2 + dt*v2 + ((1/2)*(dt^2)*(m1/norm(r2-r1)^3)*(r1-r2));
    
    for n = 2 : nt
        d1 = norm(r(1,:,n) - r(2,:,n));
        d2 = norm(r(2,:,n) - r(1,:,n));
        pos1 = r(2,:,n) - r(1,:,n);
        pos2 = r(1,:,n) - r(2,:,n);
        r(1,:,n+1) = 2*r(1,:,n) - r(1,:,n-1) + (dt^2)*(m2/d1^3)*pos1;
        r(2,:,n+1) = 2*r(2,:,n) - r(2,:,n-1) + (dt^2)*(m1/d2^3)*pos2;
        v(1,:,n) = (r(1,:,n+1) - r(1,:,n-1)) / (2*dt);
        v(2,:,n) = (r(2,:,n+1) - r(2,:,n-1)) / (2*dt);
        
        xcore(n,1) = r(1,1,n);
        
    end
    
    v(1,:,nt) = 2 * v(1,:,nt-1) - v(1,:,nt-2);
    v(2,:,nt) = 2 * v(2,:,nt-1) - v(2,:,nt-2);
end
