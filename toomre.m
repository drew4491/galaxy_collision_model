function [t, r, v, BIG] = toomre(tmax, level, r1, r2, v1, v2, m1, m2, rot1, rot2, N1, N2)
  
%initializing values to be used

    nt = (2^level) + 1;
%3d arrays
    N = N1 + N2;
    r = zeros(2, 3, nt);
    v = zeros(N+2, 3, nt);
    t = linspace(0.0, tmax, nt);
    dt = t(2) - t(1);
    
%setting up poistion and velocity vectors.    
    r(1,:,1) = r1;
    r(2,:,1) = r2;
    v(1,:,1) = v1;
    v(2,:,1) = v2;
    
%initializing position after second time step 
    r(1,:,2) = r1 + dt*v1 + (1/2)*dt^2*(m2/(norm(r1-r2)^3))*(r2-r1);
    r(2,:,2) = r2 + dt*v2 + (1/2)*dt^2*(m1/(norm(r2-r1)^3))*(r1-r2);
    
    cp1 = [[r1(1) r1(2)]; zeros(nt-2,2)];
    cp2 = [[r2(1) r2(2)]; zeros(nt-2,2)];
%stars initialization
%star 1    

    theta1 = 2*pi*rand(N1,1);
    dist1 =  4 * rand((N1),1) + 0.06*m1 + 0.55;
    
    xc1 = r1(1) + dist1.*cos(theta1);
    yc1 = r1(2) + dist1.*sin(theta1);
    
    stars1 = zeros(N1, 3, nt);
    stars1(:,1,1) = xc1;
    stars1(:,2,1) = yc1;
    stars1(:,3,1) = zeros(length(xc1),1);
    

%star 2    

    theta2 = 2*pi*rand(N2,1);
    dist2 =  4 * rand(N2,1) + 0.06*m2 + 0.55;
    
    xc2 = r2(1) + dist2.*cos(theta2);
    yc2 = r2(2) + dist2.*sin(theta2);
    
    stars2 = zeros(N2, 3, nt);
    stars2(:,1,1) = xc2;
    stars2(:,2,1) = yc2;
    stars2(:,3,1) = zeros(length(xc2),1);
    
    
    
%matrix positions of cores and stars
    BIG = [r;stars1;stars2];
    
%velocity of stars
    
    if rot1 >= 1
        vs1x = -sqrt(m1./dist1).*sin(theta1) + v1(1);
        vs1y = sqrt(m1./dist1).*cos(theta1) + v1(2);
    else
        vs1x = sqrt(m1./dist1).*sin(theta1) + v1(1);
        vs1y = -sqrt(m1./dist1).*cos(theta1) + v1(2);
    end
    
    if rot2 >= 1
        vs2x = -sqrt(m2./dist2).*sin(theta2) + v2(1);
        vs2y = sqrt(m2./dist2).*cos(theta2) + v2(2);
    else
        vs2x = sqrt(m2./dist2).*sin(theta2) + v2(1);
        vs2y = -sqrt(m2./dist2).*cos(theta2) + v2(2);

    end
    
    vstar1 = [vs1x vs1y zeros(length(vs1x),1)];
    vstar2 = [vs2x vs2y zeros(length(vs2x),1)];
    
    vall(:,:,1) = [v1; v2; vstar1; vstar2];
    
    d1s = zeros(N,1,nt);
    d2s = zeros(N,1,nt);
            
    dsc1 = zeros(N,3,nt);         
    dsc2 = zeros(N,3,nt);  
    for a = 3:N+2
        
        sep1 = norm(BIG(1,:,1) - BIG(a,:,1));
        sep2 = norm(BIG(2,:,1) - BIG(a,:,1));
        
        dist1 = BIG(1,:,1) - BIG(a,:,1);
        dist2 = BIG(2,:,1) - BIG(a,:,1);
        
        BIG(a,:,2) = BIG(a,:,1) + dt*vall(a,:,1) + (1/2)*(dt^2)*(((m1/sep1^3)*dist1) + ((m2/sep2^3)*dist2));
    end
    
    for n = 2 : nt - 1
        d1 = norm(BIG(1,:,n) - BIG(2,:,n));
        d2 = norm(BIG(2,:,n) - BIG(1,:,n));
        
        pos1 = BIG(2,:,n) - BIG(1,:,n);
        pos2 = BIG(1,:,n) - BIG(2,:,n);
        
        BIG(1,:,n+1) = 2 * BIG(1,:,n) - BIG(1,:,n-1) + (dt^2)*(m2/d1^3)*pos1;
        BIG(2,:,n+1) = 2 * BIG(2,:,n) - BIG(2,:,n-1) + (dt^2)*(m1/d2^3)*pos2;
        
        cp1(n,1) = BIG(1,1,n);
        cp1(n,2) = BIG(1,2,n);
        cp2(n,1) = BIG(2,1,n);
        cp2(n,2) = BIG(2,2,n);
        
        dsc1(:,:,n) = BIG(1,:,n) - BIG(3:end,:,n);
        dsc2(:,:,n) = BIG(2,:,n) - BIG(3:end,:,n);
        
        d1s(:,:,n) = vecnorm(dsc1(:,:,n),2,2);
        d2s(:,:,n) = vecnorm(dsc2(:,:,n),2,2);
            
        
            
        BIG(3:end,:,n+1) = 2*BIG(3:end,:,n) - BIG(3:end,:,n-1) + (dt^2)*((m1./d1s(:,:,n).^3).*dsc1(:,:,n) + (m2./d2s(:,:,n).^3).*dsc2(:,:,n));
        
    end
    
    avifilename = 'galaxycw.avi';
    aviobj = VideoWriter(avifilename);
    open(aviobj);
    
    for n = 1: nt
        if mod(n,5) == 0
            clf;
            hold on;
            axis square;
            box on;
            set(gca,'color','k');
            set(gcf,'color','k');
            xlim([-40, 40]);
            ylim([-40, 40]);
            plot(BIG(1,1,n), BIG(1,2,n),'m.', BIG(2,1,n), BIG(2,2,n), ...
                  'b.', BIG(3:(N1)+2,1,n), BIG(3:(N1)+2,2,n), 'yp', ...
                  BIG((N2)+3:N,1,n), BIG((N2)+3:N,2,n), 'cp');
            drawnow;
      
                writeVideo(aviobj, getframe(gcf));
        end
        
   
    end    
    close(aviobj);
end

% [t, BIG] = step1(40, 8, [-4,-9,0], [4,9,0], [0.9,0.4,0], [-0.9,-0.4,0], 8, 8, 6000);