[t8, r8, v8, xcore8] = test(40, 8, [1/3,0,0], [-2/3,0,0], [0,sqrt(1/6),0], [0,-sqrt(2/3),0], 1, 0.5);
[t9, r9, v9, xcore9] = test(40, 9, [1/3,0,0], [-2/3,0,0], [0,sqrt(1/6),0], [0,-sqrt(2/3),0], 1, 0.5);
[t10, r10, v10, xcore10] = test(40, 10, [1/3,0,0], [-2/3,0,0], [0,sqrt(1/6),0], [0,-sqrt(2/3),0], 1, 0.5);

t8 = t8';
t9 = t9';
t10 = t10';

figure()
plot(t8,xcore8,'k-');
hold on;
plot(t9,xcore9,'b-');
hold on;
plot(t10,xcore10,'g-');
legend('xcore8', 'xcore9', 'xcore10');
pbaspect([1 1 1])
hold off;

xcore9 = xcore9(1:2:end);
xcore10 = xcore10(1:4:end);


dxcore89 = xcore8 - xcore9;
dxcore910 = xcore9 - xcore10;


figure()
hold on;
plot(t8, dxcore89, 'r-');
plot(t8, dxcore910, 'b-');
legend('xcore8-xcore9', 'xcore9-xcore10');
pbaspect([1 1 1])
hold off;

dxcore910 = 4*dxcore910;


figure()
plot(t8, dxcore89, 'r-');
hold on;
plot(t8, dxcore910, 'b-');
legend('xcore8-xcore9', 'xcore9-xcore10');
pbaspect([1 1 1])
hold off;
