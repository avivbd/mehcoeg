function compare_with_python

[T,Y] = ode45(@odefun,[0 50],[1 2]);

subplot(2,1,1)
plot(T,Y(:,1))

subplot(2,1,2)
plot(T,Y(:,2),'g')


function dydt = odefun(t,y)
dydt = zeros(2,1);
dydt(1) = -2*y(1)+2*y(2) + 1 + 50*randn(1);
dydt(2) = -4*y(1) + 2*y(2) + 50*randn(1); 

