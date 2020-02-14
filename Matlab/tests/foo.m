function [ ] = runTests()
%TEST_SOLVER Summary of this function goes here
%   Detailed explanation goes here

%for testing
t=[0 1];
y0 = [1 1];


% [T,Ym]=ode15s(@odefun,t,y0);
% Ym = (Ym)';

% s=dsolve('Dy1 = -10*y1','Dy2 = -12*y2','y1(0)=2','y2(0)=2');
% y2an = subs(s.y2,T);
% y1an = subs(s.y1,T);
% 
% hold on
% plot(T,y2an,'sb-')
% plot(T,y1an,'^b-')
% plot(T,Y(1,:),'og-')
% plot(T,Y(2,:),'.g-')


end

function dy = odefun(t,y)

y = y(:);

dy(1) = -10*y(1);
dy(2) = -12*y(2)*y(1);

dy = dy(:);

end



