close all; clear all; clc

figure(1)
hold on
%plot points 
%point A
plot3(1,2,3,'Color','r','marker','^','DisplayName','A','MarkerFaceColor','r')
%point B
plot3(-3,2,5,'Color','b','marker','s','DisplayName','B','MarkerFaceColor','b')
%point c
plot3(pi,exp(1),-sqrt(2),'Color','m','marker','o','DisplayName','C','MarkerFaceColor','m')

x = [10 -10 -10 10]; % Generate data for x vertices
y = [10 10 -10 -10]; % Generate data for y vertices;
z = (3.90335731207141490E-002.*x + 0.36338249407750045.*y - 1)/(-7.80671462414282979E-002);
patch(x,y,z,'g','DisplayName','Plane')
legend
hold off