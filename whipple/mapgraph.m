function [x, y] = mapgraph(str)
%
% MAPGRAPH
% Written by R.W. Lessard
% Purdue University
% 
% Arguments
% str = 'log'    for log plot
%     = 'linear' for linear plot
%     = 'loglin' for log x, lin y
%     = 'linlog' for lin x, log y
% output [x,y]
%
% Use right button for last point.
%
clf
hold on
disp('click left button on top left corner')
[xtl,ytl,but] = ginput(1)
xleft = input('x = ');
ytop = input('y = ');
disp('click left button on bottom right corner')
[xbr,ybr,but] = ginput(1)
xright = input('x = ');
ybottom = input('y = ');
set(gca,'Position',[0.13 0.11 xbr*.775 ytl*.815]);
if strcmp(str,'log') == 1
   set(gca,'YScale','log');
   set(gca,'XScale','log');
end
if strcmp(str,'loglin') == 1
   set(gca,'XScale','log');
end
if strcmp(str,'linlog') == 1
   set(gca,'YScale','log');
end
axis([xleft xright ybottom ytop]);
but = 1;
x = [];
y = [];
n = 0;
while but == 1
  [xi,yi,but] = ginput(1);
  plot(xi,yi,'go')
  n = n+1;
  x(n,1) = xi;
  y(n,1) = yi;
end
