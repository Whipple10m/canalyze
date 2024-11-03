function y = rectan(cx,cy,w,h,ang)
%
ang= ang*2*pi/360.;
if ang == 0

	func = sprintf('%f',cy+(h/2));
	[x,y] = fplot(func,[cx-w/2 cx+w/2]);
	plot(x,y,'r');
	func = sprintf('%f',cy-(h/2));
	[x,y] = fplot(func,[cx-w/2 cx+w/2]);
	plot(x,y,'r');
	y=(cy-h/2):h/10:(cy+h/2);
	x=zeros(size(y));
	x=x+(cx-w/2);
	plot(x,y,'r');
	x=zeros(size(y));
	x=x+(cx+w/2);
	plot(x,y,'r');
else
	r=sqrt(h^2/4 + w^2/4);
	theta = atan(h/w);
	xur = cx+r*cos(theta-ang);
	yur = cy+r*sin(theta-ang);
	xul = cx-r*cos(theta+ang);
	yul = cy+r*sin(theta+ang);
	m = (yul - yur)/(xul - xur);
	b1 = yul-(m*xul);
	if b1 < 0
		func = sprintf('%f*x - %f',m,-b1);
	else
		func = sprintf('%f*x + %f',m,b1);
	end
	[x,y] = fplot(func,[xul xur]);
	plot(x,y,'r');
	b2 = yur-(m*xur)-(h/sin((pi/2)-ang));
	if b2 < 0
		func = sprintf('%f*x - %f',m,-b2);
	else
		func = sprintf('%f*x + %f',m,b2);
	end
	xlr = cx+r*cos(theta+ang);
	ylr = cy-r*sin(theta+ang);
	xll = cx-r*cos(theta-ang);
	yll = cy-r*sin(theta-ang);
	[x,y] = fplot(func,[xll xlr]);
	plot(x,y,'r');
	m = (yul-yll)/(xul-xll);
	b1 = yul-(m*xul);
	if b1 < 0
		func = sprintf('%f*x - %f',m,-b1);
	else
		func = sprintf('%f*x + %f',m,b1);
	end
	[x,y] = fplot(func,[xll xul]);
	plot(x,y,'r');
	b2 = yur-(m*xur);
	if b2 < 0
		func = sprintf('%f*x - %f',m,-b2);
	else
		func = sprintf('%f*x + %f',m,b2);
	end
	[x,y] = fplot(func,[xlr xur]);
	plot(x,y,'r');
end
