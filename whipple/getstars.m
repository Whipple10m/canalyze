function s = getstars(raz,decz,epoch,sidtime,maxmag,fov)
% GETSTARS Return array of stars in field of view.
%    GETSTARS(RAZ,DECZ,EPOCH,SIDTIME,MAG,FOV) returns an array of
%    star coordinates [deg], magnititudes and colors, given telescope
%    tracking coordinates [rad], sidereal time [rad,local], epoch, 
%    limiting magnitude and field of view [deg].
%
whiplat = 0.552978;
whiplong = 1.935190;
fp = fopen('yalebsc.dat','r');
[A,count] = fscanf(fp,'%f %d %f %f',[4 inf]);
ra = hhmmsstorad(A(1,:));
dec = ddmmsstorad(A(2,:));
mag = A(3,:);
b_v = A(4,:);
clear A;
%
[ra,dec] = precessfrom2000(ra,dec,epoch);
%
hourangle = sidtime - raz;
[azz, altz] = sla_de2h(hourangle,decz,whiplat);
%
hourangle = sidtime - ra;
[az, alt] = sla_de2h(hourangle,dec,whiplat);
%
[x, y, j] = sla_ds2tp(az,alt,azz,altz);
%
x = x*180/pi;
y = y*180/pi;
%
starsinfov = find(((x.^2 + y.^2) < fov^2) & (j == 0) & (mag <= maxmag));
s = [x([starsinfov])' y([starsinfov])' mag([starsinfov])' b_v([starsinfov])'];
