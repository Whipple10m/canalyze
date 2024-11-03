function y = mdiagnose(file)
%MDIAGNOSE   Performs system diagnositics
%   MDIAGNOSE('file')
%   The first argument is the file name of the output written by
%   fz2hdf [gtxxxxxxd.mat]
minfo = zeros(1,2);
azi = zeros(1,1);
eval(['load ' file]);
switch minfo(2)
case 120
   npmt = 109;
   fov = 3.5;
case 156
   npmt = 151;
   fov = 4;
case 336
   fov = 5;
   npmt = 331;
case 492
   npmt = 490;
   fov = 4;
otherwise
   npmt = 0;
end
%
% Begin Plotting
%
% Positioning
%
%
% Page 1
%
figure
set(gcf,'Position',[10 10 523 638])
set(gcf,'PaperPosition',[0.5 1.0 7.5 9.0])
%
% TOP TEXT
%
subplot('position',[0.10 0.975 0.90 0.05])
axis('off');
[status,whoami] = unix('whoami');
[status,whereami] = unix('hostname');
[status,date] = unix('date');
text(0.0,0.0,[whoami(1:size(whoami,2)-1) '@' ...
whereami(1:size(whereami,2)-1) '  ' date(1:size(date,2)-1)])
%
% RATE
%
[xx,yy] = stairs(0.5:1.0:size(rate,1)-0.5,rate./60);
subplot('position',[0.10 0.80 0.35 0.15]),plot(xx,yy,'k-');
axis([0 28 0 max(rate./60)+max(rate./60)/5])
xlabel('time (min)')
ylabel('raw rate [Hz]')
%
% DELTA T
%
[xx,yy] = stairs(0.5:1.0:499.5,gpsdt);
subplot('position',[0.10 0.55 0.35 0.20]),semilogx(xx,yy,'g-')
hold on
[xx,yy] = stairs(0.5:1.0:499.5,oscdt);
semilogx(xx,yy,'r-');
axis([0 500 0 max([gpsdt ; oscdt])+max([gpsdt ; oscdt])/20])
legend('GPS','OSC')
set(gca,'XTick',[1 10 100])
grid
hold off
xlabel('\Delta t [ms]')
ylabel('number of events')
%
% GPS-OSC
%
[xx,yy] = stairs(0.5:1:499.5,gpsoscdt(500:-1:1));
subplot('position',[0.55 0.55 0.35 0.20]),semilogx(xx,yy,'r-');
[xx,yy] = stairs(0.5:1:499.5,gpsoscdt(501:1000));
hold on
semilogx(xx,yy,'g-');
axis([0 500 0 max(gpsoscdt) + max(gpsoscdt)/20]);
xlabel('GPS - OSC time difference [ms]');
ylabel('number of events');
legend('neg','pos');
set(gca,'XTick',[1 10 100])
grid
hold off
%
% RUN/DIAGNOSTICS INFO
%
subplot('position',[0.10 0.05 0.35 0.40])
axis('off');
h = text(0.0,1.0,'Run Info');
set(h,'FontWeight','Bold');
text(0.0,0.95,['File: ' file]);
string = sprintf('Number ADCs: %3d',minfo(2));
text(0.0,0.90,string);
string = sprintf('Duration: %2d mins',minfo(1));
text(0.0,0.85,string);
string = sprintf('Live Time: %4.1f [mins]',minfo(4));
text(0.0,0.80,string);
string = sprintf('UTC Start: %9.3f [mjd]',minfo(3));
text(0.0,0.75,string);
text(0.0,0.70,['Source: ' source(1:7)]);
string = sprintf('RA:  %8.1f [hhmmss]',radtohhmmss(minfo(5)));
text(0.0,0.65,string);
string = sprintf('DEC: %8.1f [ddmmss]',radtoddmmss(minfo(6)));
text(0.0,0.60,string);
% Timing Faults
h = text(0.0,0.55,'Timing Report');
set(h,'FontWeight','Bold');
string = 'Oscillator calibrated: ';
if mfault(1) == 0.0
  string = [string 'no'];
else
  string = [string 'yes'];
end
text(0.0,0.50,string);
if npmt > 151
   string = sprintf('Average Osc. Frequency: %5.1f [ns/ns]',mfault(2)/1.0e9);
   text(0.0,0.45,string);
   string = 'GPS clock fault: ';
   if mfault(3) == 0.0
      string = [string 'no'];
   else
      string = [string 'yes'];
   end
   text(0.0,0.40,string);
else
   string = sprintf('Average Osc. Frequency: %5.1f [MHz]',mfault(2)/1000000.0);
   text(0.0,0.45,string);
   text(0.0,0.40,sprintf('Number of missed markers: %d',round(mfault(3))));
end
text(0.0,0.35,sprintf('Number of GPS errors: %d',round(mfault(4))));
text(0.0,0.30,sprintf('Number of OSC errors: %d',round(mfault(5))));
text(0.0,0.25,sprintf('Number of GPS/OSC differences: %d',round(mfault(6))));
% TUBES OFF
tubesoff = find(hvstat(1:npmt) == 0.0);
string = [];
h = text(0.0,0.20,'HV Report');
set(h,'FontWeight','Bold');
text(0.0,0.15,'Tubes off');
row = 0.10;
for i=1:size(tubesoff,1)
   string = [string sprintf('%4d',tubesoff(i))];
   if mod(i,7) == 0
      text(0.0,row,string);
      row = row - 0.05;
      string = [];
   end
end
if ~isempty(string)
   text(0.0,row,string);
end
%
% TRACKING ERROR
%
[xx,yy] = stairs(0.5:1.0:89.5,elev);
subplot('position',[0.55 0.35 0.35 0.10]),plot(xx,yy,'k-');
if max(elev) == 0
   axis([0 90 -1 1]);
else
   axis([0 90 0 max(elev)+max(elev)/10]);
end
xlabel('elevation [deg]');
ylabel('number of rcds');
[xx,yy] = stairs(-179.5:1.0:179.5,azi);
subplot('position',[0.55 0.20 0.35 0.10]),plot(xx,yy,'k-');
if max(azi) == 0
   axis([-180 180 -1 1]);
else
   axis([-180 180 0 max(azi)+max(azi)/10]);
end
xlabel('azimuth [deg]');
ylabel('number of rcds');
[xx,yy] = stairs(0.01:.02:0.99,trackerr);
subplot('position',[0.55 0.05 0.35 0.10]),plot(xx,yy,'k-');
xlabel('tracking error [deg]');
ylabel('number of rcds');
%
% FFT
%
subplot('Position',[0.55 0.80 0.35 0.15])
axis('off')
z = fft(rawtime);
z = z.*conj(z)/size(z,1);
f = (524288/1680.0*(0:262143)/524288)';
plot(f(2:size(f,1)),z(2:size(z,1)/2));
axis([0 f(size(f,1)) 0 max(z(2:size(z,1)/2))+max(z(2:size(z,1)/2))/10]);
xlabel 'frequency [Hz]'
ylabel 'Fourier power'
% 
% Page 2
%
figure
%
%Begin Plotting
%
% Positioning
%
set(gcf,'Position',[10 10 523 638])
set(gcf,'PaperPosition',[0.5 1.0 7.5 9.0])
%
% TOP TEXT
%
subplot('position',[0.10 0.975 0.90 0.05])
axis('off');
[status,whoami] = unix('whoami');
[status,whereami] = unix('hostname');
[status,date] = unix('date');
text(0.0,0.0,[whoami(1:size(whoami,2)-1) '@' ...
whereami(1:size(whereami,2)-1) '  ' date(1:size(date,2)-1)])
%
% Peds and pedvars
%
eval(['load whipplecoords_' num2str(npmt) 'RWL.mat']);
[iyy,imm,idd,fr,j] = sla_djcl(minfo(3));
if iyy < 2000
   minfo(7) = (iyy-1900)*10000+imm*100+idd;
else
   minfo(7) = (iyy-2000)*10000+imm*100+idd;
end
minfo(8) = sla_gmst(minfo(3));
p = getcpeds(file(1:8),minfo(7),'./');
%
% Plot peds
%
subplot(2,2,1);
peds = p(1:npmt);
maxp = max(peds);
hold on
for i=1:npmt
  if i > 379
    circle(whipplecoords(i,2),whipplecoords(i,3),...
	peds(i)/maxp*.125/2,[.8 .8 .8]);
  else
    circle(whipplecoords(i,2),whipplecoords(i,3),...
	peds(i)/maxp*.125/2,[.8 .8 .8]);
  end
end
xlabel 'x [deg]'
ylabel 'y [deg]'
axis([-fov/2 fov/2 -fov/2 fov/2]);
axis 'square'
title 'pedestals'
%
% Plot pedvars
%
subplot(2,2,2);
pedvars = p(minfo(2)+1:minfo(2)+1+npmt);
%pedvars = pedvars - mean(pedvars);
%pedvars(find(pedvars<0)) = 0;
maxp = max(pedvars);
hold on
for i=1:npmt
  if i > 379
    circle(whipplecoords(i,2),whipplecoords(i,3),...
          pedvars(i)/maxp*.125/2,[.8 .8 .8]);
  else
    circle(whipplecoords(i,2),whipplecoords(i,3),...
          pedvars(i)/maxp*.125/2,[.8 .8 .8]);
  end
end
xlabel 'x [deg]'
ylabel 'y [deg]'
axis([-fov/2 fov/2 -fov/2 fov/2]);
axis 'square'
title 'pedestal variances'
%
% Plot Stars
%
daysinmonth = [31 28 31 30 31 30 31 31 30 31 30 31];
dayofyear = sum(daysinmonth(1:imm-1)) + idd;
epoch = iyy + dayofyear/365;
whiplong = 1.935190;
s = getstars(minfo(5),minfo(6),epoch,minfo(8)-whiplong,6,3.0);
plotstars(s);
s = getstars(minfo(5),minfo(6),epoch,minfo(8)+(28/60*pi/12)-whiplong,6,2.0);
plotstars(s);
