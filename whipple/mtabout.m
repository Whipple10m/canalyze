function [on , off] = mtabout(files,varargin)
%MTABOUT   Tabulates supercut results.
%   [ON,OFF] = MTABOUT(file_list,'PropertyName','ProperyValue')
%   The first argument is the file name of the files list and is required
%   The next arguments are optional Properties as list below.
%
%   [Options]
%   option 'Grid', values ['RaDec','Camera'] - default 'RaDec'
%   option 'Epoch', - default 2000.0
%   option 'Mode', values ['Pair','Tracking'] - default 'Pair'
%   option 'Plot', values ['Image','Spectrum','LightCurve','Phase']
%                  - default 'Image'
%   option 'AlphaRatio',  - default flat
%   option 'ThetaRatio', - default flat
%
%   [Output]
%   ON/OFF is a structure of results: lightcurve, alfa, theta, image, phase
%      lightcurve is matrix of excess rates
%        column 1: raw
%        column 2: trigger
%        column 3: shape
%        column 4: alfa ** super cuts **
%        column 5: theta
%        column 6: tracking alfa
%        column 7: tracking theta
%      alfa,theta,image,phase are vectors of counts
%      time,solar_bary,orbit_bary,phase,energy are accumlated event parameters
%      ChiSq 
%
disp '**new** 980712 mtabout uses mcut alpha cut value for tracking analysis'
disp '**new** 990127 mtabout returns structure of results'
disp '**new** 990504 mtabout calculates ChiSq for Phase Option'
disp '**new** 000512 mtabout canalyze version 4.0 only'
warning off
%
% default options
%
grid = 'RaDec';
epoch = 2000.0;
mode = 'Pair';
plotopt = 'Image';
%
% initialize
%
alfaratio = 0; thetaratio = 0;
theta = 0; phase = 0;
onchisq = 0; offchisq = 0;
%
% options
%
for i=1:2:nargin-2
   option=varargin(i);
   optionval=varargin(i+1);
   if strcmp(option,'Grid')
      if strcmp(optionval,'RaDec')
         grid='RaDec'
      elseif strcmp(optionval,'Camera')
         grid='Camera'
      end
   elseif strcmp(option,'Epoch')
      if isstr(optionval)==0
         epoch=optionval{1}
      end
   elseif strcmp(option,'Mode')
      if strcmp(optionval,'Pair')
         mode='Pair'
      elseif strcmp(optionval,'Tracking')
         mode='Tracking'
         plotopt='LightCurve'
      end
   elseif strcmp(option,'AlphaRatio')
      if isstr(optionval)==0
         alfaratio=optionval{1}
      end
   elseif strcmp(option,'ThetaRatio')
      if isstr(optionval)==0
         thetaratio=optionval{1}
      end
   elseif strcmp(option,'Plot')
      if strcmp(optionval,'LightCurve')
         plotopt='LightCurve'
      elseif strcmp(optionval,'Spectrum')
         plotopt='Spectrum'
      elseif strcmp(optionval,'Image');
         plotopt='Image'
      elseif strcmp(optionval,'Phase');
         plotopt='Phase'
      end
   end
end 
%
% load file list
%
fp = fopen(files,'r');
files = fscanf(fp,'%s');
fclose(fp);
%
% open ascii output file
%
fp = fopen('mtabout.txt','w');
fprintf(fp,...
        '-------------------------------------------------\n');
fprintf(fp,...
        '|            |   RAW|  TRIG|    SH| ALPHA| THETA|\n');
fprintf(fp,...
        '-------------------------------------------------\n');
%
% read in files
%
nfiles = size(files,2)/12;
maxsig = 5;
i=1;
while i <= nfiles
   if i==1
%
% ON
%
      eval(['load ' files(1:12)]); onfile = files(1:12);
%
% info 
%
      source_ra = info.source_ra; source_dec = info.source_dec;
      source_date = info.utc_start;
      onduration = info.duration; tmponduration = info.duration;
%
% events passed
%
      onpassed = passed; tmponpassed = passed;
%
% distributions
%
      onalfa = alfa; tmponalfa = alfa;
      ontheta = theta; tmpontheta = theta;
      onimage = smoothed_image';
      onphased = phased; 
%
% photon data
%
      ontime = time; onsolar_bary = solar_bary;
      onorbit_bary = orbit_bary; 
      onphase = phase;
      onenergy = energy;
      if strcmp(mode,'Pair')
%
% OFF
%
         eval(['load ' files(13:24)]); offfile = files(13:24);
%
% info
%
         offduration = info.duration; tmpoffduration = info.duration;
%
% events passed
%
         offpassed = passed; tmpoffpassed = passed;
%
% distributions
%
         offalfa = alfa; tmpoffalfa = alfa;
         offtheta = theta; tmpofftheta = theta;
         offimage = smoothed_image';
         offphased = phased; 
%
% photon data
%
         offtime = time; offsolar_bary = solar_bary;
         offorbit_bary = orbit_bary; 
         offphase = phase;
         offenergy = energy;
         i=i+1;
      end
%
% determine bin sizes and distribution x axes
%
      phasebinwidth=2/size(phased,1);
      phasex = phasebinwidth/2:phasebinwidth:2-phasebinwidth/2;
%
      alfabinwidth=90/size(alfa,1);
% tracking on region 0 - alfa cut
      upperalfaon = cut.alfa/alfabinwidth;
% tracking off region 20 - 65
      loweralfaoff = 20/alfabinwidth + 1;
      upperalfaoff = 65/alfabinwidth;
% assume flat alfa as default
      if alfaratio == 0
         alfaratio = upperalfaon/(upperalfaoff - loweralfaoff + 1);
      end
      alfax = alfabinwidth/2:alfabinwidth:90-alfabinwidth/2;
%
      thetabinwidth=2/size(theta,1);
% tracking on region 0 - theta cut
      upperthetaon = round(cut.theta/thetabinwidth);
% tracking off region 0.50 - 1.50 deg
      lowerthetaoff = round(0.50/thetabinwidth + 1);
      upperthetaoff = round(1.50/thetabinwidth);
      thetax = thetabinwidth/2:thetabinwidth:2-thetabinwidth/2;
      areax = pi*([thetabinwidth:thetabinwidth:2.0].^2 - ...
                  [0.0:thetabinwidth:2-thetabinwidth].^2)';
% assume flat theta^2 - equal area
      if thetaratio == 0
         thetaratio = pi*(thetabinwidth*upperthetaon)^2/...
                      (pi*(thetabinwidth*upperthetaoff)^2 - ...
                       pi*(thetabinwidth*lowerthetaoff)^2);
      end
      i=i+1;
   else
%
% ON - file > 1
%
      eval(['load ' files((i-1)*12+1:i*12)]); onfile = files((i-1)*12+1:i*12);
%
% info 
%
      onduration = onduration + info.duration; tmponduration = info.duration;
%
% events passed
%
      onpassed = onpassed + passed; tmponpassed = passed;
%
% distributions
%
      onalfa = onalfa + alfa; tmponalfa = alfa;
      ontheta = ontheta + theta; tmpontheta = theta;
      onimage = onimage + smoothed_image';
      onphased = onphased + phased; 
      ontime = [ontime ; time]; 
      onsolar_bary = [onsolar_bary ; solar_bary];
      onorbit_bary = [onorbit_bary ; orbit_bary]; 
      onphase = [onphase ; phase]; 
      onenergy = [onenergy ; energy];
%
% OFF - file > 1
%
      if strcmp(mode,'Pair')
         eval(['load ' files(i*12+1:(i+1)*12)]);
         offfile = files(i*12+1:(i+1)*12);
         offduration = offduration + info.duration; 
	 tmpoffduration = info.duration;
         offpassed = offpassed + passed; tmpoffpassed = passed;
         offalfa = offalfa + alfa; tmpoffalfa = alfa;
         offtheta = offtheta + theta; tmpofftheta = theta;
         offimage = offimage + smoothed_image';
         offphased = offphased + phased;
         offtime = [offtime ; time]; 
	 offsolar_bary = [offsolar_bary ; solar_bary];
         offorbit_bary = [offorbit_bary ; orbit_bary]; 
         offphase = [offphase ; phase]; 
         offenergy = [offenergy ; energy];
         i=i+1;
      end
      i=i+1;
   end
%
% intermediate calculations
%
   if strcmp(mode,'Pair')
      timeratio = tmponduration/tmpoffduration;
      excess = tmponpassed - timeratio*tmpoffpassed;
      sigma = excess./sqrt(timeratio*(tmponpassed + tmpoffpassed));
      rate = excess/tmponduration;
      erate = sqrt(tmponpassed + timeratio*timeratio*tmpoffpassed)/...
              tmponduration;
   end
%
% tracking analysis
%
   tmponpassed(6:7) = [sum(tmponalfa(1:upperalfaon)) ...
                        sum(tmpontheta(1:upperthetaon))];
   tmpoffpassed(6:7) = [sum(tmponalfa(loweralfaoff:upperalfaoff)) ...
                         sum(tmpontheta(lowerthetaoff:upperthetaoff))];
   excess(6) = tmponpassed(6) - alfaratio*tmpoffpassed(6);
   excess(7) = tmponpassed(7) - thetaratio*tmpoffpassed(7);
   sigma(6) = excess(6)/sqrt(alfaratio*(tmponpassed(6) + tmpoffpassed(6)));
   sigma(7) = excess(7)/sqrt(thetaratio*(tmponpassed(7) + tmpoffpassed(7)));
   rate(6:7) = excess(6:7)/tmponduration;
   erate(6) = sqrt(tmponpassed(6) + ...
              alfaratio*alfaratio*tmpoffpassed(6))/tmponduration;
   erate(7) = sqrt(tmponpassed(7) + ...
              thetaratio*thetaratio*tmpoffpassed(7))/tmponduration;
%
% output intermediate results to mtabout.txt
%
   if strcmp(mode,'Pair')
      fprintf(fp,'|%12s|%6d|%6d|%6d|%6d|%6d|\n',onfile,...
              tmponpassed(1),tmponpassed(2),tmponpassed(3),...
              tmponpassed(4),tmponpassed(5));
      fprintf(fp,'|%12s|%6d|%6d|%6d|%6d|%6d|\n',offfile,...
              round(timeratio*tmpoffpassed(1)),...
              round(timeratio*tmpoffpassed(2)),...
              round(timeratio*tmpoffpassed(3)),...
              round(timeratio*tmpoffpassed(4)),...
              round(timeratio*tmpoffpassed(5)));
      fprintf(fp,'|%12s|%6.2f|%6.2f|%6.2f|%6.2f|%6.2f|\n',...
              'excess',...
              sigma(1),sigma(2),sigma(3),...
              sigma(4),sigma(5));
      fprintf(fp,'|%12s|%6.2f|%6.2f|%6.2f|%6.2f|%6.2f|\n',...
              'rate',...
              rate(1),rate(2),rate(3),...
              rate(4),rate(5));
      fprintf(fp,'%4s %12.6f %4s %6.2f/%6.2f %5s %6.2f %6.2f \n',...
	'UTC:',info.utc_start,...
        'DUR:',tmponduration,tmpoffduration,...
        'RATE:',rate(5),erate(5));
%
% build light curve
%
      lightcurve((i-1)/2,1) = info.utc_start;
      lightcurve((i-1)/2,2:8) = rate';
      lightcurve((i-1)/2,9:15) = erate';
   else
      fprintf(fp,'|%12s|%6d|%6d|%6d|%6d|%6d|\n',onfile,...
              tmponpassed(1),tmponpassed(2),tmponpassed(3),...
              tmponpassed(6),tmponpassed(7));
      fprintf(fp,'|%12s|%6s|%6s|%6s|%6d|%6d|\n','',...
              '','','',round(alfaratio*tmpoffpassed(6)),...
              round(thetaratio*tmpoffpassed(7)));
      fprintf(fp,'|%12s|%6s|%6s|%6s|%6.2f|%6.2f|\n','excess',...
              '','','',sigma(6),sigma(7));
      fprintf(fp,'|%12s|%6s|%6s|%6s|%6.2f|%6.2f|\n','rate',...
              '','','',rate(6),rate(7));
      fprintf(fp,'%4s %12.6f %4s %6.2f %5s %6.2f %5.2f \n',...
	'UTC:',info.utc_start,...
        'DUR:',tmponduration,...
        'RATE:',rate(6),erate(6));
      lightcurve(i-1,1) = info.utc_start;
      lightcurve(i-1,2:8) = rate;
      lightcurve(i-1,9:15) = erate;
   end
end
fprintf(fp,...
        '-------------------------------------------------\n');
%
% final calculations
%
if strcmp(mode,'Pair')
   timeratio = onduration/offduration;
   excess = onpassed - timeratio*offpassed;
   sigma = excess./sqrt(timeratio*(onpassed + offpassed));
end
%
% tracking analysis
%
onpassed(6:7) = [sum(onalfa(1:upperalfaon)) sum(ontheta(1:upperthetaon))];
offpassed(6:7) = [sum(onalfa(loweralfaoff:upperalfaoff)) ...
                   sum(ontheta(lowerthetaoff:upperthetaoff))];
excess(6) = onpassed(6) - alfaratio*offpassed(6);
excess(7) = onpassed(7) - thetaratio*offpassed(7);
sigma(6) = excess(6)/sqrt(alfaratio*(onpassed(6) + offpassed(6)));
sigma(7) = excess(7)/sqrt(thetaratio*(onpassed(7) + offpassed(7)));
rate = excess/onduration;
if strcmp(mode,'Pair')
   fprintf(fp,'|%12s|%6d|%6d|%6d|%6d|%6d|\n','Total on',...
           onpassed(1),onpassed(2),onpassed(3),...
           onpassed(4),onpassed(5));
   fprintf(fp,'|%12s|%6d|%6d|%6d|%6d|%6d|\n','Total off',...
           round(timeratio*offpassed(1)),...
           round(timeratio*offpassed(2)),...
           round(timeratio*offpassed(3)),...
           round(timeratio*offpassed(4)),...
           round(timeratio*offpassed(5)));
   fprintf(fp,'|%12s|%6.2f|%6.2f|%6.2f|%6.2f|%6.2f|\n',...
           'excess',...
           sigma(1),sigma(2),sigma(3),...
           sigma(4),sigma(5));
   fprintf(fp,'|%12s|%6.2f|%6.2f|%6.2f|%6.2f|%6.2f|\n',...
           'rate',...
           rate(1),rate(2),rate(3),...
           rate(4),rate(5));
else
   fprintf(fp,'|%12s|%6d|%6d|%6d|%6d|%6d|\n','Total on',...
           onpassed(1),onpassed(2),onpassed(3),...
           onpassed(6),onpassed(7));
   fprintf(fp,'|%12s|%6s|%6s|%6s|%6d|%6d|\n','Total off',...
           '','','',round(alfaratio*offpassed(6)),...
           round(thetaratio*offpassed(7)));
   fprintf(fp,'|%12s|%6s|%6s|%6s|%6.2f|%6.2f|\n','excess',...
           '','','',sigma(6),sigma(7));
   fprintf(fp,'|%12s|%6s|%6s|%6s|%6.2f|%6.2f|\n','rate',...
           '','','',rate(6),rate(7));
end
fprintf(fp,...
        '-------------------------------------------------\n');
fclose(fp);
if strcmp(mode,'Pair')
   imsigma = (onimage-timeratio*offimage)./sqrt(timeratio*(onimage+offimage));
   imsigma(find(isnan(imsigma))) = zeros(1,size(find(isnan(imsigma)),1));
   imsigmac = min(min(imsigma)):0.1:max([max(imsigma) maxsig]);
   z = [imsigmac' imsigmac'];
else
   imsigma = onimage;
end
imsigmax = -1.9:.1:1.9;
imsigmay = -1.9:.1:1.9;
if strcmp(mode,'Pair')
  [maxsigma,maxrow] = max(imsigma);
  [maxsigma,maxcol] = max(maxsigma);
   sprintf('Peak: SIG =  %5.2f, ON = %5d, OFF = %5d at x=%5.2f, y=%5.2f',...
   maxsigma,onimage(maxrow(maxcol),maxcol),offimage(maxrow(maxcol),maxcol),...
   imsigmax(maxcol)+0.05,imsigmay(maxrow(maxcol))+0.05)
end
%
%Begin Plotting
%
% Positioning
%
% position - [left bottom width height]
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
% SUMMARY TEXT
%
subplot('Position',[0.5 0.70 0.35 0.25])
axis('off')
text(0.0,1.0,['Source: ' source])
if strcmp(mode,'Pair')
   string = sprintf('Duration (on/off): %6.1f/%6.1f min',...
      onduration,offduration);
else
   string = sprintf('Duration (on): %6.1f min',onduration);
end
text(0.0,0.92,string)
text(0.0,0.84,'Tracking Analysis');
string = sprintf('Significance: %6.2f (alpha)',sigma(6));
text(0.05,0.76,string)
string = sprintf('Rate: %6.2f min^{-1} (alpha)',rate(6));
text(0.05,0.68,string)
string = sprintf('Significance: %6.2f (theta)',sigma(7));
text(0.05,0.60,string)
string = sprintf('Rate: %6.2f min^{-1} (theta)',rate(7));
text(0.05,0.52,string)
string = sprintf('Flux: %6.2g cm^{-2} s^{-1}',rate(6)/3.5e8/60.0);
text(0.05,0.44,string)
string = sprintf('Flux u.l.: %6.2g cm^{-2} s^{-1}',...
   gupper(onpassed(6),offpassed(6),alfaratio)/3.5e8/onduration/60.0);
text(0.05,0.36,string)
if strcmp(mode,'Pair')
   text(0.0,0.28,'ON/OFF Analysis');
   string = sprintf('Significance: %6.2f (alpha)',sigma(4));
   text(0.05,0.20,string)
   string = sprintf('Rate: %6.2f min^{-1} (alpha)',rate(4));
   text(0.05,0.12,string)
   string = sprintf('Significance: %6.2f (theta)',sigma(5));
   text(0.05,0.04,string)
   string = sprintf('Rate: %6.2f min^{-1} (theta)',rate(5));
   text(0.05,-0.04,string)
   string = sprintf('Flux: %6.2g cm^{-2} s^{-1}',rate(4)/3.5e8/60.0);
   text(0.05,-0.12,string)
   string = sprintf('Flux u.l.: %6.2g cm^{-2} s^{-1}',...
      gupper(onpassed(4),offpassed(4),timeratio)/3.5e8/onduration/60.0);
   text(0.05,-0.20,string)
end
%
% PASSED TEXT
%
subplot('Position',[0.10 0.17 0.90 0.10])
%plot([])
axis('off')
text(0.20,1.0,'RAW','HorizontalAlignment','right')
text(0.35,1.0,'TRIG','HorizontalAlignment','right')
text(0.50,1.0,'SHAPE','HorizontalAlignment','right')
text(0.65,1.0,'ALPHA','HorizontalAlignment','right')
text(0.80,1.0,'THETA','HorizontalAlignment','right')
text(0.10,0.8,'ON','HorizontalAlignment','right')
string = sprintf('%8d',onpassed(1));
text(0.20,0.8,string,'HorizontalAlignment','right')
string = sprintf('%8d',onpassed(2));
text(0.35,0.8,string,'HorizontalAlignment','right')
string = sprintf('%8d',onpassed(3));
text(0.50,0.8,string,'HorizontalAlignment','right')
if strcmp(mode,'Pair')
   string = sprintf('%8d',onpassed(4));
   text(0.65,0.8,string,'HorizontalAlignment','right')
   string = sprintf('%8d',onpassed(5));
   text(0.80,0.8,string,'HorizontalAlignment','right')
else
   string = sprintf('%8d',onpassed(6));
   text(0.65,0.8,string,'HorizontalAlignment','right')
   string = sprintf('%8d',onpassed(7));
   text(0.80,0.8,string,'HorizontalAlignment','right')
end
text(0.10,0.6,'OFF','HorizontalAlignment','right')
if strcmp(mode,'Pair')
   string = sprintf('%8.0f',timeratio*offpassed(1));
   text(0.20,0.6,string,'HorizontalAlignment','right')
   string = sprintf('%8.0f',timeratio*offpassed(2));
   text(0.35,0.6,string,'HorizontalAlignment','right')
   string = sprintf('%8.0f',timeratio*offpassed(3));
   text(0.50,0.6,string,'HorizontalAlignment','right')
   string = sprintf('%8.0f',timeratio*offpassed(4));
   text(0.65,0.6,string,'HorizontalAlignment','right')
   string = sprintf('%8.0f',timeratio*offpassed(5));
   text(0.80,0.6,string,'HorizontalAlignment','right')
else
   string = sprintf('%8.0f',alfaratio*offpassed(6));
   text(0.65,0.6,string,'HorizontalAlignment','right')
   string = sprintf('%8.0f',thetaratio*offpassed(7));
   text(0.80,0.6,string,'HorizontalAlignment','right')
end
text(0.10,0.4,'EXCESS','HorizontalAlignment','right')
if strcmp(mode,'Pair')
   string = sprintf('%8.0f',excess(1));
   text(0.20,0.4,string,'HorizontalAlignment','right')
   string = sprintf('%8.0f',excess(2));
   text(0.35,0.4,string,'HorizontalAlignment','right')
   string = sprintf('%8.0f',excess(3));
   text(0.50,0.4,string,'HorizontalAlignment','right')
   string = sprintf('%8.0f',excess(4));
   text(0.65,0.4,string,'HorizontalAlignment','right')
   string = sprintf('%8.0f',excess(5));
   text(0.80,0.4,string,'HorizontalAlignment','right')
else
   string = sprintf('%8.0f',excess(6));
   text(0.65,0.4,string,'HorizontalAlignment','right')
   string = sprintf('%8.0f',excess(7));
   text(0.80,0.4,string,'HorizontalAlignment','right')
end
text(0.10,0.2,'SIGMA','HorizontalAlignment','right')
if strcmp(mode,'Pair')
   string = sprintf('%8.2f',sigma(1));
   text(0.20,0.2,string,'HorizontalAlignment','right')
   string = sprintf('%8.2f',sigma(2));
   text(0.35,0.2,string,'HorizontalAlignment','right')
   string = sprintf('%8.2f',sigma(3));
   text(0.50,0.2,string,'HorizontalAlignment','right')
   string = sprintf('%8.2f',sigma(4));
   text(0.65,0.2,string,'HorizontalAlignment','right')
   string = sprintf('%8.2f',sigma(5));
   text(0.80,0.2,string,'HorizontalAlignment','right')
else
   string = sprintf('%8.2f',sigma(6));
   text(0.65,0.2,string,'HorizontalAlignment','right')
   string = sprintf('%8.2f',sigma(7));
   text(0.80,0.2,string,'HorizontalAlignment','right')
end

%
% CUTS TEXT
%
subplot('Position',[0.10 0.05 0.90 0.10])
%plot([])
axis('off')
text(0.10,1.0,'CUTS','HorizontalAlignment','right')
text(0.10,0.8,'Size','HorizontalAlignment','right')
string = sprintf('%8d',cut.lower_size);
text(0.20,0.8,string,'HorizontalAlignment','right')
string = sprintf('%8d',cut.upper_size);
text(0.30,0.8,string,'HorizontalAlignment','right')
text(0.10,0.6,'Trigger','HorizontalAlignment','right')
string = sprintf('%8d',cut.trigger1);
text(0.20,0.6,string,'HorizontalAlignment','right')
string = sprintf('%8d',cut.trigger2);
text(0.30,0.6,string,'HorizontalAlignment','right')
string = sprintf('%8d',cut.trigger3);
text(0.40,0.6,string,'HorizontalAlignment','right')
text(0.10,0.4,'Frac3','HorizontalAlignment','right')
string = sprintf('%8.2f',cut.frac);
text(0.20,0.4,string,'HorizontalAlignment','right')
text(0.10,0.2,'Width','HorizontalAlignment','right')
string = sprintf('%8.2f',cut.lower_width);
text(0.20,0.2,string,'HorizontalAlignment','right')
string = sprintf('%8.2f',cut.upper_width);
text(0.30,0.2,string,'HorizontalAlignment','right')
text(0.10,0.0,'Length','HorizontalAlignment','right')
string = sprintf('%8.2f',cut.lower_length);
text(0.20,0.0,string,'HorizontalAlignment','right')
string = sprintf('%8.2f',cut.upper_length);
text(0.30,0.0,string,'HorizontalAlignment','right')
text(0.60,0.8,'Distance','HorizontalAlignment','right')
string = sprintf('%8.2f',cut.lower_distance);
text(0.70,0.8,string,'HorizontalAlignment','right')
string = sprintf('%8.2f',cut.upper_distance);
text(0.80,0.8,string,'HorizontalAlignment','right')
text(0.60,0.6,'Alpha','HorizontalAlignment','right')
string = sprintf('%8.2f',cut.alfa);
text(0.70,0.6,string,'HorizontalAlignment','right')
text(0.60,0.4,'Aperture','HorizontalAlignment','right')
string = sprintf('%8.2f',cut.theta);
text(0.70,0.4,string,'HorizontalAlignment','right')
text(0.60,0.2,'Ap. Center','HorizontalAlignment','right')
string = sprintf('%8.2f',cut.center_x);
text(0.70,0.2,string,'HorizontalAlignment','right')
string = sprintf('%8.2f',cut.center_y);
text(0.80,0.2,string,'HorizontalAlignment','right')
text(0.60,0.0,'Tr. Ratio','HorizontalAlignment','right')
string = sprintf('%8.2f',alfaratio);
text(0.70,0.0,string,'HorizontalAlignment','right')
string = sprintf('%8.2f',thetaratio);
text(0.80,0.0,string,'HorizontalAlignment','right')
text(0.10,-0.2,'Length/Size','HorizontalAlignment','right')
string = sprintf('%8.2g',cut.length_size);
text(0.20,-0.2,string,'HorizontalAlignment','right')
text(0.60,-0.2,'Asymm.','HorizontalAlignment','right')
string = sprintf('%8.2g',cut.asymm);
text(0.70,-0.2,string,'HorizontalAlignment','right')
%
% ALPHA
%
[xx,yy] = stairs(alfax,onalfa);
subplot('position',[0.10 0.70 0.35 0.25]),plot(xx,yy,'k-')
if strcmp(mode,'Pair')
   hold on
   [xx,yy] = stairs(alfax,offalfa);
   plot(xx,yy,'k--')
   hold off
   axis([0 90 0 max([onalfa ; offalfa])+max([onalfa ; offalfa])/20])
else
   axis([0 90 0 max(onalfa)+max(onalfa)/20])
end
xlabel('\alpha [deg]')
ylabel('cnts');
%
%THETA
%
[xx,yy] = stairs(thetax,ontheta./areax);
subplot('position',[0.10 0.35 0.35 0.25]),plot(xx,yy,'k-')
hold on
if strcmp(mode,'Pair')
   [xx,yy] = stairs(thetax,offtheta./areax);
   plot(xx,yy,'k--')
   axis([0 2.0 0 max([ontheta./areax ; offtheta./areax])+...
        max([ontheta./areax ; offtheta./areax])/20]);
else
   axis([0 2.0 0 max(ontheta./areax)+max(ontheta./areax)/20]);
end
hold off
xlabel('theta distance from center [deg]')
ylabel('surface brightness [cnts deg^{-2}]')
%
% PlotOpt
%
% RA DEC COORDINATES
%
daysinmonth = [31 28 31 30 31 30 31 31 30 31 30 31];
[iyy,imm,idd,fr,j] = sla_djcl(source_date);
dayofyear = sum(daysinmonth(1:imm-1)) + idd;
source_epoch = iyy + dayofyear/365;
[source_ra source_dec] = precessto2000(source_ra,source_dec,source_epoch);
[source_ra source_dec] = precessfrom2000(source_ra,source_dec,epoch);
if strcmp(plotopt,'Image')
   subplot('position',[0.55 0.35 0.35 0.25]),pcolor(imsigmax,imsigmay,imsigma)
   axis('square')
   if strcmp(mode,'Pair')
      caxis([min(min(imsigma)) max([max(imsigma) maxsig])]);
   end
   colorbar
   colormap(jet)
   if strcmp(mode,'Pair')
      hold on
      contour(imsigmax,imsigmay,imsigma,...
              1:ceil(max(max(imsigma))/10):ceil(max(max(imsigma))),'k')
   end
   s = getstars(source_ra,source_dec,source_epoch,source_ra,6,3);
   plotstars(s);
   if strcmp(grid,'RaDec')
      shading interp
      axis('off')
      radeclines(source_ra,source_dec,-1.9,1.9,-1.9,1.9,0.8,'w');
   else
      xlabel('x[deg]')
      ylabel('y[deg]')
   end
   str = sprintf('RA: %.1f DEC: %.1f J%.1f',radtohhmmss(source_ra),...
                  radtoddmmss(source_dec),epoch);
   title(str);
elseif strcmp(plotopt,'LightCurve')
%
% Plot Light Curve
%
   if strcmp(mode,'Pair')
      subplot('position',[0.55 0.35 0.35 0.25]),plot(lightcurve(:,1),...
              lightcurve(:,5),'o');
      [minrate,iminrate] = min(lightcurve(:,5));
      [maxrate,imaxrate] = max(lightcurve(:,5));
      axis([min(lightcurve(:,1))-1 max(lightcurve(:,1))+1 ...
           minrate-2*lightcurve(iminrate,12) maxrate+2*lightcurve(imaxrate, ...
           12)]);
      set(gca,'YGrid','on');
      for i=1:size(lightcurve,1)
         line([lightcurve(i,1) lightcurve(i,1)],...
              [lightcurve(i,5)+lightcurve(i,12)...
              lightcurve(i,5)-lightcurve(i,12)]);
      end
   elseif strcmp(mode,'Tracking')
      subplot('position',[0.55 0.35 0.35 0.25]),plot(lightcurve(:,1),...
              lightcurve(:,7),'o');
      [minrate,iminrate] = min(lightcurve(:,7));
      [maxrate,imaxrate] = max(lightcurve(:,7));
      axis([min(lightcurve(:,1))-1 max(lightcurve(:,1))+1 ...
           minrate-2*lightcurve(iminrate,14) maxrate+2*lightcurve(imaxrate, ...
           14)]);
      set(gca,'YGrid','on');
      for i=1:size(lightcurve,1)
         line([lightcurve(i,1) lightcurve(i,1)],...
              [lightcurve(i,7)+lightcurve(i,14)...
              lightcurve(i,7)-lightcurve(i,14)]);
      end
   end
   xlabel('Modified Julian Date');
   ylabel('\gamma min^{-1}');
   str = sprintf('RA: %.1f DEC: %.1f J%.1f',radtohhmmss(source_ra),...
                  radtoddmmss(source_dec),epoch);
   title(str);
elseif strcmp(plotopt,'Phase')
   [xx,yy] = stairs(phasex,onphased);
   subplot('position',[0.55 0.35 0.35 0.25]),plot(xx,yy,'k-');
   onchisq = sum(((onphased(1:length(onphased)/2) - ...
                     mean(onphased(1:length(onphased/2)))).^2)./...
                     (onphased(1:length(onphased)/2)));
   if strcmp(mode,'Pair')
      hold on
      [xx,yy] = stairs(phasex,offphased);
      plot(xx,yy,'k--')
      offchisq = sum(((offphased(1:length(offphased)/2) - ...
                     mean(offphased(1:length(offphased/2)))).^2)./...
                     (offphased(1:length(offphased)/2)));
      hold off
      axis([0 2.0 min([onphased ; offphased]) ...
                   max([onphased ; offphased]) + ...
                   max([onphased ; offphased])/20]);
      str = sprintf('\\chi^2=%.2f/%.2f',onchisq,offchisq);
   else
      axis([0 2.0 min(onphased) max(onphased)+max(onphased)/20])
      str = sprintf('\\chi^{2}=%.2f (P=%.2g)',...
            onchisq,1-gammainc(onchisq/2,(length(onphased)/2-1)/2));
   end
   h = title(str);
   xlabel('Phase');
   ylabel('cnts');
elseif strcmp(plotopt,'Spectrum')
   energy = 10.^[2.4:.1:3.6]; % midpoints of energy bins (equal log width)
   non = hist(onenergy,energy);
   if strcmp(mode,'Pair')
      noff = hist(offenergy,energy);
      subplot('position',[0.55 0.35 0.35 0.25]),errorbar(energy,non-noff,...
         sqrt(non+noff),'k-');
   else
      subplot('position',[0.55 0.35 0.35 0.25]),errorbar(energy,non,...
         sqrt(non),'k-');
   end
   xlabel 'size [d.c.]'
   ylabel 'cnts'
end
%
% fill structure of results to return
%
on = struct('lightcurve',lightcurve,...
              'alfa',alfax,'alfacnts',onalfa,...
              'theta',thetax,'thetacnts',ontheta,...
              'phased',phasex,'phasecnts',onphased,...
              'time',ontime,'solar_bary',onsolar_bary,...
              'orbit_bary',onorbit_bary,'energy',onenergy,...
              'phase',onphase,...
              'imagex',imsigmax,'imagey',imsigmay,'image',onimage);
if strcmp(mode,'Pair')
   off = struct('lightcurve',lightcurve,...
                 'alfa',alfax,'alfacnts',offalfa,...
                 'theta',thetax,'thetacnts',offtheta,...
                 'phased',phasex,'phasecnts',offphased,...
                 'time',offtime,'solar_bary',offsolar_bary,...
                 'orbit_bary',offorbit_bary,'energy',offenergy,...
                 'phase',offphase,...
                 'imagex',imsigmax,'imagey',imsigmay,'image',offimage);
end
