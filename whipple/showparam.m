function showparam(varargin)
% SHOWPARAM  Parameter Distribution Display
%    SHOWPARAM('gtxxxxxxxxp.hdf','ON'...);
%    The first argument is a HDF parameter file name.
%        second argument is the LEGEND text string
%        repeat as needed
%    **NEW** reads HDF file directly
%
if mod(nargin,2)
   error('invalid number of arguments')
   return
end
set(gcf,'Position',[360 300 523 638])
set(gcf,'PaperPosition',[0.5 1.0 7.5 9.0])
clmap = [ 0 0 1 ; 1 0 0 ; 0 1 0; 0.6667 0 1];
for i=1:2:nargin
   infile = varargin{i};
   disp(['displaying ' infile]);
   legendstring{ceil(i/2)} = varargin{i+1};
   file_id = hdfh('open',infile,'DFACC_RDONLY',0);
   hdfv('start',file_id);
   vpar_ref = hdfvs('find',file_id,'10M Parameters');
   vpar_id = hdfvs('attach',file_id,vpar_ref,'r');
   status = hdfvs('setfields',vpar_id,...
      'GPSTIME,OSCTIME,LIVETIME,LENGTH,WIDTH,MISS,DISTANCE,AZWIDTH,FRAC,SIZE,LOC,MAX,ASYMM,XC,YC,XO,YO');
   [count,status] = hdfvs('Querycount',vpar_id);
   [par,count] = hdfvs('read',vpar_id,count);
   status = hdfvs('detach',vpar_id);
   hdfv('end',file_id);
   hdfh('close',file_id);
%
% frac2
%
   subplot(3,2,1);
   hold on;
   [xx,yy] = hist(double(par{9}(2,:)),0:.02:1);
   h = stairs(yy,xx);
   set(h,'Color',clmap(ceil(i/2),:));
   if i == 1
      maxx = max(xx);
   else
      maxx = max([max(xx) get(gca,'YLim')]);
   end
   axis([0 1 0 maxx+maxx/10]);
   legend(legendstring{:});
   xlabel 'frac2'
%
% frac 3
%
   subplot(3,2,2);
   hold on
   [xx,yy] = hist(double(par{9}(3,:)),0:.02:1);
   h = stairs(yy,xx);
   set(h,'Color',clmap(ceil(i/2),:));
   if i == 1
      maxx = max(xx);
   else
      maxx = max([max(xx) get(gca,'YLim')]);
   end
   axis([0 1 0 maxx+maxx/10]);
   legend(legendstring{:});
   xlabel 'frac3'
%
% length
%
   subplot(3,2,3);
   hold on
   [xx,yy] = hist(double(par{4}),0:.02:1);
   h = stairs(yy,xx);
   set(h,'Color',clmap(ceil(i/2),:));
   if i == 1
      maxx = max(xx);
   else
      maxx = max([max(xx) get(gca,'YLim')]);
   end
   axis([0 1 0 maxx+maxx/10]);
   legend(legendstring{:});
   xlabel 'length'
%
% width
%
   subplot(3,2,4);
   hold on
   [xx,yy] = hist(double(par{5}),0:.01:.6);
   h = stairs(yy,xx);
   set(h,'Color',clmap(ceil(i/2),:));
   if i == 1
      maxx = max(xx);
   else
      maxx = max([max(xx) get(gca,'YLim')]);
   end
   axis([0 0.6 0 maxx+maxx/10]);
   legend(legendstring{:});
   xlabel 'width'
%
% distance
%
   subplot(3,2,5);
   hold on
   [xx,yy] = hist(double(par{7}),0:.05:1.5);
   h = stairs(yy,xx);
   set(h,'Color',clmap(ceil(i/2),:));
   if i == 1
      maxx = max(xx);
   else
      maxx = max([max(xx) get(gca,'YLim')]);
   end
   axis([0 1.5 0 maxx+maxx/10]);
   legend(legendstring{:},2);
   xlabel 'distance'
%
% alpha
%
   subplot(3,2,6);
   hold on
   [xx,yy] = hist(asin(double(par{6})./double(par{7}))*180/pi,2.5:5:87.5);
   h = stairs(yy,xx);
   set(h,'Color',clmap(ceil(i/2),:));
   if i == 1
      maxx = max(xx);
   else
      maxx = max([max(xx) get(gca,'YLim')]);
   end
   axis([0 90 0 maxx+maxx/10]);
   legend(legendstring{:},2);
   xlabel 'alpha'
end
