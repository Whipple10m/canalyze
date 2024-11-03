function s = getcn2gains(runid,date,path)
% getcn2gains
% gets a vector of gain factors from database file
% usage: s = getcn2gains(runid,date,path);
%        arguments:  runid  e.g. 'gt007514'
%                    date   e.g. 971107
%                    path   e.g. '/usr/local/db/'
%        returns:    s vector of gain factors
yy = fix(date/10000);
file = [path 'hrc' sprintf('%02d',yy) '.cn2gains'];
fp = fopen(file,'r');
[x,cnt] = fscanf(fp,'%s %d %d',[3 1]);
while cnt == 3,
   id = sprintf('%c%c%c%c%c%c%c%c',x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8));
   dt = x(9);
   nadc = x(10);
   [s,cnt] = fscanf(fp,'%f',[nadc 1]);
   if strcmp(id,runid) == 1 & dt == date
      break;
   else
      s = zeros(1,nadc);
   end
   [x,cnt] = fscanf(fp,'%s %d %d',[3 1]);
end
