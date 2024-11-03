function s = getcoff(runid,date,path)
% getcoff
% gets a vector of tubesoff from database file
% usage: s = getcoff(runid,date,path);
%        arguments:  runid  e.g. 'gt007514'
%                    date   e.g. 971107
%                    path   e.g. '/usr/local/db/'
%        returns:    s vector of tubes off
yy = fix(date/10000);
file = [path 'hrc' sprintf('%02d',yy) '.coff'];
fp = fopen(file,'r');
line = fgets(fp);
while line ~= -1,
   [x,cnt] = sscanf(line,'%s',[1 1]);
   id = sprintf('%c%c%c%c%c%c%c%c',x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8));
   [s,cnt] = sscanf(line(length(x)+3:length(line)),'%d');
   if strcmp(id,runid) == 1
      break;
   else
      s = [];
   end
   line = fgets(fp);
end
