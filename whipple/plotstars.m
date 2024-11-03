function y = plotstars(s)
%
hold on
for i=1:size(s,1)
   if s(i,4) > 1
      h1 = plot(s(i,1),s(i,2),'ro');
%      h2 = plot(s(i,1),s(i,2),'r*');
   elseif s(i,4) < 0
      h1 = plot(s(i,1),s(i,2),'bo');
%      h2 = plot(s(i,1),s(i,2),'b*');
   else
      h1 = plot(s(i,1),s(i,2),'yo');
%      h2 = plot(s(i,1),s(i,2),'y*');
   end
   set(h1,'MarkerSize',-5*s(i,3)+35);
%   set(h2,'MarkerSize',-5*s(i,3)+50);
end
hold off