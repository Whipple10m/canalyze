function shownext(index)
%
global runinfo gui_handles
%
%if index == -1
%   if fseek(runinfo.fid,-2*(2*2+3*8+runinfo.nadc*2*2+4*2),0) == -1
%      fseek(runinfo.fid,-1*(2*2+3*8+runinfo.nadc*2*2+4*2),0);
%   end
%end
%
% main event loop
%
counter = 1;
while 1
   [event,count] = hdfvs('read',runinfo.vevent_id,1);
   code = double(event{1});
   adc = double(event{2});
%
% mimclean
%
   if get(gui_handles.pedopt,'Value') == 1
      adc = adc - runinfo.peds;
   end
   if get(gui_handles.toffopt,'Value') == 1
      adc(runinfo.tubesoff) = 0;
   end
   picture = zeros(runinfo.npmt,1);
   boundary = zeros(runinfo.npmt,1);
   if get(gui_handles.cleanopt,'Value') == 1
% get picture tubes
      picture(find(adc > 4.25*runinfo.pedvars)) = 1;
% get boundary tubes
      for k=1:runinfo.npmt
        if adc(k) > 2.25*runinfo.pedvars(k)
          if picture(k) == 0
             for j=2:7
                if runinfo.whippleneighs(k,j) > 0
                   if picture(runinfo.whippleneighs(k,j)) == 1
                      boundary(k) = 1; 
                   end 
                end
             end
          end 
        end 
      end
% zero tubes not in picture or boundary
      for k=1:runinfo.npmt
         if (picture(k) == 0) & (boundary(k) == 0)
            adc(k) = 0;
         end
      end  
   end 
   if get(gui_handles.gainopt,'Value') == 1
      adc = adc.*runinfo.gains;
   end
   adc(runinfo.npmt+1:runinfo.nadc) = [];
%
% parameterize
%
   if get(gui_handles.paropt,'Value') == 1
      xo = 0.0; yo = 0.0;
%
      sumsig = sum(adc);
%
      xrel = runinfo.whipplecoords(:,2) - xo; 
      yrel = runinfo.whipplecoords(:,3) - xo;
      wxbyev = xrel.*adc(1:runinfo.npmt); wybyev = yrel.*adc;
      sumxsig = sum(wxbyev); sumysig = sum(wybyev);
      sumx2sig = sum(xrel.*wxbyev); sumy2sig = sum(yrel.*wybyev);
      sumxysig = sum(xrel.*wybyev); sumx3sig = sum(xrel.*xrel.*wxbyev);
      sumx2ysig = sum(xrel.*xrel.*wybyev); sumxy2sig = sum(xrel.*yrel.*wybyev);
      sumy3sig = sum(yrel.*yrel.*wybyev);
%
      xmean = sumxsig/sumsig; x2mean = sumx2sig/sumsig;
      ymean = sumysig/sumsig; y2mean = sumy2sig/sumsig;
      xymean = sumxysig/sumsig; x3mean = sumx3sig/sumsig;
      x2ymean = sumx2ysig/sumsig; xy2mean = sumxy2sig/sumsig;
      y3mean = sumy3sig/sumsig; xmean2 = xmean*xmean;
      ymean2 = ymean*ymean; meanxy = xmean*ymean;
%
      sdevx2 = x2mean - xmean2; sdevy2 = y2mean - ymean2;
      sdevxy = xymean - meanxy;
      sdevx3 = x3mean - 3.0*xmean*x2mean + 2.0*xmean2*xmean;
      sdevy3 = y3mean - 3.0*ymean*y2mean + 2.0*ymean2*ymean;
      sdevx2y = x2ymean - x2mean*ymean - 2.0*xymean*xmean + 2.0*xmean2*ymean;
      sdevxy2 = xy2mean - y2mean*xmean - 2.0*xymean*ymean + 2.0*xmean*ymean2;
%
      d = sdevy2 - sdevx2;
      temp = d^2 + 4.0*sdevxy^2;
      z = real(sqrt(temp));
% distance
      temp = xmean2 + ymean2;
      distance = real(sqrt(temp)); 
% length
      temp = (sdevx2 + sdevy2 + z)/2.0;
      length = real(sqrt(temp));
% width
      temp = (sdevy2 + sdevx2 - z)/2.0;
      width = real(sqrt(temp));
% miss 
      if z > 0.0 
         u = 1.0 + d/z; 
         v = 2.0 - u;
         temp = (u*xmean2 + v*ymean2)/2.0 - meanxy*(2.0*sdevxy/z);
         miss = real(sqrt(temp));
      else
         miss = distance;
      end 
% alpha
      alpha = asin(miss/distance)*180.0/pi; 
% centroid
      xc = xmean; yc = ymean; 
% size and max
      size = sumsig;
      maximum = flipud(sort(adc));
% major axis
     temp = d^2 + 4.0*sdevxy^2;
     if temp > 0.0
        if  sdevxy == 0.0
           slope = 999.;
        else
           slope = (d + sqrt(temp))/2.0/sdevxy;
           yint = ymean - slope*xmean;
        end 
     end
% asymmetry
     if length == 0.0
        asymm = 0.0;
     else
        psi = atan2((d+z)*ymean+2.0*sdevxy*xmean,2.0*sdevxy*ymean-(d-z)*xmean);
        cpsi=cos(psi); spsi=sin(psi);
        asymmetry = sdevx3*cpsi^3+3.0*sdevx2y*spsi*cpsi^2+...
                3.0*sdevxy2*cpsi*spsi^2+sdevy3*spsi^3;
        if asymmetry < 0
           asymmetry = -exp(log(abs(asymmetry))/3.0)/length;
        else
           asymmetry = exp(log(abs(asymmetry))/3.0)/length;
        end
      end
% point of origin
      disp = 1.37*(1 - width/length);
      qa = slope^2 + 1;
      qb = 2*slope*yint - 2*ymean*slope - 2*xmean;
      qc = xmean^2 + ymean^2 - 2*ymean*yint + yint^2 - disp^2;
      xoa = (-qb - sqrt(qb*qb - 4.0*qa*qc))/2.0/qa;
      yoa = slope*xoa + yint;
      xob = (-qb + sqrt(qb*qb - 4.0*qa*qc))/2.0/qa;
      yob = slope*xob + yint;
   end
%
% display event
%
   if get(gui_handles.cutopt,'Value') ~= 1 | ...
      (size > get(gui_handles.sli_size_l,'Val') & ...
       size < get(gui_handles.sli_size_u,'Val') & ...
       maximum(2) > get(gui_handles.sli_max2_l,'Val') & ...
       maximum(2) < get(gui_handles.sli_max2_u,'Val') & ...
       length/size > get(gui_handles.sli_los_l,'Val') & ...
       length/size < get(gui_handles.sli_los_u,'Val') & ...
       width > get(gui_handles.sli_width_l,'Val') & ...
       width < get(gui_handles.sli_width_u,'Val') & ...
       length > get(gui_handles.sli_length_l,'Val') & ...
       length < get(gui_handles.sli_length_u,'Val') & ...
       distance > get(gui_handles.sli_distance_l,'Val') & ...
       distance < get(gui_handles.sli_distance_u,'Val') & ...
       alpha > get(gui_handles.sli_alpha_l,'Val') & ...
       alpha < get(gui_handles.sli_alpha_u,'Val') & ...
       asymmetry > get(gui_handles.sli_asymmetry_l,'Val') & ...
       asymmetry < get(gui_handles.sli_asymmetry_u,'Val') )
      figure(2);
      clf;
      axis([-runinfo.fov/2 runinfo.fov/2 -runinfo.fov/2 runinfo.fov/2]);
      axis('square');
      hold on
      color = [0 0 1];
      if runinfo.npmt > 331
         pixel_size = 0.118;
      else
         pixel_size = 0.259;
      end
      maxadc = fix(max(adc));
      for i=1:runinfo.npmt
         if picture(i) == 1
            color = [1 0 0];
         end
         if boundary(i) == 1
            color = [0 1 0];
         end
         if i > 379
            circle(runinfo.whipplecoords(i,2),runinfo.whipplecoords(i,3),...
               .233/2,[1 1 1]);
            if(adc(i) > 1)
               circle(runinfo.whipplecoords(i,2),runinfo.whipplecoords(i,3),...
                  adc(i)/maxadc*.233/2,color);
            end
         else
            circle(runinfo.whipplecoords(i,2),runinfo.whipplecoords(i,3),...
               pixel_size/2,[1 1 1]);
            if(adc(i) > 1)
               circle(runinfo.whipplecoords(i,2),runinfo.whipplecoords(i,3),...
                  adc(i)/maxadc*pixel_size/2,color);
            end
         end
      end
      if get(gui_handles.paropt,'Value') == 1
         plot([xoa xob],[yoa yob],'m*');
         fplot(sprintf('%f*x+%f',slope,yint),[-runinfo.fov/2 runinfo.fov/2]);
         ellipse(length,width,slope,xc,yc);
         figure(3)
         clf
         axis 'off'
         title('Hillas Parameters');
         text(0,1,sprintf('size=%d',size));
         text(0,.95,sprintf('max2=%d',maximum(2)));
         text(0,.90,sprintf('length/size=%.5f',length/size));
         text(0,.85,sprintf('length=%.4f',length));
         text(0,.80,sprintf('width=%.4f',width));
         text(0,.75,sprintf('distance=%.4f',distance));
         text(0,.70,sprintf('alpha=%.4f',alpha));
         text(0,.65,sprintf('asymmetry=%.4f',asymmetry));
      end
      break;
   else
      if get(gui_handles.cutopt,'Value') == 1
         figure(3)
         clf
         axis 'off'
         title('Hillas Parameters');
         text(0,1,sprintf('Failed: %d times',counter));
         counter = counter + 1;   
      end
   end
end
