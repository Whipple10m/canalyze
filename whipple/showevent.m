function showevent(infile,date,n2id,path)
warning off
%
% SHOWEVENT  Interactive Event Display
%    SHOWEVENT('gtxxxxxxxx.hdf',971107,'gt009386','/usr/local/db/')
%    The first argument is a HDF file name.
%        second argument is the date
%        third argument is the nitrogen runid
%        fourth argument is the path to the database files
%    **New version** reads HDF files directly
%
global runinfo gui_handles size_l_cv
%
% find number of adcs/pmts
%
file_id = hdfh('open',infile,'DFACC_RDONLY',0);
hdfv('start',file_id);
% runinfo
vruninfo_ref = hdfvs('find',file_id,'10M Run Information');
vruninfo_id = hdfvs('attach',file_id,vruninfo_ref,'r');
status = hdfvs('setfields',vruninfo_id,'NADC');
[inadc,count] = hdfvs('read',vruninfo_id,1);
nadc = double(inadc{1});
if nadc == 120
   npmt = 109;
   fov = 4;
elseif nadc == 156
   npmt = 151;
   fov = 5;
elseif nadc == 336
   npmt = 331;
   fov = 6;
elseif nadc == 492
   npmt = 490;
   fov = 4;
else
   disp(sprintf('**error** no valid PMT configuration for nadc = %d',nadc));
end
disp(sprintf('**info** %d pixel camera configuration, fov %d [deg]',npmt,fov));
% events
vevent_ref = hdfvs('find',file_id,'10M Event');
vevent_id = hdfvs('attach',file_id,vevent_ref,'r');
status = hdfvs('setfields',vevent_id,'CODE,ADC');
[count,status] = hdfvs('Querycount',vevent_id);
%
% load coords/peds/gains
%
eval(sprintf('load whipplecoords_%dRWL',npmt));
eval(sprintf('load whippleneighs_%dRWL',npmt));
disp '**info** loaded pixel coordinates and neighbours'
peds = getcpeds(infile(1:8),date,path);
disp '**info** loaded pedestals' 
pedvars = peds(nadc+1:nadc*2);
peds(nadc+1:nadc*2) = [];
if strcmp(n2id,'none')
   disp '**info** no gains applied'
   gains = ones(nadc,1);
else
   gains = getcn2gains(n2id,date,path);
   disp '**info** loaded gains'
end
tubesoff = getcoff(infile(1:8),date,path);
%
% runinfo storage
%
runinfo = struct('vevent_id',vevent_id,'vruninfo_id',vruninfo_id,...
		 'file_id',file_id,'nadc',nadc,'npmt',npmt,'fov',fov,...
                 'whipplecoords',whipplecoords,...
                 'whippleneighs',whippleneighs,...
                 'peds',peds,'pedvars',pedvars,'gains',gains,...
                 'tubesoff',tubesoff);
figure(3);
figure(2);
figure(1);
%
% gui
%
set(gcf,'Position',[10 10 800 400])
h = title('Whipple 10 m Event Display');
set(h,'FontWeight','bold')
axis 'off'
text(.05,.95,'Image Processing');
pedopt = uicontrol(gcf,'Style','radio','Units','normalized',...
                       'Position',[.05 .8 .2 .05],...
                       'String','Apply Pedestals',...
                       'Value',0);
toffopt = uicontrol(gcf,'Style','radio','Units','normalized',...
                        'Position',[.05 .75 .2 .05],...
                        'String','Set Tubes Off',...
                        'Value',0);
cleanopt = uicontrol(gcf,'Style','radio','Units','normalized',...
                         'Position',[.05 .70 .2 .05],...
                         'String','Apply Cleaning',...
                         'Value',0);
gainopt = uicontrol(gcf,'Style','radio','Units','normalized',...
                        'Position',[.05 .65 .2 .05],...
                        'String','Apply Gains',...
                        'Value',0);
paropt = uicontrol(gcf,'Style','radio','Units','normalized',...
                         'Position',[.05 .60 .2 .05],...
                         'String','Parameterize',...
                         'Value',0);
cutopt = uicontrol(gcf,'Style','radio','Units','normalized',...
                         'Position',[.05 .55 .2 .05],...
                         'String','Apply cuts',...
                         'Value',0);
text(.05,.35,'Event Control');
pbquit = uicontrol(gcf,'Style','push','Units','normalized',...
                       'Position',[.05 .1 .2 .05],...
                       'String','Quit','CallBack',...
                       'global runinfo;close all;status=hdfvs(''detach'',runinfo.vruninfo_id);status=hdfvs(''detach'',runinfo.vevent_id);hdfv(''end'',runinfo.file_id);hdfh(''close'',runinfo.file_id);');
pbnextevent = uicontrol(gcf,'Style','push','Units','normalized',...
                            'Position',[.05 .15 .2 .05],...
                            'String','Next Event','CallBack','shownext(1)');
pbprevevent = uicontrol(gcf,'Style','push','Units','normalized',...
                            'Position',[.05 .20 .2 .05],...
                            'String','Prev Event','CallBack','shownext(-1)');
pbprint = uicontrol(gcf,'Style','push','Units','normalized',...
                        'Position',[.05 .25 .2 .05],...
                        'String','Print Event','CallBack',...
                        'print -depsc -f2 shownext.ps');
pbsave = uicontrol(gcf,'Style','push','Units','normalized',...
                        'Position',[.05 .30 .2 .05],...
                        'String','Save Event','CallBack',...
                        'print -djpeg -f2 shownext.jpg');
text(.4,.95,'Event Selection');
%
% lower size cut
%
sli_size_l = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.40 .75 .20 .05],'Min',0,'Max',100000,...
			'Value',0,'CallBack',...
		        'global gui_handles; set(gui_handles.size_l_cv,''String'',num2str(get(gui_handles.sli_size_l,''Val'')))');
size_l_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.35 .75 .05 .05],...
			'String',num2str(get(sli_size_l,'Min')));
size_l_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.60 .75 .05 .05],...
			'String',num2str(get(sli_size_l,'Max')),...
			'Horizontal','right');
size_l_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.40 .8 .1 .05],...
			'String','size >','Horizontal','right');
size_l_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.50 .8 .1 .05],...
			'String',num2str(get(sli_size_l,'Value')),...
			'Horizontal','left');
%
% upper size cut
%
sli_size_u = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.70 .75 .20 .05],'Min',0,'Max',100000,...
			'Value',100000,'CallBack',...
		        'global gui_handles; set(gui_handles.size_u_cv,''String'',num2str(get(gui_handles.sli_size_u,''Val'')))');
size_u_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.65 .75 .05 .05],...
			'String',num2str(get(sli_size_u,'Min')));
size_u_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.90 .75 .05 .05],...
			'String',num2str(get(sli_size_u,'Max')),...
			'Horizontal','right');
size_u_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.70 .8 .1 .05],...
			'String','size <','Horizontal','right');
size_u_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.80 .8 .1 .05],...
			'String',num2str(get(sli_size_u,'Value')),...
			'Horizontal','left');
%
% lower max2 cut
%
sli_max2_l = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.40 .65 .20 .05],'Min',0,'Max',1000,...
			'Value',0,'CallBack',...
		        'global gui_handles; set(gui_handles.max2_l_cv,''String'',num2str(get(gui_handles.sli_max2_l,''Val'')))');
max2_l_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.35 .65 .05 .05],...
			'String',num2str(get(sli_max2_l,'Min')));
max2_l_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.60 .65 .05 .05],...
			'String',num2str(get(sli_max2_l,'Max')),...
			'Horizontal','right');
max2_l_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.40 .70 .1 .05],...
			'String','max2 >','Horizontal','right');
max2_l_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.50 .70 .1 .05],...
			'String',num2str(get(sli_max2_l,'Value')),...
			'Horizontal','left');
%
% upper max2 cut
%
sli_max2_u = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.70 .65 .20 .05],'Min',0,'Max',1000,...
			'Value',1000,'CallBack',...
		        'global gui_handles; set(gui_handles.max2_u_cv,''String'',num2str(get(gui_handles.sli_max2_u,''Val'')))');
max2_u_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.65 .65 .05 .05],...
			'String',num2str(get(sli_max2_u,'Min')));
max2_u_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.90 .65 .05 .05],...
			'String',num2str(get(sli_max2_u,'Max')),...
			'Horizontal','right');
max2_u_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.70 .70 .1 .05],...
			'String','max2 <','Horizontal','right');
max2_u_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.80 .70 .1 .05],...
			'String',num2str(get(sli_max2_u,'Value')),...
			'Horizontal','left');
%
% lower los cut
%
sli_los_l = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.40 .55 .20 .05],'Min',0,'Max',0.005,...
			'Value',0,'CallBack',...
	                'global gui_handles; set(gui_handles.los_l_cv,''String'',num2str(get(gui_handles.sli_los_l,''Val'')))');
los_l_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.35 .55 .05 .05],...
			'String',num2str(get(sli_los_l,'Min')));
los_l_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.60 .55 .05 .05],...
			'String',num2str(get(sli_los_l,'Max')),...
			'Horizontal','right');
los_l_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.40 .60 .1 .05],...
			'String','length/size >','Horizontal','right');
los_l_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.50 .60 .1 .05],...
			'String',num2str(get(sli_los_l,'Value')),...
			'Horizontal','left');
%
% upper los cut
%
sli_los_u = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.70 .55 .20 .05],'Min',0,'Max',0.005,...
			'Value',.005,'CallBack',...
	                'global gui_handles; set(gui_handles.los_u_cv,''String'',num2str(get(gui_handles.sli_los_u,''Val'')))');
los_u_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.65 .55 .05 .05],...
			'String',num2str(get(sli_los_u,'Min')));
los_u_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.90 .55 .05 .05],...
			'String',num2str(get(sli_los_u,'Max')),...
			'Horizontal','right');
los_u_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.70 .60 .1 .05],...
			'String','length/size <','Horizontal','right');
los_u_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.80 .60 .1 .05],...
			'String',num2str(get(sli_los_u,'Value')),...
			'Horizontal','left');
%
% lower width cut
%
sli_width_l = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.40 .45 .20 .05],'Min',0,'Max',1.0,...
			'Value',0.05,'CallBack',...
	                'global gui_handles; set(gui_handles.width_l_cv,''String'',num2str(get(gui_handles.sli_width_l,''Val'')))');
width_l_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.35 .45 .05 .05],...
			'String',num2str(get(sli_width_l,'Min')));
width_l_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.60 .45 .05 .05],...
			'String',num2str(get(sli_width_l,'Max')),...
			'Horizontal','right');
width_l_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.40 .50 .1 .05],...
			'String','width >','Horizontal','right');
width_l_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.50 .50 .1 .05],...
			'String',num2str(get(sli_width_l,'Value')),...
			'Horizontal','left');
%
% upper width cut
%
sli_width_u = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.70 .45 .20 .05],'Min',0,'Max',1.0,...
			'Value',0.12,'CallBack',...
	                'global gui_handles; set(gui_handles.width_u_cv,''String'',num2str(get(gui_handles.sli_width_u,''Val'')))');
width_u_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.65 .45 .05 .05],...
			'String',num2str(get(sli_width_u,'Min')));
width_u_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.90 .45 .05 .05],...
			'String',num2str(get(sli_width_u,'Max')),...
			'Horizontal','right');
width_u_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.70 .50 .1 .05],...
			'String','width <','Horizontal','right');
width_u_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.80 .50 .1 .05],...
			'String',num2str(get(sli_width_u,'Value')),...
			'Horizontal','left');
%
% lower length cut
%
sli_length_l = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.40 .35 .20 .05],'Min',0,'Max',1.0,...
			'Value',0.12,'CallBack',...
	                'global gui_handles; set(gui_handles.length_l_cv,''String'',num2str(get(gui_handles.sli_length_l,''Val'')))');
length_l_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.35 .35 .05 .05],...
			'String',num2str(get(sli_length_l,'Min')));
length_l_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.60 .35 .05 .05],...
			'String',num2str(get(sli_length_l,'Max')),...
			'Horizontal','right');
length_l_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.40 .40 .1 .05],...
			'String','length >','Horizontal','right');
length_l_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.50 .40 .1 .05],...
			'String',num2str(get(sli_length_l,'Value')),...
			'Horizontal','left');
%
% upper length cut
%
sli_length_u = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.70 .35 .20 .05],'Min',0,'Max',1.0,...
			'Value',0.25,'CallBack',...
                  	'global gui_handles; set(gui_handles.length_u_cv,''String'',num2str(get(gui_handles.sli_length_u,''Val'')))');
length_u_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.65 .35 .05 .05],...
			'String',num2str(get(sli_length_u,'Min')));
length_u_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.90 .35 .05 .05],...
			'String',num2str(get(sli_length_u,'Max')),...
			'Horizontal','right');
length_u_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.70 .40 .1 .05],...
			'String','length <','Horizontal','right');
length_u_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.80 .40 .1 .05],...
			'String',num2str(get(sli_length_u,'Value')),...
			'Horizontal','left');
%
% lower distance cut
%
sli_distance_l = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.40 .25 .20 .05],'Min',0,'Max',2.0,...
			'Value',0.35,'CallBack',...
	                'global gui_handles; set(gui_handles.distance_l_cv,''String'',num2str(get(gui_handles.sli_distance_l,''Val'')))');
distance_l_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.35 .25 .05 .05],...
			'String',num2str(get(sli_distance_l,'Min')));
distance_l_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.60 .25 .05 .05],...
			'String',num2str(get(sli_distance_l,'Max')),...
			'Horizontal','right');
distance_l_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.40 .30 .1 .05],...
			'String','distance >','Horizontal','right');
distance_l_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.50 .30 .1 .05],...
			'String',num2str(get(sli_distance_l,'Value')),...
			'Horizontal','left');
%
% upper distance cut
%
sli_distance_u = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.70 .25 .20 .05],'Min',0,'Max',2.0,...
			'Value',1.0,'CallBack',...
 	                'global gui_handles; set(gui_handles.distance_u_cv,''String'',num2str(get(gui_handles.sli_distance_u,''Val'')))');
distance_u_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.65 .25 .05 .05],...
			'String',num2str(get(sli_distance_u,'Min')));
distance_u_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.90 .25 .05 .05],...
			'String',num2str(get(sli_distance_u,'Max')),...
			'Horizontal','right');
distance_u_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.70 .30 .1 .05],...
			'String','distance <','Horizontal','right');
distance_u_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.80 .30 .1 .05],...
			'String',num2str(get(sli_distance_u,'Value')),...
			'Horizontal','left');
%
% lower alpha cut
%
sli_alpha_l = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.40 .15 .20 .05],'Min',0,'Max',90,...
			'Value',0,'CallBack',...
 	                'global gui_handles; set(gui_handles.alpha_l_cv,''String'',num2str(get(gui_handles.sli_alpha_l,''Val'')))');
alpha_l_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.35 .15 .05 .05],...
			'String',num2str(get(sli_alpha_l,'Min')));
alpha_l_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.60 .15 .05 .05],...
			'String',num2str(get(sli_alpha_l,'Max')),...
			'Horizontal','right');
alpha_l_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.40 .20 .1 .05],...
			'String','alpha >','Horizontal','right');
alpha_l_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.50 .20 .1 .05],...
			'String',num2str(get(sli_alpha_l,'Value')),...
			'Horizontal','left');
%
% upper alpha cut
%
sli_alpha_u = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.70 .15 .20 .05],'Min',0,'Max',90,...
			'Value',15,'CallBack',...
	                'global gui_handles; set(gui_handles.alpha_u_cv,''String'',num2str(get(gui_handles.sli_alpha_u,''Val'')))');
alpha_u_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.65 .15 .05 .05],...
			'String',num2str(get(sli_alpha_u,'Min')));
alpha_u_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.90 .15 .05 .05],...
			'String',num2str(get(sli_alpha_u,'Max')),...
			'Horizontal','right');
alpha_u_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.70 .20 .1 .05],...
			'String','alpha <','Horizontal','right');
alpha_u_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.80 .20 .1 .05],...
			'String',num2str(get(sli_alpha_u,'Value')),...
			'Horizontal','left');
%
% lower asymmetry cut
%
sli_asymmetry_l = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.40 .05 .20 .05],'Min',-2.0,'Max',2.0,...
			'Value',-2,'CallBack',...
                        'global gui_handles; set(gui_handles.asymmetry_l_cv,''String'',num2str(get(gui_handles.sli_asymmetry_l,''Val'')))');
asymmetry_l_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.35 .05 .05 .05],...
			'String',num2str(get(sli_asymmetry_l,'Min')));
asymmetry_l_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.60 .05 .05 .05],...
			'String',num2str(get(sli_asymmetry_l,'Max')),...
			'Horizontal','right');
asymmetry_l_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.40 .10 .1 .05],...
			'String','asymmetry >','Horizontal','right');
asymmetry_l_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.50 .10 .1 .05],...
			'String',num2str(get(sli_asymmetry_l,'Value')),...
			'Horizontal','left');
%
% upper asymmetry cut
%
sli_asymmetry_u = uicontrol(gcf,'Style','slider','Units','normalized',...
			'Position',[.70 .05 .20 .05],'Min',-2.0,'Max',2.0,...
			'Value',2,'CallBack',...
                        'global gui_handles; set(gui_handles.asymmetry_u_cv,''String'',num2str(get(gui_handles.sli_asymmetry_u,''Val'')))');
asymmetry_u_min = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.65 .05 .05 .05],...
			'String',num2str(get(sli_asymmetry_u,'Min')));
asymmetry_u_max = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.90 .05 .05 .05],...
			'String',num2str(get(sli_asymmetry_u,'Max')),...
			'Horizontal','right');
asymmetry_u_label = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.70 .10 .1 .05],...
			'String','asymmetry <','Horizontal','right');
asymmetry_u_cv = uicontrol(gcf,'Style','text','Units','normalized',...
			'Position',[.80 .10 .1 .05],...
			'String',num2str(get(sli_asymmetry_u,'Value')),...
			'Horizontal','left');
%
% gui handles storage
%
gui_handles = struct('pedopt',pedopt,'toffopt',toffopt,'gainopt',gainopt,...
                     'cleanopt',cleanopt,'paropt',paropt,'cutopt',cutopt,...
                     'sli_size_l',sli_size_l,'size_l_cv',size_l_cv,...
		     'sli_size_u',sli_size_u,'size_u_cv',size_u_cv,...
                     'sli_max2_l',sli_max2_l,'max2_l_cv',max2_l_cv,...
		     'sli_max2_u',sli_max2_u,'max2_u_cv',max2_u_cv,...
                     'sli_los_l',sli_los_l,'los_l_cv',los_l_cv,...
		     'sli_los_u',sli_los_u,'los_u_cv',los_u_cv,...
                     'sli_width_l',sli_width_l,'width_l_cv',width_l_cv,...
		     'sli_width_u',sli_width_u,'width_u_cv',width_u_cv,...
                     'sli_length_l',sli_length_l,'length_l_cv',length_l_cv,...
                     'sli_length_u',sli_length_u,'length_u_cv',length_u_cv,...
                     'sli_distance_l',sli_distance_l,...
		     'distance_l_cv',distance_l_cv,...
                     'sli_distance_u',sli_distance_u,...
		     'distance_u_cv',distance_u_cv,...
                     'sli_alpha_l',sli_alpha_l,'alpha_l_cv',alpha_l_cv,...
		     'sli_alpha_u',sli_alpha_u,'alpha_u_cv',alpha_u_cv,...
                     'sli_asymmetry_l',sli_asymmetry_l,...
		     'asymmetry_l_cv',asymmetry_l_cv,...
                     'sli_asymmetry_u',sli_asymmetry_u,...
		     'asymmetry_u_cv',asymmetry_u_cv);