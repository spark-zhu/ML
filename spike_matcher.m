
function varargout = spike_matcher(varargin)
% SPIKE_MATCHER MATLAB code for spike_matcher.fig
%      SPIKE_MATCHER, by itself, creates a new SPIKE_MATCHER or raises the existing
%      singleton*.
%
%      H = SPIKE_MATCHER returns the handle to a new SPIKE_MATCHER or the handle to
%      the existing singleton*.
%
%      SPIKE_MATCHER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPIKE_MATCHER.M with the given input arguments.
%
%      SPIKE_MATCHER('Property','Value',...) creates a new SPIKE_MATCHER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spike_matcher_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spike_matcher_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spike_matcher

% Last Modified by GUIDE v2.5 21-Sep-2017 19:20:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spike_matcher_OpeningFcn, ...
                   'gui_OutputFcn',  @spike_matcher_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before spike_matcher is made visible.
function spike_matcher_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spike_matcher (see VARARGIN)

% Choose default command line output for spike_matcher
handles.output = hObject;
load('store.mat') 

% store.AllData = AllData;
% store.MaskFeatureAll = MaskFeatureAll;
% store.DataForFeature=DataForFeature;
% store.lookUpTable=lookUpTable;
% store.DPS=DPS;
% store.files=files;
% store.LocForFeature=LocForFeature;
% store.TimeFeatureAll=TimeFeatureAll;
global store
store=store1;
set(handles.listbox3,'String',{store.files.name})
% global DataForFeature
% global MaskFeatureAll
% global colors
% global lookUpTable
% global unitNums
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spike_matcher wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = spike_matcher_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Add.
function Add_Callback(hObject, eventdata, handles)
% hObject    handle to Add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Groups AloneUnitList unitsToBeLooped I_Dis his unit path lookUpTable


if path==1
    if lookUpTable(unitsToBeLooped(unit),1)~=lookUpTable(AloneUnitList(I_Dis(his)),1)
        Groups{end+1}=[AloneUnitList(I_Dis(his)),unitsToBeLooped(unit)];

AloneUnitList(I_Dis(his))=[];
set(handles.listbox1,'String',displayGroups_special(Groups))
unitsToBeLooped(unit)=[];
his=0;


end

elseif path==2
    list = Groups{I_Dis(his)};
    if ~ismember(lookUpTable(unitsToBeLooped(unit),1),lookUpTable(list,1))
         list(end+1) = unitsToBeLooped(unit);
    Groups{I_Dis(his)}=list;
    set(handles.listbox1,'String',displayGroups_special(Groups))
unitsToBeLooped(unit)=[];
his=0;

end
   
    

end
% AloneUnitList(I_Dis(his))=[]; % delete from the AloneUnitList This needs to be done. 
% unitsToBeLooped(unit)=[];

% unit=unit+1;
if unit>numel(unitsToBeLooped)
    unit=unit-1;
end

clear_main(handles)

% need to determine ending condition for this session. 
uiresume



% --- Executes on button press in Delete.
function Delete_Callback(hObject, eventdata, handles)
% hObject    handle to Delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Groups AloneUnitList
listbox2value=get(handles.listbox2,'Value');
listbox1value=get(handles.listbox1,'Value');
SelectedGroup = Groups{listbox1value};
unitToBeDeleted=SelectedGroup(listbox2value);

set(handles.edit15,'String',unitToBeDeleted);
SelectedGroup(SelectedGroup==unitToBeDeleted)=[]; % not considering only one unit left alone
AloneUnitList(end+1) = unitToBeDeleted;
Groups{listbox1value}=SelectedGroup;

if numel(SelectedGroup)==1 % perform another deletion
    if listbox2value==2
        unitToBeDeleted=SelectedGroup(listbox2value-1);
    else
        unitToBeDeleted=SelectedGroup(listbox2value);

    end

set(handles.edit15,'String',unitToBeDeleted);
SelectedGroup(SelectedGroup==unitToBeDeleted)=[];
AloneUnitList(end+1) = unitToBeDeleted;
Groups(listbox1value)=[];
end
set(handles.listbox1,'Value',listbox1value-1);
set(handles.listbox1,'String',displayGroups_special(Groups))
set(handles.listbox2,'Value',1);
set(handles.listbox2,'String',displayGroups_special(num2cell(Groups{listbox1value-1})))
clear_main(handles)

% addback to Alone Unit List 

% perfomr deletion

 %AloneUnitList unitsToBeLooped I_Dis his unit



% --- Executes on button press in Senior.
function Senior_Callback(hObject, eventdata, handles)
% hObject    handle to Senior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Junior.
function Junior_Callback(hObject, eventdata, handles)
% hObject    handle to Junior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Overlay.
function Overlay_Callback(hObject, eventdata, handles)
% hObject    handle to Overlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Form.
function Form_Callback(hObject, eventdata, handles)
% hObject    handle to Form (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Dismiss.
function Dismiss_Callback(hObject, eventdata, handles)
% hObject    handle to Dismiss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiopen;
global  Groups
global AloneUnitList
global unitsToBeLooped % unitsToBeLoopedIsAvailable

Groups=result.Groups;
AloneUnitList=result.AloneUnitList;
unitsToBeLooped=result.unitsToBeLooped;
set(handles.listbox1,'String',displayGroups_special(Groups))
set(handles.listbox4,'String',displayGroups(num2cell(AloneUnitList)))
set(handles.listbox5,'String',displayGroups(num2cell(unitsToBeLooped)))



% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global store Groups
global AloneUnitList
global unitsToBeLooped % unitsToBeLoopedIsAvailable
result.store=store;
result.Groups=Groups;
result.AloneUnitList=AloneUnitList;
result.unitsToBeLooped=unitsToBeLooped;

uisave('result','result-')


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
function axes1_CreateFcn(hObject, eventdata, handles)
function axes1_DeleteFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in pushbutton15.
global his
his = his - 1;
if his<1
    his = his+1;
end
clear_main(handles)
uiresume

function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global his
his = his +1;
global AloneUnitList
if his>numel(AloneUnitList) % need to have a full record = grouops + alone unit List
his = his-1;
end
clear_main(handles)
uiresume

function plot_unit(handles,unitNums,colors,oldFlag)
global GUI_axes_map sorted_Dis his MSEList
GUI_axes_map =[1	2	3	4	5	6	7	8
9	10	11	12	13	14	15	16
17	18	19	20	21	22	23	24
25	26	27	28	29	30	31	32
];
t=[-30/20e3:1/20e3:10/20e3]*1000;
global AllData DataForFeature MaskFeatureAll lookUpTable store TimeFeatureAll LocForFeature
AllData=store.AllData;  DataForFeature=store.DataForFeature;
%MaskFeatureAll=store.MaskFeatureAll; 
lookUpTable=store.lookUpTable; TimeFeatureAll=store.TimeFeatureAll;
LocForFeature=store.LocForFeature;
for row=1:4
            for col=1:8
        
                for unit=1:numel(unitNums)
                    Ch_Map=AllData{lookUpTable(unitNums(unit),1)}.Ch_Map;
                    temp=Ch_Map(row,col)-15;
                    if Ch_Map(row,col)>0 %&& MaskFeatureAll{unitNums(unit)}(temp)>0.01
                        eval(['axes(handles.axes' num2str(GUI_axes_map(row,col)) ');' ]    );



                        plot(t,DataForFeature{unitNums(unit)}(temp,:),'Color',colors{unit},'LineWidth',2)
                        set(gca,'xticklabel',[])
                        hold on;
                        axis([t(1) t(end) -175 25])
                    end
                    
                    
                    
                end
            end
end
 for unit=1:numel(unitNums)
 FiringTimeForThisUnit = TimeFeatureAll{unitNums(unit)};
ISI_in_MS_bins = diff(double(FiringTimeForThisUnit))/20; % Intervals in miliseconds. 
[counts] = histc(ISI_in_MS_bins,0:1:50);
[counts_M] = histc(ISI_in_MS_bins,0:20:1000);
[counts_L] = histc(ISI_in_MS_bins,0:1000:20000);
Sec_bins = double(FiringTimeForThisUnit)/20E3;
[counts_5s_1] = histc(Sec_bins,0:5:720);
[counts_5s_2] = histc(Sec_bins,2.5:5:722.5);
Max_Ins_5s_FR = max(max(counts_5s_1),max(counts_5s_2))/5;
FR_avg = numel(FiringTimeForThisUnit)/720;


axes(handles.axes33)
if ~isempty(counts)
stairs(0.5:1:50.5,counts/max(counts),'Color',colors{unit})
xlim([0 51])
end
hold on;

axes(handles.axes34)
if ~isempty(counts_M)
stairs(10:20:1010,counts_M/max(counts_M),'Color',colors{unit})
xlim([0 1020])
end
hold on;        

axes(handles.axes35)
if ~isempty(counts_L)
stairs(0.5:1:20.5,counts_L/max(counts_L),'Color',colors{unit})
axis([0 21 0 1])
end
hold on;
if oldFlag==2
set(handles.edit6,'string',num2str(Max_Ins_5s_FR),'ForegroundColor','r')
set(handles.edit7,'string',num2str(FR_avg),'ForegroundColor','r')

else if unit ==1
    set(handles.edit8,'string',num2str(Max_Ins_5s_FR),'ForegroundColor','b')
set(handles.edit9,'string',num2str(FR_avg),'ForegroundColor','b')

    end
end

 end
%%
 axes(handles.axes37)
%  Ox=12.5;
%         Oy=12.5;
%         gap=5;
%         Gw=25;
%         Gh=25;
        
Ox=15;    
Oy=15;
gap=20;
        Gw=30;
        Gh=30;        

        numGx = 8;
        numGy = 4;
        Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
        Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
        [X,Y] = meshgrid(Gx,Gy);
    
    
    for r=1:4
        for c=1:8
%         scatter(X(r,c),Y(r,c),22000,[0 0 0],'s','LineWidth',2)
        pAx = X(r,c)-12.5;
        pAy = Y(r,c)-12.5;
        pCx = X(r,c)+12.5;
        pCy = Y(r,c)+12.5;
        plot([pAx pCx],[pAy pAy],'k','LineWidth',1)
        hold on
        plot([pCx pCx],[pAy pCy],'k','LineWidth',1)
        hold on 
        plot([pAx pCx],[pCy pCy],'k','LineWidth',1)
        hold on 
        plot([pAx pAx],[pAy pCy],'k','LineWidth',1)
        end
    end
    hold on 
    if oldFlag==2 % new unit
     scatter(LocForFeature{unitNums}(1),LocForFeature{unitNums}(2),30+(LocForFeature{unitNums}(3))^2,'+','LineWidth',2,'MarkerEdgeColor','r')

    else if oldFlag==1 % old unit
     scatter(LocForFeature{unitNums}(1),LocForFeature{unitNums}(2),30+(LocForFeature{unitNums}(3))^2,'+','LineWidth',2,'MarkerEdgeColor','b')
        end
    end
axis([-40 440 -40 200])

  if ~isempty(sorted_Dis) && ~isempty(his)  
 axes(handles.axes38)
 try
  plot(1:10,sorted_Dis(1:10))
 catch
     plot(1:numel(sorted_Dis),sorted_Dis)
 end
  hold on;
  scatter(his,sorted_Dis(his),20,'filled','MarkerFaceColor','r')
  axis([0 11 0 50])
  end
 
   axes(handles.axes36)
  plot(1:10,MSEList)
  hold on;
  scatter(his,MSEList(his),20,'filled','MarkerFaceColor','r')
  axis([0 11 0 50])
 



function clear_main(handles)
global GUI_axes_map
for row=1:4
            for col=1:8
        
              
            eval(['axes(handles.axes' num2str(GUI_axes_map(row,col)) ');' ]    );
            cla
                    
    
            end
end
axes(handles.axes33)
cla
axes(handles.axes34)
cla
axes(handles.axes35)
cla
axes(handles.axes37)
cla
axes(handles.axes36)
cla
axes(handles.axes38)
cla

% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global store Groups
global AloneUnitList
AloneUnitList=find(store.lookUpTable(:,1)==1);
Groups = cell(1,1);
Groups{1} = [0,0];

set(handles.listbox1,'String',displayGroups_special(Groups))
set(handles.listbox4,'String',displayGroups(num2cell(AloneUnitList)))
% test = cell(1,3);
% test{1}='1, 2';
% test{2}='2, 3';
% test{3}='3, 4';
% 
function strGroups = displayGroups(Groups)
strGroups =cell(size(Groups));
for i=1:numel(Groups)
unitList = Groups{i};
strList = ['[' num2str(unitList(1)) ','];
for unit=2:numel(unitList)
    strList=[strList num2str(unitList(unit)) ','];
end
    strList(end)=']';
    strGroups{i}=strList;
end

function strGroups = displayGroups_special(Groups)
global lookUpTable store TimeFeatureAll
TimeFeatureAll=store.TimeFeatureAll;
 lookUpTable=store.lookUpTable;
strGroups =cell(size(Groups));
label='0';
for i=1:numel(Groups)
    
    if sum(Groups{i})>0
unitList = lookUpTable(Groups{i},1);
GrptimeFeatures = giveGrpTime(Groups{i});
counts = GrptimeFeatures.counts;
label='S';
if counts(1)>23
    label='M';
end
    else
unitList =Groups{i};
        
    end
    
strList = [label '['  num2str(unitList(1)) ','];
for unit=2:numel(unitList)
    strList=[strList num2str(unitList(unit)) ','];
end
    strList(end)=']';
    strGroups{i}=strList;

end


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global store
global AloneUnitList
global his unit unitsToBeLooped I_Dis breakFlag% unitsToBeLoopedIsAvailable
global sessionNum 
sessionNum = str2num(get(handles.edit16,'String'));
unitsToBeLooped = find(store.lookUpTable(:,1)==sessionNum);% loading the second session. 
 set(handles.listbox5,'String',displayGroups(num2cell(unitsToBeLooped)))
% can display them in another list box, 
% unitsToBeLoopedIsAvailable=true(size(unitsToBeLooped));
% unit=1;
% while unit<=numel(unitsToBeLooped) && unit>=1
%     breakFlag=0;
%     set(handles.listbox5,'String',displayGroups(num2cell(unitsToBeLooped)))
%     set(handles.listbox5,'Value',unit);
%     set(handles.edit13,'String',num2str(unitsToBeLooped(unit)));
%     historyUnitLocations = cell2mat(store.LocForFeature(AloneUnitList));
%     thisUnitLocation = store.LocForFeature{unitsToBeLooped(unit)};
%     Distance = pdist2(thisUnitLocation,historyUnitLocations);
%     Two_D_Distance = pdist2(thisUnitLocation(1:2),historyUnitLocations(:,1:2));
%     [sorted_Dis,I_Dis] = sort(Distance);
%     his = 1;
%     while his<=numel(AloneUnitList) && his>=1
%             set(handles.listbox4,'String',displayGroups(num2cell(AloneUnitList(I_Dis))))
%             set(handles.listbox4,'Value',his);
%             set(handles.edit14,'String',num2str(AloneUnitList(I_Dis(his))));
%     
%         % 1. plot the history unit
%          plot_unit(handles,AloneUnitList(I_Dis(his)),cellstr('b'),1)  
%          hold on;
%          plot_unit(handles,unitsToBeLooped(unit),cellstr('r'),2) 
%          set(handles.edit11,'string',num2str(sorted_Dis(his)))
%                   set(handles.edit12,'string',num2str(Two_D_Distance(I_Dis(his))))
% 
%     % do overlay plot , currentUnit and closest history unit. 
%     % wait for user input to do one of the following and break inner loop,
%     % add to existing   
%     uiwait
%     if breakFlag==1
%     break;
%     end
%     
%     end
%     
% end


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global unit his unitsToBeLooped
his=0;
unit=unit+1;
if unit>numel(unitsToBeLooped)
    unit=unit-1;
end
clear_main(handles)
uiresume

% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global unit his
his=0;
unit=unit-1;
if unit<1
    unit=unit+1;
end
clear_main(handles)
uiresume



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
value=get(handles.listbox1,'Value');
set(handles.edit15,'String',value);
global Groups %AloneUnitList unitsToBeLooped I_Dis his unit
set(handles.listbox2,'Value',1);
set(handles.listbox2,'String',displayGroups_special(num2cell(Groups{value})))


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2
clear_main(handles)
uiresume

% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  AloneUnitList his I_Dis

% new = str2num(get(handles.edit13,'String'));
old = str2num(get(handles.edit14,'String'));
% unit = find(unitsToBeLooped==new);
his = find(AloneUnitList(I_Dis)==old);
% 
% if unit>numel(unitsToBeLooped)
%     unit=unit-1;
% end

clear_main(handles)

% need to determine ending condition for this session. 
uiresume



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  AloneUnitList unitsToBeLooped 
unitsToBeLooped=AloneUnitList;
set(handles.listbox5,'String',displayGroups(num2cell(unitsToBeLooped)))


% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4


% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox5.
function listbox5_Callback(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox5


% --- Executes during object creation, after setting all properties.
function listbox5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  unitsToBeLooped unit his

new = str2num(get(handles.edit13,'String'));
% old = str2num(get(handles.edit14,'String'));
unit = find(unitsToBeLooped==new);
% his = find(AloneUnitList(I_Dis)==old);
% 
if unit>numel(unitsToBeLooped)
    unit=unit-1;
end
his=0;
clear_main(handles)

% need to determine ending condition for this session. 
uiresume



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global his
his = his - 1;
if his<1
    his = his+1;
end
clear_main(handles)
uiresume

% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global his
his = his +1;
global Groups
if his>numel(Groups)-1 % need to have a full record = grouops + alone unit List -1 is to tackle filler(0,0);
his = his-1;
end
clear_main(handles)
uiresume

% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Groups AloneUnitList unitsToBeLooped his
listbox5value=get(handles.listbox5,'Value');
SelectedUnit= unitsToBeLooped(listbox5value);

set(handles.edit15,'String',SelectedUnit);
% not considering only one unit left alone
AloneUnitList(end+1) = SelectedUnit;
unitsToBeLooped(listbox5value)=[];

if listbox5value>numel(unitsToBeLooped)
listbox5value=listbox5value-1;
end
set(handles.listbox5,'Value',listbox5value);
set(handles.listbox5,'String',displayGroups(num2cell(unitsToBeLooped)))

if numel(unitsToBeLooped)<1
set(handles.listbox5, 'String', '');
end


clear_main(handles);
his=0;
% need to determine ending condition for this session. 
uiresume

function MSEList  = SortMSE(unitList,unitNow)
global DataForFeature 
for unit=1:numel(unitList)
[~,MSE]=Sim_index(DataForFeature{unitList(unit)},DataForFeature{unitNow});
MSEList(unit)=mean(MSE);
end

function MSEList  = SortMSE_GP(Gps,unitNow)
global DataForFeature 
for unit=1:numel(Gps)
[~,MSE]=Sim_index(giveGrpWaveform(Gps{unit}),DataForFeature{unitNow});
MSEList(unit)=mean(MSE);
end


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global store Groups path
global AloneUnitList DataForFeature
global his unit unitsToBeLooped sorted_Dis breakFlag% unitsToBeLoopedIsAvailable
global sessionNum  MSEList I_Dis
DataForFeature=store.DataForFeature;
sessionNum = str2num(get(handles.edit16,'String'));
% unitsToBeLooped = find(store.lookUpTable(:,1)==sessionNum);% loading the second session. 
% can display them in another list box, 
% unitsToBeLoopedIsAvailable=true(size(unitsToBeLooped));
unit=1;
path=2; % compare with existing groups first
if numel(Groups)==1
    path=1;
end
while unit<=numel(unitsToBeLooped) && unit>=1
    breakFlag=0;
    set(handles.listbox5,'String',displayGroups(num2cell(unitsToBeLooped)))
    set(handles.listbox5,'Value',unit);
    set(handles.edit13,'String',num2str(unitsToBeLooped(unit)));
    
    if path ==1 
    historyUnitLocations = cell2mat(store.LocForFeature(AloneUnitList));
    thisUnitLocation = store.LocForFeature{unitsToBeLooped(unit)};
    Distance = pdist2(thisUnitLocation,historyUnitLocations);
    Two_D_Distance = pdist2(thisUnitLocation(1:2),historyUnitLocations(:,1:2));
    [sorted_Dis,I_Dis] = sort(Distance);

    
      if numel(I_Dis)>=10
    MSEList  =  SortMSE(AloneUnitList(I_Dis(1:10)),unitsToBeLooped(unit));
    else
        I_Dis_mod = I_Dis;
        I_Dis_mod(numel(I_Dis)+1 :10)=I_Dis_mod(end);
            MSEList  = SortMSE(AloneUnitList(I_Dis_mod),unitsToBeLooped(unit));
    end
     % function sortedMSE  = (AloneUnitList(I_Dis),unitsToBeLooped(unit));
    his = 1;
    while his<=numel(AloneUnitList) && his>=1
            set(handles.listbox4,'String',displayGroups(num2cell(AloneUnitList(I_Dis))))
            set(handles.listbox4,'Value',his);
            set(handles.edit14,'String',num2str(AloneUnitList(I_Dis(his))));
    
        % 1. plot the history unit
         plot_unit(handles,AloneUnitList(I_Dis(his)),cellstr('b'),1)  
         hold on;
         plot_unit(handles,unitsToBeLooped(unit),cellstr('r'),2) 
         set(handles.edit11,'string',num2str(sorted_Dis(his)))
                  set(handles.edit12,'string',num2str(Two_D_Distance(I_Dis(his))))
                  
                  [XCOR,MSE]=Sim_index(DataForFeature{AloneUnitList(I_Dis(his))},DataForFeature{unitsToBeLooped(unit)});
                  set(handles.edit15,'string',num2str(mean(MSE)))
                  XCOR(isnan(XCOR))=[];
                  set(handles.edit17,'string',num2str((1-mean(XCOR))*100))

    % do overlay plot , currentUnit and closest history unit. 
    % wait for user input to do one of the following and break inner loop,
    % add to existing   
 
    uiwait
    if breakFlag==1
    break;
    end
    end
    
    end
    
    if path ==2
    % path 2 
    historyGrpLocations = GroupsLocation();% some function given Groups variable, return locations.;
    thisUnitLocation = store.LocForFeature{unitsToBeLooped(unit)};
    Distance = pdist2(thisUnitLocation,historyGrpLocations);
    Two_D_Distance = pdist2(thisUnitLocation(1:2),historyGrpLocations(:,1:2));
    [sorted_Dis,I_Dis] = sort(Distance);
    if numel(I_Dis)>=10
    MSEList  = SortMSE_GP(Groups(I_Dis(1:10)),unitsToBeLooped(unit));
    else
        I_Dis_mod = I_Dis;
        I_Dis_mod(numel(I_Dis)+1 :10)=I_Dis_mod(end);
            MSEList  = SortMSE_GP(Groups(I_Dis_mod),unitsToBeLooped(unit));
    end
    his = 1;
    while his<=numel(Groups) && his>=1
            set(handles.listbox1,'String',displayGroups_special(Groups(I_Dis))) 
            set(handles.listbox1,'Value',his);
            set(handles.edit14,'String',num2str(Groups{I_Dis(his)}));
    
        % 1. plot the history unit
         plot_Group(handles,Groups{I_Dis(his)},cellstr('b'),1)  % some function plotting all units in a group.
         hold on;
         plot_unit(handles,unitsToBeLooped(unit),cellstr('r'),2)  % nochange 
         set(handles.edit11,'string',num2str(sorted_Dis(his)))
                  set(handles.edit12,'string',num2str(Two_D_Distance(I_Dis(his))))
           
 
                  [XCOR,MSE]=Sim_index(giveGrpWaveform(Groups{I_Dis(his)}),DataForFeature{unitsToBeLooped(unit)});
                  set(handles.edit15,'string',num2str(mean(MSE)))
                  XCOR(isnan(XCOR))=[];
                  set(handles.edit17,'string',num2str((1-mean(XCOR))*100))

    % do overlay plot , currentUnit and closest history unit. 
    % wait for user input to do one of the following and break inner loop,
    % add to existing   
    uiwait
    if breakFlag==1
    break;
    end
    
    end
    end
    
end

function Locations=GroupsLocation()
global Groups store
Locations = zeros(numel(Groups),3);
Locations(1,:)= [1000 1000 1000];
for i=2:numel(Groups)
    listOfUnit = Groups{i};
    Locations(i,:)=mean(cell2mat(store.LocForFeature(listOfUnit)));
end

function Ch_Map=giveGrpCh_Map(oneGroup)
global GUI_axes_map AllData lookUpTable 
Ch_Map=zeros([size(GUI_axes_map),numel(oneGroup)]); 
for i=1:numel(oneGroup)
Ch_Map(:,:,i)=AllData{lookUpTable(oneGroup(i),1)}.Ch_Map;
end
Ch_Map=max(Ch_Map,[],3);

%global AllData DataForFeature MaskFeatureAll lookUpTable store TimeFeatureAll LocForFeature
function GrpMask = giveGrpMask(oneGroup)
global  MaskFeatureAll 
GrpMask=min(cell2mat(MaskFeatureAll(oneGroup)));

function GrpLocation = giveGrpLocation(oneGroup)
global LocForFeature
GrpLocation=mean(cell2mat(LocForFeature(oneGroup)));


function GrptimeFeatures = giveGrpTime(oneGroup)
global  TimeFeatureAll 

countsSave=zeros([size(0:1:50),numel(oneGroup)]); 
counts_MSave=zeros([size(0:20:1000),numel(oneGroup)]); 
counts_LSave=zeros([size(0:1000:20000),numel(oneGroup)]); 
Max_Ins_5s_FRSave = zeros(numel(oneGroup),1);
FR_avgSave = zeros(numel(oneGroup),1);

for unit=1:numel(oneGroup)
FiringTimeForThisUnit = TimeFeatureAll{oneGroup(unit)};
ISI_in_MS_bins = diff(double(FiringTimeForThisUnit))/20; % Intervals in miliseconds. 
[counts] = histc(ISI_in_MS_bins,0:1:50);
[counts_M] = histc(ISI_in_MS_bins,0:20:1000);
[counts_L] = histc(ISI_in_MS_bins,0:1000:20000);
Sec_bins = double(FiringTimeForThisUnit)/20E3;
[counts_5s_1] = histc(Sec_bins,0:5:720);
[counts_5s_2] = histc(Sec_bins,2.5:5:722.5);
Max_Ins_5s_FR = max(max(counts_5s_1),max(counts_5s_2))/5;
FR_avg = numel(FiringTimeForThisUnit)/720;


countsSave(:,:,unit)=counts; 
counts_MSave(:,:,unit)=counts_M; 
counts_LSave(:,:,unit)=counts_L;  
Max_Ins_5s_FRSave(unit)=Max_Ins_5s_FR;
FR_avgSave (unit)= FR_avg;
end

GrptimeFeatures.counts=mean(countsSave,3);
GrptimeFeatures.counts_M=mean(counts_MSave,3);
GrptimeFeatures.counts_L=mean(counts_LSave,3);
GrptimeFeatures.Max_Ins_5s_FR=mean(Max_Ins_5s_FRSave);
GrptimeFeatures.FR_avg=mean(FR_avgSave);

function GrpWaveformFinal  = giveGrpWaveform(oneGroup)
% 1. get collective dataStructure
% 2. for each channel, remove those without datapoints. 
% 3. average the remaining signals directly
global DataForFeature 

GrpWaveform=zeros([size(DataForFeature{1}),numel(oneGroup)]); 
for i=1:numel(oneGroup)
GrpWaveform(:,:,i)=DataForFeature{oneGroup(i)};
end
for ch=1:32
OneChAllTime = reshape(GrpWaveform(ch,:,:),size(GrpWaveform,2),size(GrpWaveform,3)); % should be numberOfTimes - by - 41 element. after transpost 41,-by number of times 
OneChAllTime(:,min(OneChAllTime)==0)=[];
if size(OneChAllTime,2)==1
OneChAllTime(:,2)=OneChAllTime(:,1);
end
if ~isempty(OneChAllTime)
GrpWaveformFinal(ch,:)=mean(OneChAllTime');
else
GrpWaveformFinal(ch,:)=zeros(size(GrpWaveform(1,:,1)));
    
end
end

function [XCOR,MSE]=Sim_index(W1,W2)
p2p_1 = max(W1')-min(W1');
p2p_2 = max(W2')-min(W2');

 [~,I] = sort(p2p_1,'descend');
chA = I(1:4);
 [~,I] = sort(p2p_2,'descend');
chB = I(1:4);

selectedCh = union(chA,chB);

XCOR = zeros(size(selectedCh));
MSE = zeros(size(selectedCh));
for i = 1:numel(selectedCh)

  XCOR(i)= max(xcorr(W1(i,:),W2(i,:),10,'coeff'));
  MSE (i) = mean((W1(i,:)-W2(i,:)).^2);

end


function plot_Group(handles,unitsInGrp,colors,oldFlag)
global GUI_axes_map sorted_Dis his MSEList 
GUI_axes_map =[1	2	3	4	5	6	7	8
9	10	11	12	13	14	15	16
17	18	19	20	21	22	23	24
25	26	27	28	29	30	31	32
];
global AllData DataForFeature MaskFeatureAll lookUpTable store TimeFeatureAll LocForFeature
AllData=store.AllData;  DataForFeature=store.DataForFeature;
%MaskFeatureAll=store.MaskFeatureAll; 
lookUpTable=store.lookUpTable; TimeFeatureAll=store.TimeFeatureAll;
LocForFeature=store.LocForFeature;

t=[-30/20e3:1/20e3:10/20e3]*1000;
GrpCh_Map=giveGrpCh_Map(unitsInGrp);
%GrpMask = giveGrpMask(unitsInGrp);
GrptimeFeatures = giveGrpTime(unitsInGrp);
GrpWaveformFinal  = giveGrpWaveform(unitsInGrp);
GrpLocation = giveGrpLocation(unitsInGrp);


for row=1:4
            for col=1:8
        % we need compressed maskFeature, Compressed
        % DataForFeature,newChannelMap(finished), newTimingInfo,
        % newLocationInfo(finished) 
                
                    
                    temp=GrpCh_Map(row,col)-15;
                    if GrpCh_Map(row,col)>0 % && GrpMask (temp)>0.01
                        eval(['axes(handles.axes' num2str(GUI_axes_map(row,col)) ');' ]    );



                        plot(t,GrpWaveformFinal(temp,:),'Color',colors{1},'LineWidth',2)
                        set(gca,'xticklabel',[])
                        hold on;
                    end
                    axis([t(1) t(end) -175 25])
                    
                    
             
            end
end

counts = GrptimeFeatures.counts;
counts_M = GrptimeFeatures.counts_M;
counts_L = GrptimeFeatures.counts_L;
Max_Ins_5s_FR = GrptimeFeatures.Max_Ins_5s_FR;
FR_avg = GrptimeFeatures.FR_avg;
unit=1;

axes(handles.axes33)
if ~isempty(counts)
stairs(0.5:1:50.5,counts/max(counts),'Color',colors{unit})
xlim([0 51])
end
hold on;

axes(handles.axes34)
if ~isempty(counts_M)
stairs(10:20:1010,counts_M/max(counts_M),'Color',colors{unit})
xlim([0 1020])
end
hold on;        

axes(handles.axes35)
if ~isempty(counts_L)
stairs(0.5:1:20.5,counts_L/max(counts_L),'Color',colors{unit})
axis([0 21 0 1])
end
hold on;
if oldFlag==2
set(handles.edit6,'string',num2str(Max_Ins_5s_FR),'ForegroundColor','r')
set(handles.edit7,'string',num2str(FR_avg),'ForegroundColor','r')

else if oldFlag ==1
    set(handles.edit8,'string',num2str(Max_Ins_5s_FR),'ForegroundColor','b')
set(handles.edit9,'string',num2str(FR_avg),'ForegroundColor','b')
    end
end


%%
 axes(handles.axes37)
%  Ox=12.5;
%         Oy=12.5;
%         gap=5;
%         Gw=25;
%         Gh=25;
        
Ox=15;    
Oy=15;
gap=20;
        Gw=30;
        Gh=30;        

        numGx = 8;
        numGy = 4;
        Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
        Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
        [X,Y] = meshgrid(Gx,Gy);
    
    
    for r=1:4
        for c=1:8
%         scatter(X(r,c),Y(r,c),22000,[0 0 0],'s','LineWidth',2)
        pAx = X(r,c)-12.5;
        pAy = Y(r,c)-12.5;
        pCx = X(r,c)+12.5;
        pCy = Y(r,c)+12.5;
        plot([pAx pCx],[pAy pAy],'k','LineWidth',1)
        hold on
        plot([pCx pCx],[pAy pCy],'k','LineWidth',1)
        hold on 
        plot([pAx pCx],[pCy pCy],'k','LineWidth',1)
        hold on 
        plot([pAx pAx],[pAy pCy],'k','LineWidth',1)
        end
    end
    hold on 
    if oldFlag==2 % new unit
     scatter(GrpLocation(1),GrpLocation(2),30+(GrpLocation(3))^2,'+','LineWidth',2,'MarkerEdgeColor','r')

    else if oldFlag==1 % old unit
     scatter(GrpLocation(1),GrpLocation(2),30+(GrpLocation(3))^2,'+','LineWidth',2,'MarkerEdgeColor','b')
        end
    end
axis([-40 440 -40 200])
  if ~isempty(sorted_Dis) && ~isempty(his)  
 axes(handles.axes38)
  plot(1:10,sorted_Dis(1:10))
  hold on;
  scatter(his,sorted_Dis(his),20,'filled','MarkerFaceColor','r')
  axis([0 11 0 50])
  end

  axes(handles.axes36)
  plot(1:10,MSEList)
  hold on;
  scatter(his,MSEList(his),20,'filled','MarkerFaceColor','r')
  axis([0 11 0 50])

   
% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
global path his
if path==2
    path=1;
elseif path==1
        path=2;
end
his=0;
clear_main(handles)
uiresume
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  Groups deleteBreak
deleteBreak =0;
while true
listbox2value=get(handles.listbox2,'Value');
listbox1value=get(handles.listbox1,'Value');
SelectedGroup = Groups{listbox1value};
unitToBeDeleted=SelectedGroup(listbox2value);
SelectedGroup(listbox2value)=[]; %  happens first.
    % path 2 
    if numel(SelectedGroup)>2
        % 1. plot the history unit
         plot_Group(handles,SelectedGroup,cellstr('b'),1)  % some function plotting all units in a group.
         hold on;
         plot_unit(handles,unitToBeDeleted,cellstr('r'),2)  % nochange 
    else
        plot_unit(handles,SelectedGroup(1),cellstr('b'),1)  % nochange 
         hold on;
         plot_unit(handles,unitToBeDeleted,cellstr('r'),2)  % nochange 
    end

    % do overlay plot , currentUnit and closest history unit. 
    % wait for user input to do one of the following and break inner loop,
    % add to existing   
    uiwait
  if deleteBreak ==1
      break
  end
    
end


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global deleteBreak
deleteBreak =0;
clear_main(handles)
uiresume


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Groups
numberOfUnit = cellfun(@numel, Groups);
figure
histogram(numberOfUnit,'BinMethod','integers')


% --- Executes on button press in pushbutton32.
function pushbutton32_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global AloneUnitList unitsToBeLooped his LocForFeature
% unitsToBeLooped(unit)
listbox5value=get(handles.listbox5,'Value');

SelectedUnit= unitsToBeLooped(listbox5value);

% set(handles.edit15,'String',SelectedUnit);
% not considering only one unit left alone

AloneUnitList(AloneUnitList == SelectedUnit)=[];
LocForFeature{SelectedUnit}=[1000 1000 1000];

unitsToBeLooped(listbox5value)=[];

if listbox5value>numel(unitsToBeLooped)
listbox5value=listbox5value-1;
end
set(handles.listbox5,'Value',listbox5value);
set(handles.listbox5,'String',displayGroups(num2cell(unitsToBeLooped)))

if numel(unitsToBeLooped)<1
set(handles.listbox5, 'String', '');
end


clear_main(handles);
his=0;
% need to determine ending condition for this session. 
uiresume


% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global store Groups path
global AloneUnitList lookUpTable
global his unit unitsToBeLooped I_Dis breakFlag% unitsToBeLoopedIsAvailable
global sessionNum 
sessionNum = str2num(get(handles.edit16,'String'));
% unitsToBeLooped = find(store.lookUpTable(:,1)==sessionNum);% loading the second session. 
% can display them in another list box, 
% unitsToBeLoopedIsAvailable=true(size(unitsToBeLooped));
unit=1;
path=2; % compare with existing groups first
if numel(Groups)==1
    path=1;
end
while unit<=numel(unitsToBeLooped) && unit>=1
    breakFlag=0;
    set(handles.listbox5,'String',displayGroups(num2cell(lookUpTable(unitsToBeLooped,3))))
    set(handles.listbox5,'Value',unit);
    set(handles.edit13,'String',num2str(unitsToBeLooped(unit)));
    
    if path ==1 
    historyUnitLocations = cell2mat(store.LocForFeature(AloneUnitList));
    thisUnitLocation = store.LocForFeature{unitsToBeLooped(unit)};
    Distance = pdist2(thisUnitLocation,historyUnitLocations);
    Two_D_Distance = pdist2(thisUnitLocation(1:2),historyUnitLocations(:,1:2));
    [sorted_Dis,I_Dis] = sort(Distance);
    his = 1;
    while his<=numel(AloneUnitList) && his>=1
            set(handles.listbox4,'String',displayGroups(num2cell(lookUpTable(AloneUnitList(I_Dis),3))))
            set(handles.listbox4,'Value',his);
            set(handles.edit14,'String',num2str(AloneUnitList(I_Dis(his))));
    
        % 1. plot the history unit
         %plot_unit_ori(handles,AloneUnitList(I_Dis(his)),cellstr('b'),1)  
        % hold on;
        plot_unit_ori(handles,unitsToBeLooped(unit),cellstr('r'),2) 
         set(handles.edit11,'string',num2str(sorted_Dis(his)))
                  set(handles.edit12,'string',num2str(Two_D_Distance(I_Dis(his))))

    % do overlay plot , currentUnit and closest history unit. 
    % wait for user input to do one of the following and break inner loop,
    % add to existing   
    uiwait
    if breakFlag==1
    break;
    end
    end
    
    end
    
    if path ==2
    % path 2 
    historyGrpLocations = GroupsLocation();% some function given Groups variable, return locations.;
    thisUnitLocation = store.LocForFeature{unitsToBeLooped(unit)};
    Distance = pdist2(thisUnitLocation,historyGrpLocations);
    Two_D_Distance = pdist2(thisUnitLocation(1:2),historyGrpLocations(:,1:2));
    [sorted_Dis,I_Dis] = sort(Distance);
    his = 1;
    while his<=numel(Groups) && his>=1
            set(handles.listbox1,'String',displayGroups_special(Groups(I_Dis))) 
            set(handles.listbox1,'Value',his);
            set(handles.edit14,'String',num2str(Groups{I_Dis(his)}));
    
        % 1. plot the history unit
         plot_Group(handles,Groups{I_Dis(his)},cellstr('b'),1)  % some function plotting all units in a group.
         hold on;
         plot_unit(handles,unitsToBeLooped(unit),cellstr('r'),2)  % nochange 
         set(handles.edit11,'string',num2str(sorted_Dis(his)))
                  set(handles.edit12,'string',num2str(Two_D_Distance(I_Dis(his))))

    % do overlay plot , currentUnit and closest history unit. 
    % wait for user input to do one of the following and break inner loop,
    % add to existing   
    uiwait
    if breakFlag==1
    break;
    end
    
    end
    end
    
end

function plot_unit_ori(handles,unitNums,colors,oldFlag)
global GUI_axes_map
color=jet(120);
color=color(5:104,:);
firstRow=[0 0 0];
color=[firstRow;color];
GUI_axes_map =[1	2	3	4	5	6	7	8
9	10	11	12	13	14	15	16
17	18	19	20	21	22	23	24
25	26	27	28	29	30	31	32
];
t=[-30/20e3:1/20e3:10/20e3]*1000;
global AllData DataForFeature MaskFeatureAll lookUpTable store TimeFeatureAll LocForFeature StdForFeature unitsToBeLooped
AllData=store.AllData;  DataForFeature=store.DataForFeature;
% MaskFeatureAll=store.MaskFeatureAll; 
lookUpTable=store.lookUpTable; TimeFeatureAll=store.TimeFeatureAll;
LocForFeature=store.LocForFeature; StdForFeature=store.StdForFeature; ValleyAcrossTime=store.ValleyAcrossTime;
for row=1:4
            for col=1:8
        
                for unit=1:numel(unitNums)
                    Ch_Map=AllData{lookUpTable(unitNums(unit),1)}.Ch_Map;
                    temp=Ch_Map(row,col)-15;

                    if Ch_Map(row,col)>0 
                    
                        eval(['axes(handles.axes' num2str(GUI_axes_map(row,col)) ');' ]    );
                  errorbar(t,DataForFeature{unitNums(unit)}(temp,:),StdForFeature{unitNums(unit)}(temp,:),'Color',[0.9 0.9 0.9])
hold on
width = abs(min(DataForFeature{unitNums(unit)}(temp,:))/ max(StdForFeature{unitNums(unit)}(temp,:)))+0.1;
%                 plot(t,DataForFeature{unitNums(unit)}(temp,:),'Color',color(floor(MaskFeatureAll{unitNums(unit)}(temp)*100)+1,:),'LineWidth',width)
                                  plot(t,DataForFeature{unitNums(unit)}(temp,:),'Color','r','LineWidth',width)
                        
set(gca,'xticklabel',[])
                        hold on;
                    end
                    axis([t(1) t(end) -200 50])
                    
                    
                end
            end
end
 for unit=1:numel(unitNums)
 FiringTimeForThisUnit = TimeFeatureAll{unitNums(unit)};
ISI_in_MS_bins = diff(double(FiringTimeForThisUnit))/20; % Intervals in miliseconds. 
[counts] = histc(ISI_in_MS_bins,0:1:50);
[counts_M] = histc(ISI_in_MS_bins,0:20:1000);
[counts_L] = histc(ISI_in_MS_bins,0:1000:20000);
Sec_bins = double(FiringTimeForThisUnit)/20E3;
[counts_5s_1] = histc(Sec_bins,0:5:720);
[counts_5s_2] = histc(Sec_bins,2.5:5:722.5);
Max_Ins_5s_FR = max(max(counts_5s_1),max(counts_5s_2))/5;
FR_avg = numel(FiringTimeForThisUnit)/720;


axes(handles.axes33)
if ~isempty(counts)
stairs(0.5:1:50.5,counts,'Color',colors{unit})
xlim([0 51])
end
hold on;

axes(handles.axes34)
if ~isempty(counts_M)
stairs(10:20:1010,counts_M,'Color',colors{unit})
xlim([0 1020])
end
hold on;        

axes(handles.axes35)
if ~isempty(counts_L)
stairs(0.5:1:20.5,counts_L,'Color',colors{unit})
axis([0 21 0 1])
end
hold on;
if oldFlag==2
set(handles.edit6,'string',num2str(Max_Ins_5s_FR),'ForegroundColor','r')
set(handles.edit7,'string',num2str(FR_avg),'ForegroundColor','r')

else if unit ==1
    set(handles.edit8,'string',num2str(Max_Ins_5s_FR),'ForegroundColor','b')
set(handles.edit9,'string',num2str(FR_avg),'ForegroundColor','b')

    end
end

 end
%%
 axes(handles.axes37)
%  Ox=12.5;
%         Oy=12.5;
%         gap=5;
%         Gw=25;
%         Gh=25;
        
Ox=15;    
Oy=15;
gap=20;
        Gw=30;
        Gh=30;        

        numGx = 8;
        numGy = 4;
        Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
        Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
        [X,Y] = meshgrid(Gx,Gy);
    
    
    for r=1:4
        for c=1:8
%         scatter(X(r,c),Y(r,c),22000,[0 0 0],'s','LineWidth',2)
        pAx = X(r,c)-12.5;
        pAy = Y(r,c)-12.5;
        pCx = X(r,c)+12.5;
        pCy = Y(r,c)+12.5;
        plot([pAx pCx],[pAy pAy],'k','LineWidth',1)
        hold on
        plot([pCx pCx],[pAy pCy],'k','LineWidth',1)
        hold on 
        plot([pAx pCx],[pCy pCy],'k','LineWidth',1)
        hold on 
        plot([pAx pAx],[pAy pCy],'k','LineWidth',1)
        end
    end
    hold on 
    
    LocMat = cell2mat(LocForFeature(unitsToBeLooped));

       scatter(LocMat(:,1),LocMat(:,2),LocMat(:,3).^2+30,'+','LineWidth',0.5,'MarkerEdgeColor','b')
       hold on;
    if oldFlag==2 % new unit
     scatter(LocForFeature{unitNums}(1),LocForFeature{unitNums}(2),30+(LocForFeature{unitNums}(3))^2,'+','LineWidth',0.5,'MarkerEdgeColor','r')

    else if oldFlag==1 % old unit
     scatter(LocForFeature{unitNums}(1),LocForFeature{unitNums}(2),30+(LocForFeature{unitNums}(3))^2,'+','LineWidth',0.5,'MarkerEdgeColor','b')
        end
    end
    
axis([-40 440 -40 200])
%  axis([-40 250 -40 150])

 axes(handles.axes36)
  scatter(double(FiringTimeForThisUnit)/20E3/60,ValleyAcrossTime{unitNums},5*ones(size(ValleyAcrossTime{unitNums})),'+')



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
