function task()
clear all;
close all;
clc

global g_modify;
prompt = {'Enter the number of jobs to be scheduled  : '};
dlg_title = 'No. of Jobs';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines);
n = str2num(answer{1});

%% Input starting time array
dat{n,2} = [];

 for i=1:n
    for j=1:2
        dat{i,j} = '';
    end
 end

columnname = {};
columneditable = {};
columnformat = {};
rowname = {};

columnname{1} = 'Starting Time';
columnformat{1} = 'numeric';


columnname{2} = 'Ending Time';
columnformat{2} = 'numeric';


for i=1:n
   rowname{i} = strcat('Job ',num2str(i)); 
end

P = figure('Name','Task Scheduling','NumberTitle','off','Position',[50 50 230 400] );

t = uitable('Parent',P,...
            'Position',[0 0 230 320],'Data', dat,... 
            'ColumnName', columnname,...
            'ColumnFormat', columnformat,...
            'ColumnEditable', true,...
            'RowName',rowname, ...
            'FontSize',12,'ForegroundColor','k', ...
            'FontName','Comic Sans MS');

% ADDING A BUTTON AT THE BOTTOM

b = uicontrol('Parent',P,...
            'Style','Pushbutton',... 
            'Units','points', ...
            'Callback',@Pushbutton1_Callback,...
         'Position',[40 250 83.1724 30.4138], ...
         'String','Click To Schedule', ...
         'Tag','checkbox1' ); 
        
 P;
 function Pushbutton1_Callback(hObject,eventdata)
if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take appropriate action
    g = get(t,'Data');
    
    t = zeros(n,2);
    for i=1:n
        for j=1:2
            t(i,j) = str2num(g{i,j});
        end
    end 
    
    
    close all;
    TASK(t)     % Fucntion Call         
end
end  % end of Pushbutton Callback    

%% TASK function starts
function TASK(t)
% Creating the index array
index = char(1,n);
d = 'A';
for i = 1:n
	index(i) = d;
	d = d+1;
end

x = zeros(1,n);
% Check if end time is less than start time
for i = 1:n
    temp1 = t(i,1);
	temp2 = t(i,2);
	[temp1 temp2] = check(temp1,temp2,i,t);	
          t(i,1) = temp1;
          t(i,2) = temp2;
    	x(i) = t(i,2);
end

for i=1:n
    for j=1:2
        g_modify{i,j} = num2str(t(i,j));
    end
end

% Sorting by start time starts here
for i = 1:n
	for j = 1:(n-1)
		if t(j,1) > t(j+1,1) 
			temp     = t(j,1);
			t(j,1)   = t(j+1,1);
			t(j+1,1) = temp;
			
			temp     = t(j,2);
			t(j,2)   = t(j+1,2);
			t(j+1,2) = temp;
			
			temp = index(j);
			index(j) = index(j+1);
			index(j+1) = temp; 
		end
	end
end

maximum = max(x);
final = char(100,maximum);
limit = zeros(1,100);
m = 1;  %% Machines used
p = 1;  %% symbols to distinguish jobs
count = n; % No. of jobs left
machines = 0;
mac = 0;
i = 1;

for f = 1:10
	for d = 1:maximum
		final(f,d) = num2str(0);
	end
end
count_index = 1;

% Task allotment process begins here
while (1)	
		if t(i,1) >= limit(m)
			for j = t(i,1):t(i,2)
				if final(m,j) == num2str(0) 
					final(m,j) = index(count_index) ;
				else
					final(m,j) = '*';
				end	
			end
			count_index = count_index +1;
			limit(m) = t(i,2);       %% New limit for the machine
			p = p+1;                 %% Next symbol
			m = 1;                   %% Reset m from start 
			mac = mac + 1;
			count = count - 1;
			if count <=0 
				final_print(final,mac,maximum,index,machines);
				return
			end
			i = i + 1;
		else
			m = m+1;
        end
end

end
%% Sub function check 
function [temp1 temp2] = check(temp1,temp2,i,t)	
	if temp2 < temp1
		str = strcat('Job ',num2str(i));
        waitfor(msgbox('Ending time is less than starting time !! Re-Enter the value',str,'warn'));
        
        prompt = {'Enter the Starting Time for above Job : ','Enter the Ending Time for above Job : '};
        dlg_title = str;
        num_lines = 1;
        def = {'0','0'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
		temp1 = str2num(answer{1});
        temp2 = str2num(answer{2});
		[temp1 temp2] = check(temp1,temp2,i,t);
	end
end

%% Final Print Function 	
function final_print(final,mac,maximum,index,machines)
	d = 0;
	for g=1:mac
 		for  r = 1:maximum 
			if final(g,r) == num2str(0)
				d = d + 1; 
			end
		end	
		
		if d ~= maximum
			machines = machines + 1;
			d = 0;
		end
		d = 0;
	end
	
	d = 'A';
	for i = 1:n
		index(i) = d;
		d = d+1;
	end
	
    new = {};

    for g=1:machines+1
        check = maximum;
        for r=1:maximum
            if final(g,r)=='0'
                check = check - 1;
            end    
        end    

        if check ==0
           machines = machines-1;
        end    
    end 
    
    rowname1 = {};
    columnname1 = {};
    for g=1:machines+1
 		for  r = 1:maximum 
			new{g,r} = final(g,r);		
        end	
        rowname1{g} = strcat('Machine ',num2str(g));
    end
    
    for g=1:maximum
        columnname1{g} = strcat('Hour ',num2str(g));
    end

    close all;
    T = figure('Position',[100 100 600 500],'Name','Job Schedule','NumberTitle','off',...
                'DeleteFcn',@Figure_Close_Callback);  
    
    y1 = uicontrol('Parent',T,...
    'Style','text',...
    'Position', [20 470 560 30], ...
    'FontSize',15,'ForegroundColor','r', ...
    'FontWeight','bold', ...
    'FontName','Times New Roman', ...
    'String','NOTE :');

    y2 = uicontrol('Parent',T,...
    'Style','text',...
    'Position', [20 440 560 30], ...
    'FontSize',14,'ForegroundColor','b', ...
    'FontName','Times New Roman', ...
    'String','0 = Free, A = Job1, B = Job 2, C= Job 3 ... and so on');

    y3 = uicontrol('Parent',T,...
    'Style','text',...
    'Position', [20 390 560 50], ...
    'FontSize',14,'ForegroundColor','b', ...
    'FontName','Times New Roman', ...
    'String',' * indicates that the neighbouring two jobs share the same start and finish times');

    v = uitable('Parent',T,'Position',[0 0 600 370],'Data', new,...
                'FontSize',14,'ForegroundColor','k', ...
                 'FontName','Times New Roman', 'FontWeight','Bold', ...
                'ColumnEditable', false,'Rowname',rowname1,'Columnname',columnname1);
    T;
end

function Figure_Close_Callback(hObject,eventdata)
        %          Construct a questdlg with two options
            choice = questdlg('Would you like to solve a new problem ?', ...
             'Solution Achieved. Try another ?', ...
             'Yes','No','View/Modify Input Data','No');
            % Handle response
            switch choice
                case 'Yes'
                    close all; clear all;
                    start();
                    
                case 'No'
                    return  
                    
                case 'View/Modify Input Data'   
                   task_modify(g_modify,n);
            end       
end

end
	
