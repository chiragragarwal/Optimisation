function floyd_warshall_modify(dist_modify,n,q)
% Floyd - Warshall Algorithm to find the shortest path
close all;
clc
global adj;
global dist;
global edge;
global dist_given;
global dist_new;

temp3 = n;
temp2 = q;
[temp3 temp2] = check1(temp3,temp2);

n = temp3;
q = temp2;   

columnname = {};
columnformat = {};
rowname = {};

for i=1:q
    columnformat{i} = 'numeric';
end    

columnname{1} = 'Node 1' ;
columnname{2} = 'Node 2' ;
columnname{3} = 'Distance' ;

for i=1:q
    rowname{i} = strcat('Edge ',num2str(i));    
end

dist = dist_modify;

P = figure('Name','Enter the Edge Information','NumberTitle','off','Position',[50 50 320 400] );
t = uitable('Parent',P,...
            'Position',[0 0 320 320],'Data', dist,...
            'ColumnName', columnname,...
            'ColumnEditable', true,...
            'FontSize',12,'ForegroundColor','k', ...
            'FontName','Comic Sans MS', ...
            'RowName',rowname);


b = uicontrol('Parent',P,...
            'Style','Pushbutton',... 
            'Units','points', ...
            'Callback',@Pushbutton1_Callback,...
         'Position',[120 250 83.1724 30.4138], ...
         'String','Click To Solve', ...
         'Tag','checkbox1' );
     
 str1 = strcat('Vertices = ',num2str(n));
 
 y1 = uicontrol('Parent',P,...
    'Style','text',...
    'Position', [30 340 110 30.4138], ...
    'FontSize',10,'ForegroundColor','k', ...
    'FontWeight','bold', ...
    'FontName','Arial Black', ...
    'String',str1);    
 P;

function Pushbutton1_Callback(hObject,eventdata)
if (get(hObject,'Value') == get(hObject,'Max'))

g = get(t,'Data');

 for i=1:q
        for j=1:2
            temp1 = g{i,j};
            temp1 = check(temp1,i,j,g);
            g{i,j} = temp1;
        end
 end

 dist_modify = g;
 
new = {};
for i=1:q
    for j=1:3
        new{i,j} = str2num(g{i,j});
    end
end

A = {};
adj = {};
for i=1:n
    for j=1:n
        A{i,j} = 0;
        adj{i,j}=0;
    end
end

for i=1:q  
        A{new{i,1},new{i,2}} = new{i,3};
        adj{new{i,1},new{i,2}} = 1;
        A{new{i,2},new{i,1}} = new{i,3};
        adj{new{i,2},new{i,1}} = 1;
end

dist = A;
dist_given = A;

for k=1:n
	for i=1:n
		for j=1:n
			if (dist{i,k}*dist{k,j} ~= 0) && (i ~= j)
			
				if (dist{i,k} + dist{k,j} < dist{i,j}) || (dist{i,j} == 0)
				
					dist{i,j} = dist{i,k} + dist{k,j};
			
                end	
               
            end            
		end
	end
end

plot_graph(adj,n,dist_given,dist)  

end
end

%% Final Plot Function
function plot_graph(dist,n,weight,dist_new) 
close all;
Q = figure('Name','Floyd Warshall Shortest Path Matrix','NumberTitle','off','Position',[50 50 700 300],...
            'DeleteFcn',@Figure_Close_Callback);  
v = uitable('Parent',Q,'Position',[0 0 700 300],'Data', dist_new,...
            'ColumnEditable', false);         
v;    
theta = 2*pi/n;

vertices = [];
r = 450;
c = 1;

for x=0:theta:(2*pi)
    xcord = r*cos(x);
    ycord = r*sin(x);
    vertices(c) = (xcord) + 1i*(ycord);
   c = c+1;
end    

for i=1:n
    X(i) = real(vertices(i));
    Y(i) = imag(vertices(i));
end

T = figure('Name','Floyd Warshall : Given Graph','NumberTitle','off','Position',[50 50 1200 600],...
            'DeleteFcn',@Figure_Close_Callback);
plot(real(vertices),imag(vertices),'ob','LineWidth',4,...
                'MarkerSize',10)        

for i=1:n
    if X(i)>=0 && Y(i)>=0
        text(X(i)+20,Y(i)+20,num2str(i),'FontSize',18,'Color','r'); 
    elseif X(i)<=0 && Y(i)>=0
        text(X(i)-30,Y(i)+30,num2str(i),'FontSize',18,'Color','r');
    elseif X(i)<=0 && Y(i)<=0
        text(X(i)-30,Y(i)-30,num2str(i),'FontSize',18,'Color','r');
    elseif X(i)>=0 && Y(i)<=0    
        text(X(i)+20,Y(i)-20,num2str(i),'FontSize',18,'Color','r');
    end    
end     

% Displaying weights
for i=1:n
    for j=1:n
          g = [X(i) X(j)];
          h = [Y(i) Y(j)];
          if X(i) <= X(j) 
              Wx(i,j) = X(i) + abs(X(j)-X(i))/2;
              if Y(i) <= Y(j)
                  Wy(i,j) = Y(i) + abs(Y(j)-Y(i))/2;
              else
                  Wy(i,j) = Y(j) + abs(Y(i)-Y(j))/2;
              end
          else 
               Wx(i,j) = X(j) + abs(X(i)-X(j))/2;
              if Y(i) <= Y(j)
                  Wy(i,j) = Y(i) + abs(Y(j)-Y(i))/2;
              else
                  Wy(i,j) = Y(j) + abs(Y(i)-Y(j))/2;
              end
          end    
          if dist{i,j} == 1
             line(g,h,'LineWidth',1.5)
          end   
    end
end

for i=1:n
    for j=1:n
         if X(i)~=X(j) && Y(i)~=Y(j)
             if weight{i,j} ~= 0
                text(Wx(i,j),Wy(i,j),num2str(weight{i,j}),'FontSize',18,'Color','m'); 
             end   
         end     
    end    
end    

end 
%% Check Function For Bounds
 function temp1 = check(temp1,i,j,g)
                if str2num(temp1) > n
                        str = strcat('Edge ',num2str(i)','--Node',num2str(j));
                        waitfor(msgbox('Vertex Number exceeds Total number of vertices ! Re-Enter the Value.',str,'warn'));
                        
                        prompt = {'Enter the Vertex Number :'};
                        dlg_title = str;
                        num_lines = 1;
                        def = {'1'};
                        answer = inputdlg(prompt,dlg_title,num_lines,def);
                        temp1 = answer{1};                     
                        temp1 = check(temp1,i,j,g);
                end            
 end     
%% Check Function for Max Edges
function [temp3 temp2] = check1(temp3,temp2)
                temp_new = temp3*(temp3-1)/2;
                if temp2 > temp_new
                        str = 'ERROR !';
                        str1 = strcat('Number of Edges cannot exceed --> ',num2str(temp_new),' <-- for -->',num2str(temp3),' <-- Vertices');
                        waitfor(msgbox(str1,str,'warn'));
                        
                        prompt = {'Enter the number of Vertices : ','Enter the number of Edges : '};
                        dlg_title = str;
                        num_lines = 1;
                        def = {'4','5'};
                        answer = inputdlg(prompt,dlg_title,num_lines,def);
                        temp3 = str2num(answer{1});    
                        temp2 = str2num(answer{2});
                    
                        [temp3 temp2] = check1(temp3,temp2);
                end            
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
                   floyd_warshall_modify(dist_modify,n,q);
            end       
end

end    % floyd
