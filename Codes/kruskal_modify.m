function kruskal_modify(edge_modify,n,q)
close all;
clc;

temp3 = n;
temp2 = q;
[temp3 temp2] = check1(temp3,temp2);

n = temp3;
q = temp2;

dist = edge_modify;
columnname = {};
columnformat = {};
rowname = {};

for i=1:q
    columnformat{i} = 'numeric';
end    

columnname{1} = 'Node 1' ;
columnname{2} = 'Node 2' ;
columnname{3} = 'Edge Weight' ;

for i=1:q
    rowname{i} = strcat('Edge ',num2str(i));    
end

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
    new = {};
    
      for i=1:q
        for j=1:2
            temp1 = g{i,j};
            temp1 = check(temp1,i,j,g);
            g{i,j} = temp1;
        end
      end
    
      edge_modify = g;
      
    for i=1:q
        for j=1:3
            new{i,j} = str2num(g{i,j});
        end
    end
    
    edge = new;

    kruskal_solve(edge,q,n)
    end  
end

%% Kruskal_Solve Function
function kruskal_solve(edge,q,n)
 father = (1:n);
 edge = edge' ;
 
% Sorting by weight begins here    
    for i=1:q-1
        for j=1:q-1
            if edge{3,j} > edge{3,j+1}
                temp(1) = edge{1,j};
                temp(2) = edge{2,j};
                temp(3) = edge{3,j};
                    edge{1,j} = edge{1,j+1};
                    edge{2,j} = edge{2,j+1};
                    edge{3,j} = edge{3,j+1};
                edge{1,j+1}=temp(1);
                edge{2,j+1}=temp(2);
                edge{3,j+1}=temp(3);
            end
        end
    end
% sorting ends here

% Kruskal Algorithm starts here	
   T = {};     
   k = 1;
   for i=1:q
       if father(edge{1,i}) ~= father(edge{2,i})
           s = father(edge{2,i});
           for j=1:n
               if father(j) == s
                   father(j) = father(edge{1,i});
               end
           end    
               T{1,k} = edge{1,i};
               T{2,k} = edge{2,i};
               T{3,k} = edge{3,i};
               k = k+1;
           
       end
   end  %% end of for loop

   adj_out = {};        %% Output adjacency Matrix
   
   for i=1:n            %% Initialising to zero
       for j=1:n
           adj_out{i,j} = 0;
       end
   end    
   
   for i=1:k-1
      adj_out{T{1,i},T{2,i}} = 1;
	  adj_out{T{2,i},T{1,i}} = 1; 
   end    
   
   plot_graph(adj_out,n)  
end   

%% FINAL PLOT_GRAPH function
function plot_graph(dist,n)

close all;    
theta = 2*pi/n;

vertices = [];
r = 450;
c = 1;

for x=0:theta:(2*pi)
    xcord = r*cos(x);
    w1 = xcord/power(10,-10);
    if abs(w1)<1
         xcord=0;  
    end   
    ycord = r*sin(x);
    w2 = ycord/power(10,-10);
    if abs(w2)<1
         ycord=0;  
    end 
    vertices(c) = (xcord) + 1i*(ycord);
   c = c+1;
end       

for i=1:n
    X(i) = real(vertices(i));
    Y(i) = imag(vertices(i));
end

T = figure('Name','Minimal Spanning Tree','NumberTitle','off','Position',[50 50 1200 600],...
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
            
for i=1:n
    for j=1:n
          g = [X(i) X(j)];
          h = [Y(i) Y(j)];
          if(dist{i,j} == 1)
             line(g,h,'LineWidth',1.5,'Color','k')
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
                   kruskal_modify(edge_modify,n,q);
            end       
end

end