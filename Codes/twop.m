function twop()
clear all
close all
clc

global c;
global phase_new;
prompt = {'Enter the number of variables : ','Enter the number of constraints :'};
dlg_title = 'Constraints/Variables';
num_lines = 1;
def = {'2','2'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
v = str2num(answer{1});
cons = str2num(answer{2});

global cons_org;
global v_org;
v_org = v;
global g_modify;

global daddy_dat;
global dat1;
global daddy_row;

daddy_dat = [];         %% Final solution table
dat1 = [];              %% each iteration tables
daddy_row = [];         %% Final rows

global rowname1;
rowname1 = [];

dat{cons+4,v+2} = [];

 for i=1:cons+4
    for j=1:v+2
        dat{i,j} = '';
    end
 end

 for i=1:v
     dat{cons+2,i} = '0';
     dat{cons+3,i} = 'Inf';
     dat{cons+4,i} = 'No';
 end
 
columnname = {};
columneditable = {};
columnformat = {};

rowname = {};

for i=1:v
    str1 = strcat('x',num2str(i));
    columnname{i} = str1 ;
   columnformat{i} = 'numeric';
   columneditable{i} = true;
end

columnname{i+1} = 'Enter <= , >= or = ';
columnformat{i+1} = {'<=' '>=' '='};
columneditable{i+1} = true;
columneditable{i+2} = true;
columnname{i+2}  = 'R.H.S.';
columnformat{i+2}= 'numeric';

rowname = {};
rowname{1} = 'Objective Function';

for i=2:cons+1
   rowname{i} = strcat('Constraint ',num2str(i-1));
 
end

rowname{i+1} = 'Lower Bound';
rowname{i+2} = 'Upper Bound';

rowname{i+3} = 'Unrestricted ?';

P = figure('Name','Two Phase Simplex Table Values','NumberTitle','off','Position',[50 50 800 400] );

t = uitable('Parent',P,...
            'Position',[0 0 800 300],'Data', dat,... 
            'ColumnName', columnname,...
            'ColumnFormat', columnformat,...
            'ColumnEditable', true,...
            'FontSize',10,'ForegroundColor','k', ...
            'FontName','Comic Sans MS', ...
            'RowName',rowname );

% ADDING A BUTTON AT THE BOTTOM

b = uicontrol('Parent',P,...
            'Style','Pushbutton',... 
            'Units','points', ...
            'Callback',@Pushbutton1_Callback,...
         'Position',[230 260 83.1724 30.4138], ...
         'String','Click To Solve', ...
         'Tag','checkbox1' );
     
% Create the button group.
h = uibuttongroup('Parent',P,'visible','on','Position',[0 0 83.1724 30.4138]);
% Create two radio buttons in the button group.
u0 = uicontrol('Style','Radio','String','Maximize',...
    'pos',[100 350 100 30],'parent',h,'HandleVisibility','on');
u1 = uicontrol('Style','Radio','String','Minimize',...
    'pos',[200 350 100 30],'parent',h,'HandleVisibility','on');
  set(h,'Visible','on');     
     
function Pushbutton1_Callback(hObject,eventdata)
if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take appropriate action
    g = get(t,'Data');
    g_modify = g;
            for i=1:cons+1
                for j=1:v
                    sz = size(g{i,j});
                    if sz(2)==0
                       g{i,j} = '0'; 
                    end    
                end
            end
      
             for i=2:cons+1
                    sz = size(g{i,v+2});
                    if sz(2)==0
                       g{i,v+2} = '0'; 
                    end      
             end
            
            for i=2:cons+1
                    sz = size(g{i,v+1});
                    if sz(2)==0
                       g{i,v+1} = '<='; 
                    end      
            end
                
    A = zeros(cons,v);
    for i=2:cons+1
        for j=1:v
            A(i-1,j) = str2num(g{i,j});
        end
    end
    
    b = zeros(1,cons);
    for i=2:cons+1
            b(i-1) = str2num(g{i,v+2});   
    end
   
    c = zeros(1,v);
    for j=1:v
            c(j) = str2num(g{1,j});   
    end
    
    inq = zeros(1,cons);
    k=1;
    
    for i=2:(cons+1)
        if strcmp('<=',g{i,v+1}) == 1
            inq(k) = -1;
            k = k+1;
        end
        
        if strcmp('=',g{i,v+1}) == 1
            inq(k) = 0;
            k = k+1;
        end
        
        if strcmp('>=',g{i,v+1}) == 1
            inq(k) = 1;
            k = k+1;
        end
    end
    cons_org = cons;
    
           % Check if lower bound is greater than upper bound
            for er=1:v
                temp1 = g{cons+2,er};
                temp2 = g{cons+3,er};
                   [temp1 temp2] = check(temp1,temp2,g,er);
                g{cons+2,er} = temp1;
                g{cons+3,er} = temp2;
            end
            
          % Check for upper bounds
            j1=0;
            for i=1:v
               if str2num(g{cons+3,i}) ~= Inf
                  temp = zeros(1,v);
                  temp(i) = 1;
                  A = vertcat(A,temp);
                  inq = horzcat(inq,-1);
                  b = horzcat(b,str2num(g{cons+3,i}));
                  j1 = j1+1;
               end    
            end    
         
            % Check for lower bounds
            j2=0;
            for i=1:v
               if str2num(g{cons+2,i}) ~= 0
                  temp = zeros(1,v);
                  temp(i) = 1;
                  A = vertcat(A,temp);
                  inq = horzcat(inq,1);
                  b = horzcat(b,str2num(g{cons+2,i}));
                  j2 = j2+1;
               end    
            end    
            
            cons = cons + j1 + j2;
            
            % Check for Unrestricted
            j = 0;
            U_indices = 0;
            U_indices_new = 0;
            U_id = 0;
            
            g_modify = g;
            for i=1:v
               if g{cons_org+4,i} == 'y'
                    if str2num(g{cons_org+3,i}) ~= Inf
                        waitfor(msgbox('A variable cannot be both unrestricted and have bounds simultaneously !','Error !','warn'));
                        twop_modify(g_modify,v_org,cons_org);
                    end    
                   Atemp = -A(:,i+j);
                   Apart1 = A(:,(1:i+j));
                   Apart2 = A(:,(i+j+1:v+j));
                   
                   Ctemp = -c(:,i+j);
                   Cpart1 = c(:,(1:i+j));
                   Cpart2 = c(:,(i+j+1:v+j));
                   
                   U_indices(i) = 1;
                   U_indices_new(i+j) = 1;
                   U_indices_new(i+j+1) = 2;
                   U_id(i+j) = i;
                   U_id(i+j+1) = i;
                   
                   A = horzcat(Apart1,Atemp,Apart2);
                   c = horzcat(Cpart1,Ctemp,Cpart2);
                   j = j+1;
               else
                   U_indices(i) = 0;
                   U_indices_new(i+j) = 0;
                   U_id(i+j) = i;
               end    
            end    
            U_indices = horzcat(U_indices,zeros(1,j));
            v = v+j;
            
 %% Check for maximize or minimize
 
 var = get(u0,'Value');
 if var == 1
    minimize = 0;
 else
     new = get(u1,'Value');
     if new == 1
         minimize = 1;
     end
 end       
    close all;
    twop_solve(c,b,A,inq,minimize,U_indices,U_indices_new,U_id)     % Function Call         
end
end  % end of Pushbutton Callback

%% Start function TWOPHASE_SOLVE
function twop_solve(c,b,A,inq,minimize,U_indices,U_indices_new,U_id)
    
global phase_check;
phase_check = 0;

global art_check;
art_check = 0;

global var_art_global;
var_art_global = 0;

global var_art;
var_art = 0;

global art_inf;
art_inf = 0;

global n;
global m;
global A2;
global c2;
global ent_var_index;
% global phase_new;
        [m,n] = size(A);
        cols = n;
        b = b' ;
        
        % Check for negative RHS, then reverse the inequality
        for i=1:m
            if b(i)<0
                inq(i) = -inq(i);
                b(i) = -b(i);
                for j=1:n
                    A(i,j) = -A(i,j);
                end
            end
        end
        global u;
        global f;   %% Final optimal solution
        global slack;
        global surp;
        global art;
        
        
        phase_new =1;
        f = 0;
        slack = 0;
        surp  = 0;
        art   = 0;
        u = 1;
        global art_indices;
        art_indices = [];           %% Columwise (Along the vertical length)
        global k_art;
        k_art = 1;                  %% count for art_indices
        
        M = zeros(m,3);         %% IDENTIFIER matrix     1-slack, 2-surplus, 3-artificial
        
        % Calculate number of slack, surplus and artificial variables
        for i=1:m
            switch inq(i)
                case -1                    %% <=
                    slack = slack + 1;
                    M(i,1) = 1;
                case 0                     %% =
                    art = art + 1;
                    M(i,3) = 1;
                case 1                       %% >=
                    art = art + 1;
                    surp = surp + 1;
                    M(i,2) = 1;
                    M(i,3) = 1;
            end
        end
        
        global total;
        total = n + slack + surp + art ;            %%  Total Columns
        for i=v+1:v+total
               U_id(i) = i; 
        end 
        
        % First loop for surplus variables
        for i=1:m
            if M(i,2)==1
                a = zeros(m,1);
                a(i) = -1;
                A(:,cols+1) =  a;
                c(:,cols+1) =  0;
                cols = cols+1;
            end
        end
        
        % Second loop for artificial and slack variables (giving identity matrix)
        for i=1:m
            if M(i,3)==1    % art
                art_indices(k_art) = n+surp+i;
                k_art = k_art + 1;
                
                a = zeros(m,1);
                a(i) = 1;
                A(:,cols+1) =  a;
                cols = cols+1;
                
            elseif M(i,1) == 1 %slack
                a = zeros(m,1);
                a(i) = 1;
                A(:,cols+1) =  a;
                c(:,cols+1) =  0;
                cols = cols+1;
            end
        end
        
        k_art = k_art-1;
        A1 = A;
        c1 = c;
        global iter;
        iter = 0;
        
        global iter_total;
        
        B_indices = (n+surp+1:n+slack+surp+art);          %% indices of Basic Variables
        V_indices = (1:n+surp);                          %% indices of Non-Basic Variables
        
        if minimize==1               %% For minization problem
            c1 = -c1 ;
            c = c1;
        end

A2 = A1;

% Initiating Phase 1 % ****************************************************     
z1 = zeros(1,total);
for i=1:total
   for j=1:k_art
       if i==art_indices(j)
           z1(i)=1;
           c(i) = 0;
       end
   end    
end
% Original value of objective function is stored in c
c1 = z1;
%% Since Phase 1 has minization problem
 c1 = -c1 ;
 
%% Begin While Loop
        while (1==1)
            iter = iter + 1;
            B = A1(:,B_indices);
            Aj = A1(:,V_indices);
            
            if u == -999
                return
            end
            
            Cb = c1(:,B_indices);
            Cj = c1(:,V_indices);
            f = Cb*(B\b);
            bstar = B\b;          %% B\b = inv(B)*b
            
            c2 = c1;
            
            cols = total;
            
            for i=1:cols
                size_B = size(find(B_indices == i));
                if size_B(2) ~= 0              %% IF BASIC VARIABLE
                    index_B = find(B_indices ==  i);
                    A2(:,i) = 0;
                    A2(index_B,i) = 1;
                    c2(1,i) = 0;
                else                           %% IF NON-BASIC VARIABLE
                    A2(:,i) = B\A1(:,i) ;
                    index_V = find(V_indices == i);
                    c2(i) = Cb*(B\Aj(:,index_V)) - Cj(index_V) ;
                end
            end
            
            if minimize == 0
                c2 = -c2;
                f  = -f;
            end
            
             %****************************** Storing values in table *******************************************

                rowname1{iter}{1} = 'z(min)';
            
            for i=2:m+1
                if B_indices(i-1) <= v      %% Decision Variable
                    if U_indices_new(B_indices(i-1)) == 0
                        rowname1{iter}{i} = strcat('x',num2str(U_id(B_indices(i-1))));
                    elseif U_indices_new(B_indices(i-1)) == 1
                        rowname1{iter}{i} = strcat('x',num2str(U_id(B_indices(i-1))),'+'); 
                    elseif U_indices_new(B_indices(i-1)) == 2
                        rowname1{iter}{i} = strcat('x',num2str(U_id(B_indices(i-1))),'-');     
                    end   
                else
                    if (B_indices(i-1) >= v+1) && (B_indices(i-1) <= v+surp)        %% Surplus Variable
                        rowname1{iter}{i} = strcat('Sx',num2str(B_indices(i-1)));
                    else    
                        var = size(find(art_indices == i-1+n+surp));
                        if var(2) == 0
                            rowname1{iter}{i} = strcat('sx',num2str(B_indices(i-1)));
                        else
                            rowname1{iter}{i} = strcat('Rx',num2str(B_indices(i-1)));
                        end    
                            
                    end    
                end
            end
            
            
            dat1{iter} = [];
            dat1{iter}{1,1} = '1';
            
            for i=2:m+1
                dat1{iter}{i,1} = '0';
            end
            
            for i=2:total+1
                if minimize == 1
                    dat1{iter}{1,i} = num2str(-c2(i-1));
                else
                    dat1{iter}{1,i} = num2str(c2(i-1));
                end    
            end
            
            dat1{iter}{1,total+2} = num2str(f);
            
            for i=1:total               %% Printing A2
                for j=1:m
                    dat1{iter}{j+1,i+1} = num2str(A2(j,i)) ;
                end
            end
            
            for i=1:m                 %% Printing bstar
                dat1{iter}{i+1,total+2} = num2str(bstar(i));
            end
                        
            %****************************** Done storing values in table *************************************************
            
            ent_var_mat = Cb*(B\Aj) - Cj;
            
            if min(ent_var_mat) >= 0                                      %% Optimality condition satisfied
                
                % Check for INFEASIBLE SOLUTION
                for i=1:m
                    for j=1:art
                        if (B_indices(i)==art_indices(j)) && (bstar(i)~=0)
                            waitfor(msgbox('The given problem has INFEASIBLE SOLUTION !'));
                            phase_check = 1;
                            final_print(f,minimize,A2,c2,bstar,m,n,v,cons,B_indices,iter,rowname1,slack,surp,art,total,art_indices,iter_total,phase_check,U_indices,U_indices_new,U_id)
                            return
                        end
                    end
                end
                
               if f~=0
                  waitfor(msgbox('The given problem has INFEASIBLE SOLUTION !'));
                  phase_check = 1;
                  final_print(f,minimize,A2,c2,bstar,m,n,v,cons,B_indices,iter,rowname1,slack,surp,art,total,art_indices,iter_total,phase_check,U_indices,U_indices_new,U_id)
                  return               
               end    
                final_print(f,minimize,A2,c2,bstar,m,n,v,cons,B_indices,iter,rowname1,slack,surp,art,total,art_indices,iter_total,phase_check,U_indices,U_indices_new,U_id)
                second_phase(A1,B,bstar,ent_var_index,B_indices,V_indices,m, u,k_art,art_indices,n,surp,phase_new,f,c,c1,minimize,b,Aj,art,A2,c2,rowname1,total,slack,iter,iter_total,daddy_dat,daddy_row,dat1,var_art,art_check,var_art_global,art_inf,U_indices,U_indices_new,U_id)         %% Call Second_phase function
                return
            else
                for i=1:total+2                 %% Printing BLANK ROW
                    dat1{iter}{m+2,i} = '' ;
                end
                
                rowname1{iter}{m+2} = strcat('******** Iteration No. ',num2str(iter+1),'  *********');
                
                ent_var_index = find(ent_var_mat == min(ent_var_mat));
                ent_var_index = V_indices(ent_var_index(1));                %% Entering Variable
                
                
                [A1 B bstar ent_var_index B_indices V_indices m u k_art art_indices phase_new rowname1 var_art art_check] = subroutine(A1,B,bstar,ent_var_index,B_indices,V_indices,m, u,k_art,art_indices,n,surp,phase_new,rowname1,iter_total,var_art,art_check,var_art_global);    %% Calling the SUBROUTINE function
            end
            
        end         %% End of while loop
end        
%% Final PRINT Phase 1 function
    function final_print(f,minimize,A2,c2,bstar,m,n,v,cons,B_indices,iter,rowname1,slack,surp,art,total,art_indices,iter_total,phase_check,U_indices,U_indices_new,U_id)    
       
          % ************** PRINTING OPTIMAL TABLE BEGINS HERE **********************************************************************
        T = figure('Name','Two Phase Method : PHASE I Table Values','NumberTitle','off','Position',[50 50 1200 600],...
            'DeleteFcn',@Figure_Close_Callback);
        columnname = {};
        columnformat = {};
        
        columnname{1} = 'z';    %% First column name is 'z'
      
        
         j=2;
        for i=2:v+1               %% Columnames for decision variables
            if U_indices(i-1) == 0
                str1 = strcat('x',num2str(i-1));
                columnname{j} = str1 ;
                j=j+1;
            else
                str1 = strcat('x',num2str(i-1),'+');
                str2 = strcat('x',num2str(i-1),'-');
                columnname{j} = str1 ;
                columnname{j+1} = str2 ;
                j = j+2;
            end    
        end
        
        for i=v+2:v+2+surp-1           %% Columnnames for surplus variables
            str1 = strcat('Sx',num2str(i-1));
            columnname{i} = str1 ;
        end
        
        for i=(v+surp+2):total+1
            var = size(find(art_indices == i-1));
            if var(2) == 0
                str1 = strcat('sx',num2str(i-1));
                columnname{i} = str1 ;
            else
                str1 = strcat('Rx',num2str(i-1));
                columnname{i} = str1 ;
            end    
        end
               
        columnname{total+2} = 'Solution';
       
        if mod(iter,2) == 0     %% If even iterations
            for i=1:2:iter
                xyz = vertcat(dat1{i},dat1{i+1});
                daddy_dat = vertcat(daddy_dat,xyz);
                
                lol = horzcat(rowname1{i},rowname1{i+1});
                daddy_row = horzcat(daddy_row,lol);
            end
            
        else
             if iter ~= 1
                for i=1:2:iter-1
                    xyz = vertcat(dat1{i},dat1{i+1});
                    daddy_dat = vertcat(daddy_dat,xyz);

                    lol = horzcat(rowname1{i},rowname1{i+1});
                    daddy_row = horzcat(daddy_row,lol);
                end
  
                daddy_dat = vertcat(daddy_dat,dat1{iter});
                daddy_row = horzcat(daddy_row,rowname1{iter});
             else
                daddy_dat = dat1{1};
                daddy_row = rowname1{1};    
             end      
        end
              
        q = uitable('Parent',T,...
            'Position',[0 0 1200 600],'Data', daddy_dat,...
            'ColumnName', columnname,...
            'ColumnEditable', false,...
            'FontSize',10,'ForegroundColor','k', ...
             'FontName','Comic Sans MS', ...
            'RowName',daddy_row );
        
       if phase_check == 1 
            T
       end     
        %*************** PRINTING OPTIMAL TABLE ENDS HERE **********************************************************************
        
        if phase_check == 1
            k = 1;
        else   
            return    
        end    
    end

%% ALTERNATE CHECK function
    function [B_indices B bstar iter A2 c2 f rowname1 iter_total dat1 U_indices U_indices_new U_id u] = alternate_check(A1,B, bstar, B_indices, V_indices, m,ent_var_mat,u,n,A2,minimize,iter,c2,f,b,c1,ent_var_index,art_indices,surp,rowname1,iter_total,total,dat1,U_indices,U_indices_new,U_id)
        var = size(ent_var_mat);
        length = var(2);
        
        for temp=1:length
            if ent_var_mat(temp)==0
                x = U_id(V_indices(temp));
                q = find(U_id(B_indices)==x);
                q = size(q);
                if q(2) == 0
                ent_var_index = V_indices(temp);
                
                waitfor(msgbox('Optimal Solution Found ! '));
                waitfor(msgbox('ALTERNATE SOLUTION EXISTS !! Click OK to proceed '));
                
                iter_total = iter_total + 1;
                
                for i=1:total+2                 %% Printing BLANK ROW
                    dat1{iter_total}{1,i} = '' ;
                end
                
                rowname1{iter_total}{1} = strcat('ALTERNATE OPTIMAL SOLUTION ');
                
                iter_total = iter_total + 1;
                
                [A1 B bstar ent_var_index B_indices V_indices m u art_indices] = subroutine_alternate(A1,B,bstar,ent_var_index,B_indices,V_indices,m,u,art_indices,n,surp);
                
                if u == -999
                    iter_total = iter_total - 1;
                     rowname1{iter_total}{1} = strcat('ALTERNATE UNBOUNDED SOLUTION ');
                    return
                end    
                
                B  = A1(:,B_indices);
                Aj = A1(:,V_indices);
                Cb = c1(:,B_indices);
                Cj = c1(:,V_indices);
                
                bstar = B\b;            %% B\b = inv(B)*b
                for y=1:m
                    w=bstar(y)/power(10,-10);
                    if abs(w)<1
                       bstar(y)=0;  
                    end                          
                end
               
                c2 = c1;
                
                f = Cb*(B\b);                                             %% Calculate Final Alternate Optimal Solution
                if minimize == 1
                    f = -f;
                end
                
                cols = total;
                
                for i=1:cols
                    size_B = size(find(B_indices == i));
                    if size_B(2) ~= 0              %% IF BASIC VARIABLE
                        index_B = find(B_indices ==  i);
                        A2(:,i) = 0;
                        A2(index_B,i) = 1;
                        c2(1,i) = 0;
                    else                           %% IF NON-BASIC VARIABLE
                        A2(:,i) = B\A1(:,i) ;
                        index_V = find(V_indices == i);
                        
                        c2(i) = Cb*(B\Aj(:,index_V)) - Cj(index_V) ;
                    end
                end   
                              
                %****************************** Storing values in table *******************************************               
                                
                if minimize == 1
                    rowname1{iter_total}{1} = 'z(min)';
                else
                    rowname1{iter_total}{1} = 'z(max)';
                end
                
                for i=2:m+1
                    if B_indices(i-1) <= v
                        if U_indices_new(B_indices(i-1)) == 0
                        rowname1{iter}{i} = strcat('x',num2str(U_id(B_indices(i-1))));
                    elseif U_indices_new(B_indices(i-1)) == 1
                        rowname1{iter}{i} = strcat('x',num2str(U_id(B_indices(i-1))),'+'); 
                    elseif U_indices_new(B_indices(i-1)) == 2
                        rowname1{iter}{i} = strcat('x',num2str(U_id(B_indices(i-1))),'-');     
                        end                      
                    else
                        rowname1{iter_total}{i} = strcat('sx',num2str(B_indices(i-1)));
                    end
                end
                
                dat1{iter_total}{1,1} = '1';
                for i=2:m+1
                    dat1{iter_total}{i,1} = '0';
                end
                
                for i=2:m+n+1
                    dat1{iter_total}{1,i} = num2str(c2(i-1));
                end
                
                dat1{iter_total}{1,m+n+2} = num2str(f);
                
                for i=1:m+n               %% Printing A2
                    for j=1:m
                        dat1{iter_total}{j+1,i+1} = num2str(A2(j,i)) ;
                    end
                end
                
                for i=1:m                 %% Printing bstar
                    dat1{iter_total}{i+1,m+n+2} = num2str(bstar(i));
                end
                
                
                %****************************** Done storing values in table*************************************************
                return
                end  
            end
        end
        return
 end

%% SUBROUTINE function
    function [A1 B bstar ent_var_index B_indices V_indices m u k_art art_indices phase_new rowname1 var_art art_check] = subroutine(A1,B,bstar,ent_var_index,B_indices,V_indices,m,u,k_art,art_indices,n,surp,phase_new,rowname1,iter_total,var_art,art_check,var_art_global)
        piv_col = B\A1(:,ent_var_index);
     
        ratio = zeros(m,1);
        for i=1:m
            if piv_col(i)>0                              %% Denominator positive
                if bstar(i)>=0                           %% Numerator non-negative
                    ratio(i) = bstar(i)/piv_col(i);      %% Finding minimum ratio
                else
                    ratio(i) = inf;
                end
            else
                ratio(i) = inf;
            end
        end
        
        k = 0;
        for i=1:m
            if (ratio(i) == inf)       %% If cannot find ratio
                k = k+1;
            end
        end
        
        if k == m
            if phase_new == 1
                waitfor(msgbox('The given problem has INFEASIBLE SOLUTION !'));
                rowname1{iter}{m+2} = 'INFEASIBLE SOLUTION';
                u = -999;
                return
            else
                waitfor(msgbox('The given problem has UNBOUNDED SOLUTION !'));
                rowname1{iter_total}{m+2} = 'UNBOUNDED SOLUTION';
                u = -999;
                return
            end    
        end
        %   If ratio exists, find the minimum ratio
        min_ratio = min(ratio);
        if numel(min_ratio) >= 2
            for i=1:numel(min_ratio)
                for j=1:k_art
                    if min_ratio(i) == art_indices(j)
                        ratio(art_indices(j)-n-surp) = -Inf ;       %% Make ratio -Inf (least) if artificial variable present
                    end
                end
            end
        end
        min_ratio = min(ratio);
        min_ratio_index = find(ratio == min(ratio));        %% Leaving Variable Index
        
        leave_var_index = B_indices(min_ratio_index(1));
        
        if art_check > 0                 %% If Artificial basic Variable after Phase I 
         if bstar(var_art) == 0           %% If ZERO Artificial Variable after Phase I
            leave_var_index = B_indices(var_art);
         end    
        end
        temp = find(B_indices == leave_var_index);
        B_indices(temp) = ent_var_index;
        temp1 = find(V_indices == ent_var_index);
        V_indices(temp1) = leave_var_index;
        return
    end

%% INTIALISING SECOND PHASE ***********************************************************************************************
% SECOND PHASE FUNCTION
    function second_phase(A1,B,bstar,ent_var_index,B_indices,V_indices,m, u,k_art,art_indices,n,surp,phase_new,f,c,c1,minimize,b,Aj,art,A2,c2,rowname1,total,slack,iter,iter_total,daddy_dat,daddy_row,dat1,var_art,art_check,var_art_global,art_inf,U_indices,U_indices_new,U_id)
     
         iter_total = iter+1;
         
         for i=1:total+2                 %% Printing BLANK ROW
             dat1{iter_total}{1,i} = '' ;
         end
                
         rowname1{iter_total}{1} = strcat('******** PHASE II , Iteration No. ',num2str(iter+1),'  *********');
                
    % Blocking the Artificial Variables for phase 2
        phase_new = phase_new + 1 ;
        c1 = c;
        
        for i=1:m
            for j=1:art
                if B_indices(i) == art_indices(j)
                   art_check = art_check + 1
                    var_art_global = B_indices(i)
                    var_art = i
                end    
            end
        end    
        
       
       for i=1:k_art
                c1(art_indices(i)) = -Inf ;
                c(art_indices(i)) = -Inf ;

       end    
        
       if art_check > 0                 %% If Artificial basic Variable after Phase I 
           if bstar(var_art) == 0           %% If ZERO Artificial Variable after Phase I
               c1(var_art_global) = 0 ;
               c(var_art_global)  = 0 ;
                art_inf = 555;
           end  
       else
               for x=1:k_art           
                    A2(:,art_indices(x)) = 0;
               end             
       end
       if minimize==1               %% For minization problem
            c1 = -c1 ;
       end
    % Done blocking
    
    % Making z co-eff of basic variables zero  
    A2_temp = A2;
    for i=1:m
       var = B_indices(i);
       A2_temp(i,:) = -c2(var)*A2_temp(i,:);    
    end
    
    for i=1:m
       c2 = c2 + A2_temp(i,:);     
    end
    
    % Done making z co-eff of basic variables zero
    
    A1 = A2;
    c1 = c;
    b = bstar;

    %% Begin While Loop
        while (1==1)       
            iter_total = iter_total + 1;
            B = A1(:,B_indices);
            Aj = A1(:,V_indices);
           
           
            if u == -999
                iter_total = iter_total - 1 ;
                final_print_ph2(f,minimize,A2,c2,bstar,m,n,v,cons,B_indices,iter,rowname1,slack,surp,art,total,art_indices,iter_total,daddy_dat,daddy_row,dat1,U_indices,U_indices_new,U_id)
                return
            end
    
            if art_inf == 1000
                c1(var_art_global) = -Inf; 
            end

            if art_inf == 555
                art_inf = 1000;
            end    

            Cb = c1(:,B_indices);
            Cj = c1(:,V_indices);
            f = Cb*(B\b);
            bstar = B\b;            %% B\b = inv(B)*b
            
              if minimize == 1
                    f = -f;
              end
                
            c2 = c1;
            
            cols = total;
            
            for i=1:cols
                size_B = size(find(B_indices == i));
                if size_B(2) ~= 0              %% IF BASIC VARIABLE
                    index_B = find(B_indices ==  i);
                    A2(:,i) = 0;
                    A2(index_B,i) = 1;
                    c2(1,i) = 0;
                else                           %% IF NON-BASIC VARIABLE
                    A2(:,i) = B\A1(:,i) ;
                    index_V = find(V_indices == i);
                    
                    c2(i) = Cb*(B\Aj(:,index_V)) - Cj(index_V) ;
                end
            end
            
             %****************************** Storing values in table *******************************************
            if minimize == 1
                rowname1{iter_total}{1} = 'z(min)';
            else
                rowname1{iter_total}{1} = 'z(max)';
            end
            
            for i=2:m+1
                if B_indices(i-1) <= v      %% Decision Variable
                    if U_indices_new(B_indices(i-1)) == 0
                        rowname1{iter_total}{i} = strcat('x',num2str(U_id(B_indices(i-1))));
                    elseif U_indices_new(B_indices(i-1)) == 1
                        rowname1{iter_total}{i} = strcat('x',num2str(U_id(B_indices(i-1))),'+'); 
                    elseif U_indices_new(B_indices(i-1)) == 2
                        rowname1{iter_total}{i} = strcat('x',num2str(U_id(B_indices(i-1))),'-');     
                    end   
                else
                    if (B_indices(i-1) >= v+1) && (B_indices(i-1) <= v+surp)        %% Surplus Variable
                        rowname1{iter_total}{i} = strcat('Sx',num2str(B_indices(i-1)));
                    else    
                        var = size(find(art_indices == i-1+n+surp));
                        if var(2) == 0
                            rowname1{iter_total}{i} = strcat('sx',num2str(B_indices(i-1)));
                        else
                            rowname1{iter_total}{i} = strcat('Rx',num2str(B_indices(i-1)));
                        end    
                            
                    end    
                end
            end
            
            
            dat1{iter_total} = [];
            dat1{iter_total}{1,1} = '1';
            
            for i=2:m+1
                dat1{iter_total}{i,1} = '0';
            end
            
            for i=2:total+1
                if minimize == 1
                    dat1{iter_total}{1,i} = num2str(-c2(i-1));
                else
                    dat1{iter_total}{1,i} = num2str(c2(i-1));
                end    
            end
            
            dat1{iter_total}{1,total+2} = num2str(f);
            
            for i=1:total               %% Printing A2
                for j=1:m
                    dat1{iter_total}{j+1,i+1} = num2str(A2(j,i)) ;
                end
            end
            
            for i=1:m                 %% Printing bstar
                dat1{iter_total}{i+1,total+2} = num2str(bstar(i));
            end          
        
            %****************************** Done storing values in table *************************************************
            
            ent_var_mat = Cb*(B\Aj) - Cj;
            
            if min(ent_var_mat) >= 0                                      %% Optimality condition satisfied
                             
              [B_indices B bstar iter A2 c2 f rowname1 iter_total dat1 U_indices U_indices_new U_id u] = alternate_check(A1,B, bstar, B_indices, V_indices, m,ent_var_mat,u,n,A2,minimize,iter,c2,f,b,c1,ent_var_index,art_indices,surp,rowname1,iter_total,total,dat1,U_indices,U_indices_new,U_id);           %% Check for ALTERNATE SOLUTIONS
                
                if u == -999
                    iter_total = iter_total - 1 ;           
%                     return
                end
              
                final_print_ph2(f,minimize,A2,c2,bstar,m,n,v,cons,B_indices,iter,rowname1,slack,surp,art,total,art_indices,iter_total,daddy_dat,daddy_row,dat1,U_indices,U_indices_new,U_id)        %% Call final_print function
                return
            else
                
                for i=1:total+2                 %% Printing BLANK ROW
                    dat1{iter_total}{m+2,i} = '' ;
                end
                
                rowname1{iter_total}{m+2} = strcat('******** Iteration No. ',num2str(iter_total),'  *********');
                
                ent_var_index = find(ent_var_mat == min(ent_var_mat));
                ent_var_index = V_indices(ent_var_index(1));                %% Entering Variable
                
                art_inf = art_inf + 1;
                [A1 B bstar ent_var_index B_indices V_indices m u k_art art_indices phase_new rowname1 var_art art_check] = subroutine(A1,B,bstar,ent_var_index,B_indices,V_indices,m, u,k_art,art_indices,n,surp,phase_new,rowname1,iter_total,var_art,art_check);    %% Calling the SUBROUTINE function
            end
            
          
        end         %% End of while loop
    
       
    end

%% SUBROUTINE for ALTERNATE function
    function [A1 B bstar ent_var_index B_indices V_indices m u k_art art_indices] = subroutine_alternate(A1,B,bstar,ent_var_index,B_indices,V_indices,m,u,k_art,art_indices,n,surp)
        piv_col = B\A1(:,ent_var_index);
        
        ratio = zeros(m,1);
        for i=1:m
            if piv_col(i)>0                              %% Denominator positive
                if bstar(i)>=0                           %% Numerator non-negative
                    ratio(i) = bstar(i)/piv_col(i);      %% Finding minimum ratio
                else
                    ratio(i) = inf;
                end
            else
                ratio(i) = inf;
            end
        end
        k = 0;
        for i=1:m
            if (ratio(i) == inf)       %% If cannot find ratio
                k = k+1;
            end
        end
        
        if k == m         
                waitfor(msgbox('The given problem has UNBOUNDED SOLUTION !'));
                u = -999;
                return         
        end
        
        %   If ratio exists, find the minimum ratio
        min_ratio = min(ratio);
        if numel(min_ratio) >= 2
            for i=1:numel(min_ratio)
                for j=1:k_art
                    if min_ratio(i) == art_indices(j)
                        ratio(art_indices(j)-n-surp) = -Inf ;       %% Make ratio -Inf (least) if artificial variable present
                    end
                end
            end
        end
        
        min_ratio = min(ratio);
        min_ratio_index = find(ratio == min(ratio));        %% Leaving Variable Index
        
        leave_var_index = B_indices(min_ratio_index(1));
        
        temp = find(B_indices == leave_var_index);
        B_indices(temp) = ent_var_index;
        temp1 = find(V_indices == ent_var_index);
        V_indices(temp1) = leave_var_index;
        return
    end

%% Final PRINT Phase 2 function
function final_print_ph2(f,minimize,A2,c2,bstar,m,n,v,cons,B_indices,iter,rowname1,slack,surp,art,total,art_indices,iter_total,daddy_dat,daddy_row,dat1,U_indices,U_indices_new,U_id)

          % ************** PRINTING OPTIMAL TABLE BEGINS HERE **********************************************************************
        T = figure('Name','Two Phase Method Final Table Values','NumberTitle','off','Position',[50 50 1200 600],...
            'DeleteFcn',@Figure_Close_Callback);
        columnname = {};
        columnformat = {};
        
        columnname{1} = 'z';    %% First column name is 'z'
      
        
         j=2;
        for i=2:v+1               %% Columnames for decision variables
            if U_indices(i-1) == 0
                str1 = strcat('x',num2str(i-1));
                columnname{j} = str1 ;
                j=j+1;
            else
                str1 = strcat('x',num2str(i-1),'+');
                str2 = strcat('x',num2str(i-1),'-');
                columnname{j} = str1 ;
                columnname{j+1} = str2 ;
                j = j+2;
            end    
        end
        
        for i=v+2:v+2+surp-1           %% Columnnames for surplus variables
            str1 = strcat('Sx',num2str(i-1));
            columnname{i} = str1 ;
        end
        
        for i=(v+surp+2):total+1
            var = size(find(art_indices == i-1));
            if var(2) == 0
                str1 = strcat('sx',num2str(i-1));
                columnname{i} = str1 ;
            else
                str1 = strcat('Rx',num2str(i-1));
                columnname{i} = str1 ;
            end    
        end
               
        columnname{total+2} = 'Solution';
        
        iter_primal = iter_total;
                   
        if mod(iter_primal,2) == 0     %% If iter_primal is even 
            if mod(iter,2) == 0        %% If iter is even
                for i=iter+1:2:iter_primal
                    xyz = vertcat(dat1{i},dat1{i+1});
                    daddy_dat = vertcat(daddy_dat,xyz);

                    lol = horzcat(rowname1{i},rowname1{i+1});
                    daddy_row = horzcat(daddy_row,lol);
                end
            else                       %% If iter is odd
                for i=iter+1:2:(iter_primal - 1)
                    xyz = vertcat(dat1{i},dat1{i+1});
                    daddy_dat = vertcat(daddy_dat,xyz);

                    lol = horzcat(rowname1{i},rowname1{i+1});
                    daddy_row = horzcat(daddy_row,lol);
                    
                end
                daddy_dat = vertcat(daddy_dat,dat1{iter_primal});
                daddy_row = horzcat(daddy_row,rowname1{iter_primal});
            end    
        else                                    %% If iter_primal is Odd
            if mod(iter,2) == 0        %% If iter is even
                for i=iter+1:2:(iter_primal - 1)
                    xyz = vertcat(dat1{i},dat1{i+1});
                    daddy_dat = vertcat(daddy_dat,xyz);

                    lol = horzcat(rowname1{i},rowname1{i+1});
                    daddy_row = horzcat(daddy_row,lol);
                end
               
                daddy_dat = vertcat(daddy_dat,dat1{iter_primal});
                daddy_row = horzcat(daddy_row,rowname1{iter_primal});
            else                        %% If iter is odd
                for i=iter+1:2:(iter_primal - 1)
                    xyz = vertcat(dat1{i},dat1{i+1});
                    daddy_dat = vertcat(daddy_dat,xyz);

                    lol = horzcat(rowname1{i},rowname1{i+1});
                    daddy_row = horzcat(daddy_row,lol);
                end
            end    
        end
        
        q = uitable('Parent',T,...
            'Position',[0 0 1200 600],'Data', daddy_dat,...
            'ColumnName', columnname,...
            'ColumnEditable', false,...
            'FontSize',10,'ForegroundColor','k', ...
            'FontName','Comic Sans MS', ...
            'RowName',daddy_row );
        
        
        %*************** PRINTING OPTIMAL TABLE ENDS HERE **********************************************************************
            
end
%% Check Function For Bounds
 function [temp1 temp2] = check(temp1,temp2,g,er)
                if str2num(temp1) > str2num(temp2)
                        str = strcat('x',num2str(er));
                        waitfor(msgbox('The Lower bound is greater than Upper Bound ! Re-Enter the Values.',str,'warn'));
                        
                        prompt = {'Enter the Lower Bound : ','Enter the Upper Bound :'};
                        dlg_title = str;
                        num_lines = 1;
                        def = {'0','Inf'};
                        answer = inputdlg(prompt,dlg_title,num_lines,def);

                        temp1 = answer{1};
                        temp2 = answer{2};
                        
                        [temp1 temp2] = check(temp1,temp2,g,er);
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
%                     close all; clear all;
                    return  
                    
                case 'View/Modify Input Data'   
                    twop_modify(g_modify,v_org,cons_org);
            end       
 end

end
