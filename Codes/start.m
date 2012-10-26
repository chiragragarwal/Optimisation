function start()
clear all
close all
clc

choice = menu('SELECT AN ALGORITHM', ...
              'All Slack Primal Simplex', ...
              'Big M Method Simplex',... 
              'Two Phase Method Simplex',...
              'Dual Simplex','Generalized Simplex',...
              'Floyd Warshall Algorithm',...
              'Kruskal Algorithm',...
              'Complete Graph',...
              'Graph from Adjacency Matrix',...
              'Fractional Knapsack Algorithm',...
              'Task Scheduling Algorithm');

switch choice
    case 1
        primal();
        return
    case 2
        bigm();
        return
    case 3
        twop();
        return
    case 4
        dual();
        return
    case 5
        generalized();
        return
    case 6
        floyd_warshall();
        return    
    case 7
        kruskal();
        return
    case 8
        graph_complete();
        return
    case 9
        graph_adj();
        return    
    case 10
        knap();
        return
    case 11
        task();
        return    
end    
end

