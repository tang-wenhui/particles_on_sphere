%% Read .xlsx files
    %delete('defect fraction.xlsx');
    xlsxfiles = dir('*.xlsx');
    xlsxfiles = {xlsxfiles.name}';
    
    Npts_ = [];        % cell number
    Hexagon_ = [];
    Largerthansix_=[];
    Lessthansix_=[];
    Hexagon=[];
    Lessthansix=[];
    Largerthansix=[];
    
    t = zeros(length(xlsxfiles),1);
    for i = 1:length(xlsxfiles)
        % initialize variables 
        
        file_dir = xlsxfiles(i);
        %% import data
        data = readtable(file_dir{1});
        t(i) = string(strtok(file_dir{1},'.'));
        % basic parameters
        Hexagon{i} = data(:,4);
        Hexagon{i} = cell2mat(table2cell(Hexagon{i}));
        Hexagon_ave(i) = mean(Hexagon{i});
        Hexagon_std(i) = std(Hexagon{i});
   
        Npts = data(:,9);
        Npts = cell2mat(table2cell(Npts));
        Npts_ave(i) = mean(Npts);
        
        Lessthansix{i} = data(:,10);
        Lessthansix{i} = cell2mat(table2cell(Lessthansix{i}));
        Lessthansix_ave(i) = mean(Lessthansix{i});
        Lessthansix_std(i) = std(Lessthansix{i});
        
        Largerthansix{i} = data(:,11);
        Largerthansix{i} = cell2mat(table2cell(Largerthansix{i}));
        Largerthansix_ave(i) = mean(Largerthansix{i});
        Largerthansix_std(i) = std(Largerthansix{i});
        
        
    end
    
    
%     polygon = zeros(length(Hexagon_),4);
%     polygon(:,1) = Hexagon_;
%     polygon(:,2) = Lessthansix_
%     polygon(:,3) = Largerthansix_
%     polygon(:,4) = Npts_;
%  
%     T1 = table(polygon);
%     writetable(T1,'defects fraction_various interaction.xlsx');
    
    
    % average vertices number
    %V_ave = 3*Triangle_+4*Rectangle_+5*Pentagon_+6*Hexagon_+7*Heptagon_+8*Octagon_+9*Nonagon_+10*Decagon_;
    % plot ploygon fraction as a function of cell number N
%     figure(1);
%     s1 = scatter(Npts_, polygon(:,1),14,'o','MarkerFaceColor',[0.9290, 0.6940, 0.1250],'MarkerEdgeColor',[0.9290, 0.6940, 0.1250]); hold on
%     s2 = scatter(Npts_, polygon(:,2),14,'o','MarkerFaceColor',[0.4660, 0.6740, 0.1880],'MarkerEdgeColor',[0.4660, 0.6740, 0.1880]); hold on
%     s3 = scatter(Npts_, polygon(:,3),14,'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on
%     s1.MarkerFaceAlpha = 0.4;
%     s2.MarkerFaceAlpha = 0.4;
%     s3.MarkerFaceAlpha = 0.4;
%     s1.MarkerEdgeAlpha = 0.4;
%     s2.MarkerEdgeAlpha = 0.4;
%     s3.MarkerEdgeAlpha = 0.4;
%     xlabel('cell number N','FontSize',22); ylabel('polygon fraction','FontSize',22);
%     set(gca,'FontSize',18);
%     legend('z<6','z=6','z>6');
    
    
%     figure(2);
%     v = [];
%     for i = 4:1:1000
%         v(i-3) = 6-12/i;
%     end
%     s6 = scatter(4:1:1000,v,14,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);hold on
%     s7 = scatter(Npts_,V_ave,14,'o','MarkerFaceColor',[0.6350, 0.0780, 0.1840],'MarkerEdgeColor',[0.6350, 0.0780, 0.1840]);
%     s6.MarkerFaceAlpha = 0.1;
%     s6.MarkerEdgeAlpha = 0.1;
%     xlabel('cell number N','FontSize',25); ylabel('Average vertex number','FontSize',25);
%     set(gca,'FontSize',17);
%     legend('ideal','Alveolospheres');

%% boxplot
        figure(1);
        s1 = boxplot([Hexagon{1},Hexagon{2},Hexagon{3},Hexagon{4},Hexagon{5},Hexagon{6},Hexagon{7},Hexagon{8},Hexagon{9}]);
        %s1 = boxplot([Hexagon{1},Hexagon{2},Hexagon{3}]);
        ylabel('polygon fraction');title('N = 800 Z=6');
        ylim([0.2 0.5]);
        
        figure(2);
        s2 = boxplot([Lessthansix{1},Lessthansix{2},Lessthansix{3},Lessthansix{4},Lessthansix{5},Lessthansix{6},Lessthansix{7},Lessthansix{8},Lessthansix{9}]); 
        %s2 = boxplot([Lessthansix{1},Lessthansix{2},Lessthansix{3}]);
        ylabel('polygon fraction');title('N = 800 Z<6');
        ylim([0.2 0.5]); 
        
        figure(3);
        s3 = boxplot([Largerthansix{1},Largerthansix{2},Largerthansix{3},Largerthansix{4},Largerthansix{5},Largerthansix{6},Largerthansix{7},Largerthansix{8},Largerthansix{9}]);
        %s3 = boxplot([Largerthansix{1},Largerthansix{2},Largerthansix{3}]);
        ylabel('polygon fraction');title('N = 800 Z>6');
        ylim([0.2 0.5]);
        