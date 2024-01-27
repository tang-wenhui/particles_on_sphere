%% Read .xlsx files
clc;clear all;clear;

xlsxfiles = dir('*.xlsx');
xlsxfiles = {xlsxfiles.name}';
t = zeros(length(xlsxfiles),1);
Largerthansix = [];
Lessthansix = [];
Hexagon_std = [];
Largerthansix_std = [];
Lessthansix_std = [];
Npts_=[];
Hexagon_=[];

for i = 1:length(xlsxfiles)
        % initialize variables 
        file_dir = xlsxfiles(i);
        %% import data
        data = readtable(file_dir{1});
        t(i) = string(strtok(file_dir{1},'.'));
        % basic parameters
        Npts = data(:,9); 
        Npts = cell2mat(table2cell(Npts));
        Npts_ = [Npts_;mean(Npts)];
        
        Hexagon = data(:,4);
        Hexagon = cell2mat(table2cell(Hexagon));
        Hexagon_ = [Hexagon_;mean(Hexagon)];
        Hexagon_std = [Hexagon_std;std(Hexagon)];
        
        Largerthan_six = data(:,11);
        Largerthan_six = cell2mat(table2cell(Largerthan_six));
        Largerthansix = [Largerthansix;mean(Largerthan_six)];
        Largerthansix_std = [Largerthansix_std;std(Largerthan_six)];
        
        Lessthan_six = data(:,10);
        Lessthan_six = cell2mat(table2cell(Lessthan_six));
        Lessthansix = [Lessthansix;mean(Lessthan_six)];
        Lessthansix_std = [Lessthansix_std;std(Lessthan_six)];
end
 
 figure(1);
 s1=errorbar(Npts_, Hexagon_, Hexagon_std/2, Hexagon_std/2,'o','MarkerSize',4,'linewidth',1.5); hold on
 s1.Color = [0 0 0]; 
 s2=errorbar(Npts_, Lessthansix, Lessthansix_std/2, Lessthansix_std/2, 'o','MarkerSize',4,'linewidth',1.5); hold on
 s2.Color = [1 0 0];
 s3=errorbar(Npts_, Largerthansix, Largerthansix_std/2, Largerthansix_std/2, 'o','MarkerSize',4,'linewidth',1.5); 
 s3.Color = [0 0 1];
 xlabel('Cell number N'); ylabel('Polygon fraction');
 xlim([0 1000]);ylim([0 1]);
 set(gca,'FontSize',18);
 %set(gca, 'YScale', 'log','XScale', 'log')
 legend('z=6 Simulation','z<6 Simulation','z>6 Simulation');
 print -depsc -tiff -r300 -painters polygonfraction_reproducealveolosphereusingsimulation.eps
    
 