% this is a sample script to create the Voronoi tessellation on a sphere
% for the purpose of demonstration, we will generate random points on a
% sphere to use as the centers of the Voronoi tessellation 
clc; clear all;clear;
% If analyzing experimental data, the input should be the nuclei positions on a unit sphere

% ---------------Initialization----------------------------- %
nframe = 10; % average over 30 frames
vertices = zeros(nframe,3);                                                % creat vertices array
WRITE_ = 0;                                                                % whether write video
particleArea = [];
AREA_ = 1;


for n = [1000]
    Npts = n; % number of cells on the sphere
    % position of Voronoi centers in spherical coords
    % restrictions: the distance of two particles cannot be closer than a
    % threshold (minDist = sqrt(4*pi*R^2/Npts*2/1.732)*x)) (x% of the average distance)
    % minDist = sqrt(1/Npts)*8*0.4;  % minimum distance between two particles
    mm=0;
    for xmin = [0.1]
    mm = mm+1;
    ave_cell_size = sqrt(4*pi*1^2/Npts*2/1.732);
    minDist = ave_cell_size*xmin;
    if WRITE_
        writerObj = VideoWriter(strcat(num2str(Npts),'_voronoi tessellation'));
        writerObj.FrameRate = 1;
        open(writerObj);
    end
    
    for i = 1:nframe
        c = 0; %Counter
        voro_theta_phi = nan(Npts,2);
        % The main idea is to start off with a set of random coordinates and continually 
        % replace coordinates that are too close to another coordinate until either 1) 
        % all distances are less than the minimum distance requirement or 2) 10 attempts 
        % were made per coordinate.
        while c<(100*Npts) && any(isnan(voro_theta_phi(:)))
        % position of Voronoi centers in spherical coords
            nanlist_num = find(isnan(voro_theta_phi(:,1)));
            for k = 1:length(nanlist_num)
                voro_theta_phi(nanlist_num(k),:) = [rand(1,1)*2*pi,acos(2*rand(1,1)-1)]; 
            end
            % converting to cartesian coordinates (assuming the sphere has radius of 1)
            voro_xyz = [cos(voro_theta_phi(:,1)).*sin(voro_theta_phi(:,2)),sin(voro_theta_phi(:,1)).*sin(voro_theta_phi(:,2)),cos(voro_theta_phi(:,2))];
            % Identify rows that are too close to another point
            [~,isTooClose] = find(triu(squareform(pdist(voro_xyz)) < minDist,1));
            % Replace too-close coordinates with NaN
            voro_theta_phi(isTooClose,:) = NaN; 
            c = c+1;
            clear nanlist_num isTooClose
        end
        
        % Throw error if the loop had to quit and missing values remain
        if any(isnan(voro_theta_phi(:)))
            error('The while-loop gave up. There are %d coordinates with missing values.',sum(isnan(voro_theta_phi(:,1))))
        end
        % Display number of attempts 
        fprintf('%d number of attempts.\n', c)
        % Show the minimum distance 
        distances{i} = triu(squareform(pdist(voro_theta_phi))); 
        fprintf('Min distance = %.2f\n', min(distances{i}(distances{i}~=0)));
        % show particle distance matrix after loop, mean and std
        ave_size = sqrt(4*pi*1^2/Npts*2/1.732);
    
   
    


    %% ---------------Compuate Voronoi-------------------------- %
    [Vertices_xyz, Cell_List, ~] = voronoisphere(voro_xyz');
    %Output content : Vertices_xyz = the list of vertex positions of the Voronoi polygons
    %                 Cell_List =  contains the indices of the voronoi cell
    %                              vertices correspond to each cell
    sa = vcell_solidangle(Vertices_xyz, Cell_List); % solid angle of voronoi polygons
    area{i} = sa;
    particleArea = [particleArea;area{i}];

    % Neighbor count for each cell
    Z_List = cellfun('length',Cell_List);
    % calculate polygon percentage in each frame
    Triangle(i) = length(find(Z_List == 3))/Npts;
    Rectangle(i) = length(find(Z_List == 4))/Npts;
    Pentagon(i) = length(find(Z_List == 5))/Npts;
    Hexagon(i) = length(find(Z_List == 6))/Npts;
    Heptagon(i) = length(find(Z_List == 7))/Npts;
    Octagon(i) = length(find(Z_List == 8))/Npts;
    Nonagon(i) = length(find(Z_List == 9))/Npts;
    Decagon(i) = length(find(Z_List == 10))/Npts;
    Hendecagon(i) = length(find(Z_List == 11))/Npts;
    Dodecagon(i) = length(find(Z_List == 12))/Npts;
    defect = 6-Z_List;
    defect = sum(defect);
    vertices(i,1) = mean(Z_List); 
    vertices(i,2) = std(Z_List);
    vertices(i,3) = Npts;
    % Video
    fprintf(['-> Video playing...','\n']);
    f = figure(1);
    clf(f);
    set(f,'Renderer','zbuffer');
    ax = axes('Parent', f);
    hold(ax, 'on');
    axis(ax,'equal');

    % plot3(ax, voro_xyz(1,:),voro_xyz(2,:),voro_xyz(3,:),'wo');
    mymap1 = [0.9 0.0780 0.1840
    0.9 0.0780 0.1840
    0.9 0.0780 0.1840
    1 1 1
    0 0 1
    0 0 1
    0 0 1
    0 0 1
    0 0 1
    0 0 1
    0 0 1
    0 0 1];
    mymap2 = parula(9);
    mymap3 = [
            1 0 0.1294
            1 0.6 0.672
            1 0.85 0.952
            1 1 1
            0.9 0.9 1
            0.8 0.8 1
            0.5 0.5 1
            0.3 0.3 1
            0 0 1];
    clmap = colormap(mymap3);
    ncl = size(clmap,1);
    zlist = zeros(Npts,1);
    cl = [ 0.5843    0.8157    0.9882];
    for k = 1:Npts
        zlist(k) = size(Cell_List{k},1);
        X = Vertices_xyz(:,Cell_List{k});
        if length(X) >= 12
            cl = clmap(9,:);
        else
            cl = clmap(mod(length(X),12)-2,:);
        end

        fill3(X(1,:),X(2,:),X(3,:),[0.5 0.5 0.5],'Parent',ax,'EdgeColor','k','linewidth',0.01,'FaceColor',cl); 
        set(0,'defaultfigurecolor','w');
    end
    colorbar;caxis([3,12]);
    title(strcat('defect =  ',num2str(defect),'  cell number =  ', num2str(Npts)),'FontSize',18);
    view([1,1,1])
    axis(ax,'equal');
    set(gca,'linewidth',2)
    set(0,'defaultfigurecolor','w');
    %scatter3(voro_xyz(:,1),voro_xyz(:,2),voro_xyz(:,3),200,'k','filled')

    axis(ax,[-1 1 -1 1 -1 1]);
    grid on
    set(gca,'xtick',-1:0.5:1)
    set(gca,'ytick',-1:0.5:1)
    set(gca,'ztick',-1:0.5:1)
    set(gca, 'Color', 'none');
    %colorbar; caxis([3,12]);
    set(gca,'visible','off');
    set(gca,'xtick',[]);
    if WRITE_
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end



    if WRITE_
        close(writerObj);
    end
    end
    
end

    
    polygon = zeros(nframe,11);
    polygon(:,1) = Triangle;
    polygon(:,2) = Rectangle;
    polygon(:,3) = Pentagon;
    polygon(:,4) = Hexagon;
    polygon(:,5) = Heptagon;
    polygon(:,6) = Octagon;
    polygon(:,7) = Nonagon;
    polygon(:,8) = Decagon;
    polygon(:,9) = Npts;
    polygon(:,10) = Triangle + Rectangle + Pentagon;
    polygon(:,11) = ones(nframe,1)-polygon(:,10)-Hexagon';
    
  
%    T2 = table(polygon);
%    writetable(T2,strcat(string(Npts),'_polygon fraction_reproduce.xlsx'));


end

%% histogram of area
%particle-particle distance/average_cell_size
%     figure;
%     histogram(sqrt(particleArea)/ave_cell_size,'BinWidth',0.02,'Normalization','pdf','FaceColor','k','EdgeColor','k','FaceAlpha',0.3); hold on
%     %set(gca,'XScale','log');
%     xlabel('\fontsize{25} particle-particle distance/average cell size'); ylim([0 6]); xlim([0 2]);
%     ylabel('\fontsize{25} PDF');
%     title(strcat(' Average {\itN} =  ', '  ', num2str(mean(Npts))));
%     box on;
%     set(gca,'linewidth',1.5,'fontsize',16);
%     print -depsc -tiff -r300 -painters P2P_N=800_x=0.6.eps

% if AREA_  % solid angle
%     figure;
%     histogram(particleArea,'BinWidth',0.001,'Normalization','pdf','FaceColor','k','EdgeColor','k','FaceAlpha',0.3); hold on
%     %set(gca,'XScale','log');
%     xlabel('\fontsize{25} solid angle'); xlim([0 0.04]);ylim([0 250]);
%     ylabel('\fontsize{25} PDF');
%     title(strcat(' Average {\itN} =  ', '  ', num2str(mean(Npts))));
%     box on;
%     set(gca,'linewidth',1.5,'fontsize',16);
%     print -depsc -tiff -r300 -painters solidangle_N=1600_simulation_x=0.2.eps
% end 