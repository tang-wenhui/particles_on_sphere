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

%for n = [800]
for n = [1200]
    Npts = n; % number of cells on the sphere
    % position of Voronoi centers in spherical coords
    % restrictions: the distance of two particles cannot be closer than a
    % threshold (minDist = sqrt(4*pi*R^2/Npts*2/1.732)*x)) (x% of the average distance)
    % minDist = sqrt(1/Npts)*8*0.4;  % minimum distance between two particles
    mm=0;
    for xmin = [0]
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
    
        %% compute pair correlation function (average over all cells)
        
        g_r_ave = [];
        angle = nan(Npts,Npts);
        %% method1: calculate radial distribution function g_r
%         % calculate the angles between every two particles
%         for j = 1:Npts-1
%             for k = j+1:Npts
%                 angle(j,k) = acos(dot(voro_xyz(j,:),voro_xyz(k,:)));
%             end
%         end
%         clear j k
%         for j = 1:Npts
%             % define the size of delta_theta
%             size_theta = 600;
%             delta_theta = pi/size_theta;
%             center_xyz = voro_xyz(j,:);
%             for ii = 1: size_theta-1 %loop for every r
%                 %count the number of cells within delta_theta
%                 count2(ii) = sum((angle(j,:) < ii*delta_theta & angle(j,:)>(ii-1)*delta_theta)) + sum((angle(:,j) < ii*delta_theta & angle(:,j)>(ii-1)*delta_theta));
%                 g_r(j,1) = 0;   
%                 g_r(j,ii+1) = 2*count2(ii)/(abs(sin(ii*delta_theta))*delta_theta*Npts);
% %                   g_r: local density/bulk density
% %                   local density: count/(2*pi*sin(ii*delta_theta)*delta_theta)
% %                   bulk density: 
%             end
%         end
%         g_r_(i,:) = mean(g_r); % average over N

    %% method2: calculate radial distribution function
    % calculate for all cells not one by one
    cc = voro_xyz;
    angle = acos(cc*cc');
    size_theta = 600;
    delta_theta = pi/size_theta;
    for ii = 1:size_theta-1
        [index1,index2] = find(angle(:,:) < ii*delta_theta & angle(:,:)>(ii-1)*delta_theta);
        count = length(index1);
        g_r(1) = 0;
        g_r(1+ii)= 2*count/(abs(sin(ii*delta_theta))*delta_theta*Npts^2); % local density/bulk density and averaged over cell number N
    end
    g_r_(i,:) = g_r;
    
    % plot the corresponding pair interaction potential U(s)/KBT
%     figure(2);
%     U = -log(g_r_ave); %U(s)/kBT
%     plot(s,U,'-','Linewidth',2.5,'Color',color_); hold on
%     xlabel('Geodesic distance r\theta');ylabel('U(r\theta)/k_BT');
%     title('N=',n); 
%     xlim([0 pi/8]); 
%     %ylim([0 3]);
%     box on;
%     set(gca,'fontsize',16);
%     xticks([0 pi/32 pi/16 3*pi/32 pi/8]);
%     xticklabels({'0','pi/32','pi/16','3pi/32','\pi/8'});


     %% plot the structure factor
%     figure(2);
%     Y = fft(g_r_ave);
%     T = x(2)-x(1); % sampling period
%     Fs = 1/T;      % sampling frequency
%     L = size_theta; % length of signal
%     f = Fs*(1:(L/2))/L; %frequency
% 
%     P2 = abs(Y/L);
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
% 
%     plot(f, P1(1:end-1), '-','Linewidth',1.5); hold on
%     xlabel('f (Hz)');ylabel('S(q)');
%     title('N=',n);
%     xlim([0 60]); 
%     box on;
%     set(gca,'fontsize',16);

    %% need to put this outside of for i=1:nframe loop
    %legend('xmin = 0','0.1','0.2','0.3','0.4','0.5','0.6','0.7');
    %print -depsc -tiff -r300 -painters UkT_Xmin.eps


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

%         if mod(length(X),ncl) >= 2
%             cl = clmap(mod(length(X),ncl)-2,:);    % length(X) is the number of cells around one cell
%         else
%             cl = clmap(9,:);
%         end
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

%     function theta_phi_r = cartesian_to_sphere(xyz)
%         r = sqrt(sum(xyz.^2,2));
%         theta = mod(atan2(xyz(:,2),xyz(:,1)),2*pi);
%         phi = acos(xyz(:,3)./r);
%         theta_phi_r = [theta, phi,r];
% 
%     end

    if WRITE_
        close(writerObj);
    end
    end
    % average radial distribution function over frames
    g_r_ave = mean(g_r_,1); % average over frames
    % plot radial distribution function
    color_ = [0.9290-0.9290/double(mm) 0.6940-0.247/double(mm) 0.1250+0.616/double(mm)]; 
    figure(2);
    s =(0:size_theta-1)*delta_theta;
    g_r_ave = smoothdata(g_r_ave, 'Gaussian',2);
    plot(s,g_r_ave,'-','Linewidth',2.5,'Color',color_); hold on
    xlabel('Geodesic distance r');ylabel('g(r\theta)');
    title('N=',n);
    xlim([0 pi/8]); 
    %ylim([0 3]);
    box on;
    set(gca,'fontsize',16);
    xticks([0 pi/32 pi/16 3*pi/32 pi/8]);
    xticklabels({'0','pi/32','pi/16','3pi/32','\pi/8'});
    %print -depsc -tiff -r300 -painters radialdistributionfunction.eps
    
end
% save polygon data in polygon
 %remove the zeros terms(blank frames) in Npts
%     Triangle(find(Triangle==0)) = [];
%     Rectangle(find(Rectangle==0)) = [];
%     Pentagon(find(Pentagon==0)) = [];
%     Hexagon(find(Hexagon==0)) = [];
%     Heptagon(find(Heptagon==0)) = [];
%     Octagon(find(Octagon==0)) = [];
%     Nonagon(find(Nonagon==0)) = [];
%     Decagon(find(Decagon==0)) = [];
    %Npts(Npts==0) = [];
    
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
  
   T2 = table(polygon);
   writetable(T2,strcat(string(Npts),['_polygon fraction_xmin=0.45.xlsx']));


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