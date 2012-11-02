function [zones,cluster_paths] =  ExtractROIs_Kmeans(images_aligned,BW,m,detrend,Full_analysis)

global scaling

numclusters = 2;
mask = BW.*(~m);
[pixel_coords(:,2),pixel_coords(:,1)] = find(mask>0);
data = images_aligned.*repmat(reshape(mask,[1,size(BW)]),[19,1,1]);
data_whole_vec = reshape(data,19,[]);
[indices] = find(data_whole_vec(3,:)>0); % Only consider pixels in mask. 3rd frame because it is always populated
data_vec = data_whole_vec(:,indices);
full_frames = find(sum(data_vec,2) ~= 0);
data_vec(sum(data_vec,2) == 0,:) = [];  % Exclude "missing frames"
if detrend == 'T'
    data_vec = data_vec - repmat(mean(data_vec),size(data_vec,1),1);
end

% K-means
[idx,cluster_paths,Jk] = kmeans(data_vec',numclusters,'Distance','sqEuclidean','replicates',10);
% [idx,cluster_paths] = GM_EM(data_vec',numclusters,'T',20,100,0.001);
% Use median for cluster path
cluster_paths(1,:) = median(data_vec(:,idx==1),2);
cluster_paths(2,:) = median(data_vec(:,idx==2),2);

% Ward's group-average agglomerative method
% Z = linkage(data_vec','ward','euclidean');
% idx = cluster(Z,'maxclust',2);
% cluster_paths(1,:) = mean(data_vec(:,idx == 1)');
% cluster_paths(2,:) = mean(data_vec(:,idx == 1)');



zones = repmat(zeros(size(BW)),[1,1,3]);
for j = 1:numclusters
    for i = 1:length(pixel_coords(idx==j,1))
        zone_pixels = pixel_coords(idx==j,:);
        if j < 4
            zones(zone_pixels(i,2), zone_pixels(i,1),j) = 255;
        elseif j == 4
            zones(zone_pixels(i,2), zone_pixels(i,1),1) = 255;
            zones(zone_pixels(i,2), zone_pixels(i,1),2) = 255;
        end
    end
end

if Full_analysis=='N'
    figure('Position',[100,100,800,400])
else
    figure('Position',[100,100,800,800])
end

%Automatic selection
for i = 1:numclusters
    zonemean(i) = mean(mean(data_vec(:,idx == i)));
end
[~,o] = sort(zonemean);
o = fliplr(o);
for i = 1:numclusters
    zone(i).idx = o(i);
end

%Re-shuffle idx
idx2 = idx;
for i = 1:numclusters
    idx2(idx == zone(i).idx) = i;
end
idx = idx2;

cluster_paths2 = cluster_paths;
for i = 1:numclusters
    cluster_paths2(i,:) = cluster_paths(zone(i).idx,:);
end
cluster_paths = cluster_paths2;


%Recolor zones
zones = repmat(zeros(size(BW)),[1,1,3]);
for j = 1:numclusters
    for i = 1:length(pixel_coords(idx==j,1))
        zone_pixels = pixel_coords(idx==j,:);
        zones(zone_pixels(i,2), zone_pixels(i,1),j) = 255;
    end
end

% Put back into 19xN vector
cluster_paths2 = zeros(numclusters,19);
cluster_paths2(:,full_frames') = cluster_paths;
cluster_paths = cluster_paths2;

% Mark zones
hold off

if Full_analysis == 'N'
    subplot(1,2,1)
else
    subplot(3,3,1); 
end
I = squeeze(images_aligned(3,:,:));
I2 = repmat(I,[1,1,3]);
I3 = I2; I3(zones>0) = I3(zones>0)*0.1.*(zones(zones>0));
imshow(I3);

if Full_analysis == 'N'
    subplot(1,2,2)
else
    subplot(3,3,2); % imshow(zones)
end

hold off
plot(scaling./data_vec(:,idx == 1),'r'); hold on
plot(scaling./data_vec(:,idx == 2),'g');
h1=plot(scaling./cluster_paths(1,:),'k--','Linewidth',2);
h2=plot(scaling./cluster_paths(2,:),'k-','Linewidth',2);
if Full_analysis == 'N'
    legend([h1,h2],'Cluster 1','Cluster 2');
end
ylabel('R2*')

if Full_analysis == 'Y'
    
    % Significant test for zones
    % Generate entirely random image
    random_zones = repmat(zeros(size(BW)),[1,1,3]);
    rand_seq = round(rand(length(pixel_coords(:,1)),1));
    random_zones_1_pixels = pixel_coords(rand_seq==0,:);
    random_zones_2_pixels = pixel_coords(rand_seq==1,:);
    for i = 1:length(rand_seq(rand_seq == 0))
        random_zones(random_zones_1_pixels(i,2),random_zones_1_pixels(i,1),1) = 255;
    end
    for i = 1:length(rand_seq(rand_seq == 1))
        random_zones(random_zones_2_pixels(i,2),random_zones_2_pixels(i,1),2) = 255;
    end
    % imshow(random_zones)
    
    % %Compute correlation (ON/OFF) matrix
    % %Random zones
    dist_matrix = squareform(pdist(pixel_coords));
    col = repmat(rand_seq,1,length(rand_seq));
    row = repmat(rand_seq',length(rand_seq),1);
    corr_matrix = ~xor(col,row); %Similarity matrix
    corr_pixels = corr_matrix.*(dist_matrix > 0).*(dist_matrix < 1.5);  % Only immediate neighbours
    num_neighbours = sum((dist_matrix < 1.5).*(dist_matrix > 0));
    bin_test = sum(corr_pixels)./num_neighbours;
    all_bin_test = sum(sum(corr_pixels))/sum(num_neighbours);
    
    subplot(3,3,3)
    [i1,~] = find(corr_pixels > 0);
    corr_image = zeros(size(BW));
    for i = 1:size(pixel_coords,1)
        corr_image(pixel_coords(i,2),pixel_coords(i,1)) = bin_test(i);
    end
    corr_image_2D = reshape(corr_image,size(BW));
    surf(flipud(corr_image_2D));
    colormap jet
    shading interp
    view(2)
    axis('tight')
    title('prob. of success (random assignment)')
    colorbar('EastOutside')
    
    subplot(3,3,5)
    [n_out,~] = hist(bin_test,7); hist(bin_test,7); hold on
    stem(all_bin_test,max(n_out)*1.1,'r')
    xlabel('prob. of success')
    ylabel('frequency')
    
    %Actual zones
    col = repmat(idx-1,1,length(idx));
    row = repmat(idx'-1,length(idx),1);
    corr_matrix = ~xor(col,row);
    corr_pixels = corr_matrix.*(dist_matrix > 0).*(dist_matrix < 1.5);  % Only immediate neighbours
    num_neighbours = sum((dist_matrix < 1.5).*(dist_matrix > 0));
    bin_test = sum(corr_pixels)./num_neighbours;
    all_bin_test = sum(sum(corr_pixels))/sum(num_neighbours);
    
    subplot(3,3,4)
    [i1,~] = find(corr_pixels > 0);
    corr_image = zeros(size(BW));
    for i = 1:size(pixel_coords,1)
        corr_image(pixel_coords(i,2),pixel_coords(i,1)) = bin_test(i);
    end
    corr_image_2D = reshape(corr_image,size(BW));
    surf(flipud(corr_image_2D));
    colormap jet
    shading interp
    view(2)
    axis('tight')
    title('prob. of success (Kmeans)')
    colorbar('EastOutside')
    
    subplot(3,3,6)
    [n_out,~] = hist(bin_test,7); hist(bin_test,7); hold on
    stem(all_bin_test,max(n_out)*1.1,'r')
    xlabel('prob. of success')
    ylabel('frequency')
    
    pout = myBinomTest(sum(sum(corr_pixels)),sum(num_neighbours),0.5,'Greater');
    text(0.1,0.9*max(n_out),strcat('Sig. greater (P=',num2str(pout) ,')'))
    
    
    for i =1:12
        [~,~,Jk,D] = kmeans(data_vec',i,'Distance','sqEuclidean','replicates',10);
        D_min = min(D,[],2);
        sumJk(i) = sum(Jk);
        %     sumJk(i) = sum(D_min.^2);
    end
    subplot(3,3,[7])
    plot(sumJk)
    xlabel('num. of clusters')
    ylabel('total intra-cluster distance')
    axis('tight')
    
    % Use BIC
    N = size(data_vec,2);
    dim = size(data_vec,1);
    for i =1:12
        [idx,~,Jk,D] = kmeans(data_vec',i,'Distance','sqEuclidean','replicates',10);
        D_min = min(D,[],2);       % Find distance to closest cluster
        sigma2 =  1/(N - i)*sum(D_min);
        for j = 1:i
            Rn = length(idx(idx == j));
            l(j) =  -Rn/2*log(2*pi) - dim*Rn/2*log(sigma2) - (Rn - i)/2 + Rn*log(Rn) - Rn*log(N);
        end
        l_tot = sum(l);
        p = i*(dim) + 1;  % all free parameters, i*dim mean coordinates and i spherical variances
        BIC(i) = l_tot-p/2*log(N);
    end
    subplot(3,3,8); plot(BIC);
    xlabel('num. of clusters')
    ylabel('BIC')
    axis('tight')
    
    % Percentage of variance explained
    % Use Explained variance
    N = size(data_vec,2);
    dim = size(data_vec,1);
    % var_tot = sum(sum(data_vec.^2))/(N-1);
    mean_tot = mean(data_vec')';
    var_tot = sum(sum((data_vec - repmat(mean_tot,1,size(data_vec,2))).^2))/(N-1);
    var_tot = sum(var(data_vec'));
    clear sigma2
    for i =1:12
        [idx,cluster_c,Jk,D] = kmeans(data_vec',i,'Distance','sqEuclidean','replicates',10);
        D_min = min(D,[],2);
        var_kmeans = sum(D_min)/(N-1);
        var_expl(i) = (var_tot - var_kmeans)/var_tot;
    end
    subplot(3,3,9); plot(var_expl);
    xlabel('num. of clusters')
    ylabel('Explained variance')
    axis('tight')
    
    var_diff = diff(var_expl);
    gain = var_diff(2:end)./var_diff(1:end-1);
    num_clusters = min(find(gain < 0.5)) + 1;
    text(4,0.1,strcat('Numclusters = ',num2str(num_clusters)))
    
    pause(2)
    
    
end