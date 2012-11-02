%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kidney_Analysis.m
%
% Sample code for "Anatomically unbiased analysis of renal BOLD magnetic 
%  resonance images reveals disruption of corticomedullary gradient during 
%  chronic angiotensin II infusion" by Robert I. Menzies, Andrew Zammit-Mangion, 
%  Lyam Hollis, Ross Lennen, Maurits A. Jansen, David J Webb, John J Mullins, 
%  James W. Dear, Guido Sanguinetti and Matthew A. Bailey
%
% Author:  Andrew Zammit-Mangion
% Date: 2 November 2012
% Requirements: MATLAB 7.14.0.739 (R2012a)
%               Statistical Toolbox
%               Image Processing Toolbox 
%               ./Functions
%
% Note: im_reg_MI.m and dependencies of this file were originally written 
% by Kateryna Artyushkova and were edited by Andrew Zammit-Mangion for this work.
%
% Description: This program considers only a subset of the entire data set;
% one rat pre-ANGII infusion (-3 days) and post-ANGII infusion (+3 days).
% The scans can be viewed in raw form, aligned and analysed post
% clustering. Since this program only considers a subset of the data
% outlier detection is not carried out. Moreover the quadrant is assumed
% to have been selected a priori.
%
% Please answer 'Y' or 'N' to the questions on starting this program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

% The scaling was applied to the T2* maps to improve conditioning. This is
% re-applied post-analysis.
global scaling
scaling = 40;       

addpath('./Functions')

disp('Analyzing K-means ROI sections on quadrant. Answer Y/N')
disp('Raw view of aligned scans?')
viewscans = upper(input('','s'));
disp('Use aligned images(''N'' will start re-alignment)?')
aligned = upper(input('','s'));

load('cropped_images'); 
     
% Raw view of scans
if viewscans == 'Y'
     disp('Showing scans for rat 2 day 3...')
     ratnum = 2; daynum = 3;
     Animate(squeeze(images_cropped(ratnum,daynum,:,:,:)))
end

% Alignment
if aligned == 'Y'
    disp('Loading aligned images...')
    for ratnum = 2
        for daynum = [2,3]
            load(['./Alignment/Set',num2str(ratnum),num2str(daynum)],'Aligned_set')  % Contains all aligned impage
            images_cropped(ratnum,daynum,:,:,:) = Aligned_set;
        end
    end
else
    disp('Commencing alignment to third scan...')
    for ratnum = 2
        for daynum = 2:3
            load('cropped_images');
            [images_aligned,h,theta,xshift,yshift] = Alignimages(images_cropped,ratnum,daynum);
            Aligned_set = images_aligned(ratnum,daynum,:,:,:);
            images_cropped(ratnum,daynum,:,:,:) = Aligned_set; %#ok<*SAGROW>
            save(['./Alignment/Set',num2str(ratnum),num2str(daynum)],'Aligned_set','theta','xshift','yshift');
        end
    end
end

% Clustering
for ratnum = 2
     disp('Commencing clustering...')
     for daynum = [2,3]
           load(['./Segmentation/Partial_masks',num2str(ratnum),num2str(daynum)]);
           [~, ~] = ExtractROIs_Kmeans(squeeze(images_cropped(ratnum,daynum,:,:,:)),BW,m,'F','Y');
           [zone_image, centroids] = ExtractROIs_Kmeans(squeeze(images_cropped(ratnum,daynum,:,:,:)),BW,m,'F','N');
           title(strcat('Rat ',num2str(ratnum),' Day ',num2str(daynum)));
      end
 end
 



