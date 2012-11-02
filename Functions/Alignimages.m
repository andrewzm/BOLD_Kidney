function [images_aligned,h,theta,xshift,yshift] = Alignimages(images_cropped,ratnum,daynum)

h=0;
theta=0;
images_aligned = images_cropped;
cIm1 = squeeze(images_cropped(ratnum,daynum,3,:,:));
cIm1(cIm1>1) = 1;
cIm1(cIm1<0) = 0;
for i = [1,2,4:19]
    cIm2 = squeeze(images_cropped(ratnum,daynum,i,:,:));
    cIm2(cIm2>1) = 1;
    cIm2(cIm2<0) = 0;
    if sum(sum(cIm2)) > 0
        [h,images_aligned(ratnum,daynum,i,:,:), theta(i),xshift(i),yshift(i),~,~] = im_reg_MI(cIm1*255,padarray(cIm2*255,[3,3]),[-3:0.1:3],1,'T');
        xshift(i) = xshift(i) -3 -1; yshift(i) = yshift(i) -3 -1; % because of padding
        %theta
    end
    close all
end
images_aligned(ratnum,daynum,3,:,:) = cIm1*255;
images_aligned(ratnum,daynum,:,:,:) = images_aligned(ratnum,daynum,:,:,:)/255;
theta
