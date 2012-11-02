function Animate(images)

for i = 1:19
    I = squeeze(double(images(i,:,:)));
    surf(flipud(I))
    caxis([0.2 1.4])
    colormap jet
    view(2)
    shading interp
    axis('tight')
    title(['scan number ',num2str(i)])
    
    pause(1)
    
end


