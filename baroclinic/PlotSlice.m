% cd('./run9')
    writerObj = VideoWriter('BaroclinicSlice.avi');
    writerObj.FrameRate=5;
    open(writerObj);
    h=figure(1);
for t=1:size(rv1_rot,3)
    subplot(3,1,1)
    imagesc(real(squeeze(rv1(:,:,tp+t-1))'));
    title('Dynamics in Full State Space')
    subplot(3,1,2)
    imagesc(real(squeeze(rv1_rot(:,:,t))'));
    title('Dynamics in the Slice')
    subplot(3,1,3)
    plot(real(condition));
    title('Transversality condition??')
    ymax=1.1*max(real(condition));
    ymin=1.1*min(real(condition));
    line([t t],[ymin ymax]);
    axis([0 size(rv1_rot,3) ymin ymax]);
    a=1;
    %Create movie
    frame=getframe(h);
    writeVideo(writerObj,frame);   
end
close(writerObj);