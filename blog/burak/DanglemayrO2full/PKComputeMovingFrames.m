clear
clc

xhatp = [1; 1; 0; 0]; %Template point

MovingFrames(xhatp);

load movingframes.mat;

xhatGS = LieElement(pi/4,xhat);

save('gramschmidt.mat','xhatGS')

