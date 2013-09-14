function chartborder = FindChartBorder(xhatp)

load('generator.mat'); %Load the Lie element generator

%Template tangent
tp = T * xhatp;

