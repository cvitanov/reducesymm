function result = ChartBorderCondition(xhatstar, xhatp)

load('generator.mat'); %Load the Lie element generator

%Template tangent
tp = T * xhatp;
%Group tangent for xstar
tstar = T * xstar;
