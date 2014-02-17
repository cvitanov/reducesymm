function y = KSm(t, x)
          %y = [10 * (x(2) - x(1));
               %x(1) * (28 - x(3));
               %x(1) * x(2) - 8/3 * x(3)];

xstr = mat2str(x); %Convert x vector to string
xstr = strrep(xstr, ';', ','); %Replace ';  ', with ', '

command = strcat('python -c "import KS; import numpy; numpy.set_printoptions(precision=15, linewidth=18); print KS.vfullssp(', xstr, ', 0);"');

[x,y]=system(command);

y = str2num(y(2:length(y)-2))';
N = size(y,1)*size(y,2);
y = reshape(y,1,N)';
