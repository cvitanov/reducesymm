a = size(x,1);
b = size(y,1);
c = size(z,1);

for i = 1:a
    for j = 1:b
        for k = 1:c
            plot3(x(i),z(k),y(j))
            hold on
        end
    end
end