nx = 15;
ny = 10;
nz = 15;
c = [1 1 1];
J = [];
for i = 0:nx
    for j = 0:ny
        for k = 0:nz      
            M{i+1,j+1,k+1} = c;
            J = [J;c];
            c = c + 1;
        end
    end 
end

tic
u1 = [];
for i = 0:nx
    for j = 0:ny
        for k = 0:nz
           u1 = [u1; M{i+1,j+1,k+1}];
        end
    end
end
toc

tic
count = 1;
u2 = [];
for i = 0:nx
    for j = 0:ny
        for k = 0:nz
           u2 = [u2; J(count,:)];
           count = count + 1;
        end
    end
end
toc