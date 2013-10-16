function returnmap(filename)

%filename: file including the Poincare section data ps(1:size(ps,1)-1,:)
%		   ordered in time t=ps(size(ps,1), :)
%V: holds the Gram-Schmidt basis for the Poincare section, where
%V(:,size(V,2)) = hyperplane normal

load(filename)
xhatp = [xhatp;	
		0]; %Attach a 0 so that 

n = size(ps,1) - 1; %get the dimensions

psectrelative = ps(1:n+1, :) - xhatp;

psectrelativedistances(i) = norm(ps(1:n,:));

[dummy, index] =  sort(psectrelative(3,:));
psectsortlength = psectrelative(:,index)

s(1) = 0; %Arclength from initial position at the curve;

for i = 2:size(psectsortlength, 2)
	
	dx = psectsortlength(1:n,i) - psectsortlength(1:n,i-1);
	
	s(i)=s(i-1)+sqrt(dot(dx,dx));
	
end

psectsortlength = [psectsortlength;
				   s];

[dummy, index] = sort(psectsortlength(n+1,:));
psectsorttime = psectsortlength(:, index);

sn=psectsorttime(n+2,1:size(psectsorttime,2)-1);
snplus1=psectsorttime(n+2,2:size(psectsorttime,2));

plot(sn,snplus1,'.', 'markersize', 5);
xlabel('$s_n$')
ylabel('$s_{n+1}$')

save('-append', filename, 'sn', 'snplus1');

hold on

plot(min(sn):max(sn)/20:max(sn), min(sn):max(sn)/20:max(sn), 'k')
