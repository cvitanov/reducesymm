load ks22f90h25angl.mat

fig1=figure(); hold on;

minangl=ppo(1).angl(1);
minangl_pos= [1, 1];

non_prob=0;

for ipo=1:size(ppo,2),
    if not(ppo(ipo).angl_prob),
        if ipo~=3 && ipo ~=19. % exclude known repeats (delete from dataset?)
            for i=1:size(ppo(ipo).angl)
                if imag(ppo(ipo).angl) ==0
                    if ppo(ipo).angl(i)<minangl,
                        minangl=ppo(ipo).angl(i);
                        minangl_pos= [ipo, i];
                    end
                    plot(ppo(ipo).T, ppo(ipo).angl(i),'b.');
                    non_prob=non_prob+1;
                end
            end
        end
    end
end

xlabel('T');
ylabel('min. angle');

hold off;

fig2=figure(); hold on;

for ipo=1:size(ppo,2),
    if not(ppo(ipo).angl_prob),
        if ipo~=3 && ipo ~=19.
            for i=1:size(ppo(ipo).angl)
                if imag(ppo(ipo).angl) ==0
                    plot(ipo, ppo(ipo).angl(i),'b.');
                end
            end
        end
    end
end

xlabel('i');
ylabel('min. angle');

