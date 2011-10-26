% load ks22f90h25angl.mat

fig1=figure(); hold on;

minangl=ppo(1).angl(1);
minangl_pos= [1, 1];

non_prob=0;

Ttab=[];
anglTab=[];

for ipo=1:size(ppo,2),
    if not(ppo(ipo).angl_prob),
        if ipo~=3 && ipo ~=19. % exclude known repeats (delete from dataset?)
            for i=1:size(ppo(ipo).angl,2)
                if imag(ppo(ipo).angl(i))==0
                    if ppo(ipo).angl(i)<minangl,
                        minangl=ppo(ipo).angl(i)/pi;
                        minangl_pos= [ipo, i];
                    end
%                     plot(ppo(ipo).T, ppo(ipo).angl(i)/pi,'k.','MarkerSize',6);
                    Ttab =[Ttab, ppo(ipo).T];
                    anglTab = [anglTab, ppo(ipo).angl(i)/pi];
                    non_prob=non_prob+1;
                end
            end
        end
    end
end

plot(Ttab,anglTab,'k.','MarkerSize',6);
axis([0 1.01*max(Ttab) 0.9*min(anglTab) 1.1*max(anglTab)]);
% axis tight;
set(gca,'box','on', 'Layer','top');
set(gca, 'fontsize', 12)
xlabel('T_p','FontSize',14);
ylabel('\theta_{min} / \pi','FontSize',14);

saveTightFigure(fig1, 'ks22ppo_min_angle.pdf');
saveas(fig1, 'ks22ppo_min_angle.eps');

hold off;

fig2=figure();
angl_edges=0:0.0001:0.32;
% angl_edges=[0,0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5];
n_elements = histc(anglTab, angl_edges)/size(anglTab,2);
loglog(angl_edges,n_elements,'k-','MarkerSize', 6);
c_elements = cumsum(n_elements);
% axis tight;

set(gca, 'fontsize', 12);

ylabel('PDF','FontSize',14);
xlabel('\theta_{min} / \pi','FontSize',14);



saveTightFigure(fig2, 'ks22ppo_min_angle_pdf.pdf');
saveas(fig2, 'ks22ppo_min_angle_pdf.eps');


% fig3=figure(); hold on;
% 
% for ipo=1:size(ppo,2),
%     if not(ppo(ipo).angl_prob),
%         if ipo~=3 && ipo ~=19.
%             for i=1:size(ppo(ipo).angl)
%                 if imag(ppo(ipo).angl) ==0
%                     plot(ipo, ppo(ipo).angl(i),'b.');
%                 end
%             end
%         end
%     end
% end
% 
% xlabel('i');
% ylabel('min. angle');

