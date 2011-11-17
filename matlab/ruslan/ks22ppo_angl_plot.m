
clear;
load ks22f90h25angl.mat


minangl=ppo(1).angl(1);
minangl_pos= [1, 1];

non_prob=0;

Ttab=[];
iTab=[]; % use it to store index ipo.
anglTab=[];
floqTab=[];
TtabProb=[];
floqTabProb=[];

knownRepeatsPPO=[3, 19, 50, 162, 576, 1592, 209];
knownRepeatsRPO=[2, 11, 31, 111, 305, 1016];

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
                    iTab = [iTab, ipo];
                    floqTab= [floqTab, log(abs(ppo(ipo).e(1)))/ppo(ipo).T];
                    anglTab = [anglTab, ppo(ipo).angl(i)/pi];
                    non_prob=non_prob+1;
                end
            end
        end
    else
        TtabProb =[TtabProb, ppo(ipo).T];
        floqTabProb= [floqTabProb, log(abs(ppo(ipo).e(1)))/ppo(ipo).T];
    end
end

save('ks22ppo_angl_T.dat', 'Ttab', '-ascii','-double','-tabs');
save('ks22ppo_angl_floq.dat', 'floqTab', '-ascii','-double','-tabs');
save('ks22ppo_angl_angl.dat', 'anglTab', '-ascii','-double','-tabs');
save('ks22ppo_angl_i.dat', 'iTab', '-ascii','-double','-tabs');

fig1=figure(); hold on;

plot(Ttab,anglTab,'k.','MarkerSize',6);
axis([0 1.01*max(Ttab) 0.9*min(anglTab) 1.1*max(anglTab)]);
% axis tight;
set(gca,'box','on', 'Layer','top');
set(gca, 'fontsize', 12)
xlabel('T_p','FontSize',14);
ylabel('\theta_{min} / \pi','FontSize',14);

period_edges=[min(Ttab):15:max(Ttab) max(Ttab)];
n_elements = histc(Ttab, period_edges);
c_elements = cumsum(n_elements);
anglTabC=zeros(size(c_elements,2),1);
anglTabC(1)=mean(anglTab(1:c_elements(1)));
for i=2:size(anglTabC), %
    anglTabC(i)=mean(anglTab(c_elements(i-1)+1:c_elements(i)));
end
stairs(period_edges,anglTabC,'r-','LineWidth', 2);
% plot(period_edges,anglTabC,'r.','MarkerSize', 6);

saveas(fig1, 'ks22ppo_min_angle.png');
saveTightFigure(fig1, 'ks22ppo_min_angle.pdf');
saveas(fig1, 'ks22ppo_min_angle.eps');

hold off;

%%% Only plot orbits with small period, small angle, 1st Floquet exponent
%%% close to Lyapunov exp.
fig11=figure();

compTab=[Ttab' anglTab' floqTab'];

excl1=Ttab>110;
excl2 = anglTab>0.01;
excl3 = or(floqTab>0.06,floqTab<0.035);
excl= [excl1'  excl2'  excl3'];

compTab(any(excl,2),:)=[];

plot3(compTab(:,1),compTab(:,2),compTab(:,3),'k.','MarkerSize',6);
xlabel('T_p','FontSize',14);
ylabel('\theta_{min} / \pi','FontSize',14);
zlabel('\mu^{(p,1)}');

%%%%

fig1a=figure(); hold on;

plot(floqTab,anglTab,'k.','MarkerSize',6);
axis([0 1.05*max(floqTab) 0.9*min(anglTab) 1.1*max(anglTab)]);
% axis tight;
set(gca,'box','on', 'Layer','top');
set(gca, 'fontsize', 12)
xlabel('\mu_p^{(1)}','FontSize',14);
ylabel('\theta_{min} / \pi','FontSize',14);

saveas(fig1a, 'ks22ppo_min_angle_floq.png');
saveTightFigure(fig1a, 'ks22ppo_min_angle_floq.pdf');
saveas(fig1a, 'ks22ppo_min_angle_floq.eps');

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


fig3=figure(); hold on;

plot(TtabProb,floqTabProb,'.');

figure();

edges=0:0.01:0.25;
n_elements_prob = histc(floqTabProb, edges);
hist(floqTabProb,edges);

fig4=figure(); hold on;

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

