clf();

plotgramschmidt(0,100,1,2,3,'b')

ah = gca();%Generate a handler for the current Axis object
set(ah,'box','off', 		%Remove the box.
	   'xtick',[0:-0.5:-2],
	   'xticklabel',{'0','','-1','','-2'},
	   'ytick',[-0.5:0.5:3],
	   'yticklabel',{'','0','','1','','2','','3'},
	   'ztick',[-0.2:0.2:0.8],
	   'zticklabel',{'','0','','0.2','','0.4','','0.8'}) 

print('rpoGS.tex','-S650,450',
'-depslatexstandalone',
'-tight',
'-F:Helvetica:10'
)
