  function wysiwyg()
%  function wysiwyg() - What You See Is What You Get!

  units=get(gcf,'units');
  set(gcf,'units',get(gcf,'PaperUnits'));
  set(gcf,'Position',get(gcf,'PaperPosition'));
  set(gcf,'Units',units);
