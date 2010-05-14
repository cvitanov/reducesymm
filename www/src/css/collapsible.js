/*COLLAPSIBLE TEXT
  Use css and javascript to show/hidde text
  based on M. C. Matti's  http://www.webreference.com/programming/css_content/index.html

    <!-- javascript part    -->*/

    function dsp(loc){
      if(document.getElementById){
	  var foc=loc.firstChild;
	  foc=loc.firstChild.innerHTML?
	    loc.firstChild:
	    loc.firstChild.nextSibling;
	  foc.innerHTML=foc.innerHTML=='+'?'-':'+';
	  foc=loc.parentNode.nextSibling.style?
	    loc.parentNode.nextSibling:
	    loc.parentNode.nextSibling.nextSibling;
	  foc.style.display=foc.style.display=='block'?'none':'block';}}  

    if(!document.getElementById)
      document.write('<style type="text/css"><!--\n'+
	  '.dspcont{display:block;}\n'+
	  '//--></style>');
