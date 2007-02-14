To reproduce poster.pdf file:

1. First latex text.tex to produce all text in the poster seperated in pages. Use:

 ./update

2. Split text.dvi to seperate eps documents with dvi_split script:

 ./dvi_split

In figs/ directory we should now have text$i.eps files that contain text fragments.

3. Open xfig:

 xfig poster.fig &

4. Export the poster in pdf format using xfig menus.
