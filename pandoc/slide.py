#!/usr/bin/env python

"""
Pandoc filter to convert plain text line "&slideshow prefix suffix n"
into html script for photo gallery,  accompanied with javascript slide

pandoc --filter slideshow.py --from markdown-markdown_in_html_blocks+raw_html --wrap=auto --mathml --mathjax -c pandoc.css -s --toc -N --toc-depth=3 $prefix.md --include-after-body slide  -o $prefix.html

the javascript should be

~~~~
<script>
 var slideIndex = [1,1];
 var slideId = ["mySlides1", "mySlides2"]
 showSlides(1, 0);
 showSlides(1, 1);

 function plusSlides(n, no) {
   showSlides(slideIndex[no] += n, no);
 }

 function showSlides(n, no) {
   var i;
   var x = "mySlides"+str(no+1); //document.getElementsByClassName(slideId[no]);
   if (n > x.length) {slideIndex[no] = 1}
   if (n < 1) {slideIndex[no] = x.length}
   for (i = 0; i < x.length; i++) {
     x[i].style.display = "none";
   }
   x[slideIndex[no]-1].style.display = "block";
 }
 function currentSlide(n, no) {
   slideIndex[no] = n;
   showSlides(n, no);
 }
</script>
~~~~

"""

from pandocfilters import toJSONFilter, Emph, Para, RawBlock


def behead(key, value, format, meta):
    if key == 'Para' :
        stop = False
        ids = 0
        markid = -1
        lv = len(value)
        while(not stop and ids < lv and markid <0):
            if (value[ids]['t'] == 'Str'):
                if (value[ids]['c'] == '&slideshow'):
                    markid = ids
            elif (value[ids]['t'] != 'Space'):
                stop = True
            ids += 1
        if (markid>=0):
            lv = len(value)
            prefix = value[markid+2]['c']
            suffix = value[markid+4]['c']
            n = int(value[markid+6]['c'])
            return RawBlock("html", create_slide_div(prefix, suffix, n))


def create_slide_div(prefix, suffix, n):
    global slide_counter
    slide_counter += 1
    headblock = ''' <div class="slideshow-container"> \n'''
    block1 = '''    <div class="mySlides{counter} fade">
      <div class="numbertext">{i} / {n}</div>
      <img src={prefix}{i}.{suffix} style="width:100%">
      <div class="text"> </div>
    </div> \n'''
    block2 = ''' <a class="prev" onclick="plusSlides(-1, {count})">❮</a>
    <a class="next" onclick="plusSlides(1, {count})">❯</a>
    </div>
    <div style="text-align:center"> \n'''
    block3 = '''<span class="dot" onclick="currentSlide({i}, {count})"></span> '''

    newstring = headblock
    for i in range(n):
        block = block1.format(i=i+1, n=n, prefix=prefix,
                suffix=suffix, counter=slide_counter)
        newstring = newstring + block

    newstring = newstring + block2.format(count=slide_counter-1)
    for i in range(n):
        block = block3.format(i=i+1, count=slide_counter-1)
        newstring = newstring + block
    newstring = newstring + '''
</div>
</div>
'''
    return newstring

if __name__ == "__main__":
  slide_counter=0
  toJSONFilter(behead)


