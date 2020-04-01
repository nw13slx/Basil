#!/usr/bin/env python

"""
Pandoc filter to convert plain text line "&slideshow prefix suffix n"
into html script for photo gallery,  accompanied with javascript slide

pandoc --filter slideshow.py --from markdown-markdown_in_html_blocks+raw_html --wrap=auto --mathml --mathjax -c pandoc.css -s --toc -N --toc-depth=3 $prefix.md --include-after-body slide  -o $prefix.html

the javascript comes from enxaneta @ https://stackoverflow.com/a/52677453/8974431

~~~~
<script>
var slideIndex;

// the array of all containers
let containers = Array.from(document.querySelectorAll(".w3-display-container"));

// for each conteiner
containers.forEach(c =>{
  // get the array of images in this container
  let images = Array.from(c.querySelectorAll(".mySlides"));
  //the left button for this container
  let left = c.querySelector(".w3-display-left");
  //the right button for this container
  let right = c.querySelector(".w3-display-right");

  slideIndex = 0;//set the first slide
  plusDivs(0,images);


  // events for the this left and right buttons
  left.addEventListener("click", ()=>{plusDivs(-1,images)})
  right.addEventListener("click", ()=>{plusDivs(1,images)})
})


function showDivs(x) {
  if (slideIndex > x.length-1) {slideIndex = 0}
  if (slideIndex < 0) {slideIndex = x.length-1}

  //All the slides are display="none"
  for (let i = 0; i < x.length; i++) {
     x[i].style.display = "none";
  }
  // the current slide is display = "block";
  x[slideIndex].style.display = "block";

}

function plusDivs(n,x) {
  // increment the value for the slideIndex and show the slide
  slideIndex += n;
  showDivs(x);
}
</script>
~~~~

"""

from pandocfilters import toJSONFilter, Emph, Para, RawBlock

def slideshow(key, value, format, meta):
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
            nlist = []
            if (len(value)>=(markid+8)):
                for idx in range(markid+6, lv, 2):
                    nlist += [(value[idx]['c'])]
                n = -1
            return RawBlock("html", create_slide_div(prefix, suffix, n, nlist))


def create_slide_div(prefix, suffix, n, nlist):

    fout = open("local_write", "a+")
    print("start", file=fout)
    print(nlist, n, file=fout)

    if (n==0):
        return " "

    headblock = '''<div class="w3-display-container"> \n'''
    block1 = ''' <img class="mySlides" src="{prefix}{i}{suffix}" style="width:80%"> \n'''
    block2 = ''' <button class="w3-display-left">&#10094;</button>
    <button class="w3-display-right">&#10095;</button>
    </div>
    <div style="text-align:center"> \n </div> \n</div>'''

    if (n>0):
        nlist = range(1, n+1)

    print(nlist, n, file=fout)

    newstring = headblock
    for i in nlist:
        block = block1.format(i=i, prefix=prefix,
                suffix=suffix)
        newstring = newstring + block

    newstring = newstring + block2
    return newstring

if __name__ == "__main__":
  toJSONFilter(slideshow)


