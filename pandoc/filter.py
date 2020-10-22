#!/usr/bin/env python

from pandocfilters import toJSONFilter, Image, Para, RawBlock, Str

def slideshow(key, value, format, meta):
    """
    Pandoc filter to convert plain text line "&slideshow prefix suffix n"
    into html script for photo gallery,  accompanied with javascript slide

    pandoc --filter slideshow.py --from markdown-markdown_in_html_blocks+raw_html --wrap=auto --mathml --mathjax -c pandoc.css -s --toc -N --toc-depth=3 $prefix.md --include-after-body slide  -o $prefix.html

    the javascript comes from enxaneta @ https://stackoverflow.com/a/52677453/8974431

    ~~~~
    <script>

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

      plusDivs(0,images);


      // events for the this left and right buttons
      left.addEventListener("click", ()=>{plusDivs(-1,images)})
      right.addEventListener("click", ()=>{plusDivs(1,images)})
    })


    function showDivs(x, idx) {
      if (idx > x.length-1) {idx = 0}
      if (idx < 0) {idx = x.length-1}

      //All the slides are display="none"
      for (let i = 0; i < x.length; i++) {
         x[i].style.display = "none";
      }
      // the current slide is display = "block";
      x[idx].style.display = "block";

    }

    function plusDivs(n,x) {

      //find the current one
      var idx = 0;
      for (let i = 0; i < x.length; i++) {
          if ( x[i].style.display == "block"){
              idx = i
          }
      }
      // increment the value for the slideIndex and show the slide
      idx += n;
      showDivs(x, idx);
    }
    </script>
    ~~~~

    """

    if key in ['Para', 'Plain'] :
        stop = False
        ids = 0
        markid = -1
        lv = len(value)
        case = {'&slideshow':1, '&loopimage':2}
        caseid = -1
        while(not stop and ids < lv and markid <0):
            if (value[ids]['t'] == 'Str'):
                if (value[ids]['c'] in case):
                    markid = ids
                    caseid = case[value[ids]['c']]
            elif (value[ids]['t'] != 'Space'):
                stop = True
            ids += 1
        if (markid>=0):
            lv = len(value)
            prefix = value[markid+2]['c']
            suffix = value[markid+4]['c']
            if (caseid == 1):
                percent = value[markid+6]['c']
                nlist = []
                if (len(value)>=(markid+10)):
                    for idx in range(markid+8, lv, 2):
                        nlist += [(value[idx]['c'])]
                    n = -1
                else:
                    n = int(value[markid+8]['c'])
                return RawBlock("html", create_slide_div(prefix, suffix, percent, n, nlist))
            elif (caseid == 2):
                nlist = []
                if (len(value)>=(markid+8)):
                    for idx in range(markid+6, lv, 2):
                        nlist += [(value[idx]['c'])]
                    n = -1
                else:
                    n = int(value[markid+6]['c'])
                return Para(create_image_div(prefix, suffix, n, nlist))

def create_image_div(prefix, suffix, n, nlist):

    if (n==0):
        return " "
    elif (n>0):
        nlist = range(1, n+1)
        alln = n
    else:
        alln = "-"

    value = []
    for i in nlist:
        value += [Image(["",[],[]],[Str(f"{i}")],[f"{prefix}{i}{suffix}","fig:"])]
    return value


def create_slide_div(prefix, suffix, percent, n, nlist):

    fout = open("local_write", "a+")
    print("start", file=fout)
    print(nlist, n, file=fout)

    if (n==0):
        return " "

    headblock = '''<div class="w3-display-container"> \n'''
    block1 = '''<div class="mySlides"> <figure> '''\
            '''<img src="{prefix}{i}{suffix}" style="width:{percent}%">'''\
            '''<figcaption> {i}/{n} </figcaption>'''\
            ''' </figure></div>\n'''
    block2 = ''' <button class="w3-display-first">&#10094;&#10094;</button>
    <button class="w3-display-left">&#10094;</button>
    <button class="w3-display-plays">&#9654;</button>
    <button class="w3-display-right">&#10095;</button>
    <button class="w3-display-last">&#10095;&#10095;</button>
    </div>
    <div style="text-align:center"> \n </div> \n</div>'''

    if (n>0):
        nlist = range(1, n+1)
        alln = n
    else:
        alln = "-"

    print(nlist, n, file=fout)

    newstring = headblock
    for i in nlist:
        block = block1.format(i=i, prefix=prefix,
                suffix=suffix, percent=percent,
                n=alln)
        newstring = newstring + block

    newstring = newstring + block2
    return newstring

if __name__ == "__main__":
  toJSONFilter(slideshow)


