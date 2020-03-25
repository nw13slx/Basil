#!/bin/env python

from bs4 import BeautifulSoup
from bs4.element import Tag
from os import getcwd, remove
from plotly.offline import plot

def iterprint(item, header=''):
    """
    print the tree structure of the html
    """

    if (item.name is None):
        return
    print(f"{header}-", item.name)
    runchild = False
    if (item.children is not None):
        if (len(list(item.children))>1):
            runchild = True
        elif (len(list(item.children))>0):
            child = list(item.children)[0]
            if (child.name is not None):
                runchild = True

    if (runchild):
        for child in item.children:
                iterprint(child, header=header+'| ')
        print(f"{header}*", item.name)

def html2markdown(fig, basename):
    """
    convert a html output from plotly to
    java scripts and lines that can be compile with pandoc

    1.   generate plotly.graph_objects.Figure object

    2.   and then convert it

    | plotlyhtml2markdown(figure, basename)

    3.   paste the output lines to markdown file and compile with pandoc

    | pandoc --from markdown-markdown_in_html_blocks+raw_html --wrap=auto \
            --mathjax --mathml -c pandoc.css -s --toc --toc-depth=3 $prefix.md -o $prefix.html

    Args:
        fig (plotly.graph_objects.Figure): plotly figure
        basename (str): prefix of output js scripts
    """

    plot(fig, filename=f'{basename}.html')
    with open(f"{basename}.html") as fin:
        soup = BeautifulSoup(fin, features="lxml")

    # iterprint(soup)

    nscript = 0
    print("<html>")
    for item in soup.body.div.children:
        if (item.name == "script"):
            # print(list(item.__dict__.keys()))
            if (len(item.contents[0])>50):
                with open(f"{basename}_{nscript}.js", "w+") as fout:
                    print(item.contents[0], file=fout)
                item.contents[0]=" "
                item.attrs["src"]=f"{getcwd()}/{basename}_{nscript}.js"
                nscript += 1
                print(item)
            else:
                print(item)
        elif (item.name is not None):
            print(item)
    print("</html>")

    remove(f'{basename}.html')

# layout setting for vline hover

vlinehover = dict(xaxis=dict(
                   showline=True,
                   showgrid=True,
                   showticklabels=True,
                   linewidth=2,
                   mirror=True,
                   showspikes = True,
                   spikemode  = 'across',
                   spikesnap = 'data',
                   spikecolor = 'black',
                   spikethickness = 0.5,
                   spikedash='dot'
               ),
               yaxis=dict(
                           showline=True,
                   showgrid=True,
                   showticklabels=True,
                   linewidth=2,
                   mirror=True,
               ),
               showlegend=True,
               hovermode='x',
               spikedistance =  -1,
               font=dict(
                   family="STIXGeneral",
                   size=18
                   )
              )
