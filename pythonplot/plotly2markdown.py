#!/bin/env python
"""
convert a html output from plotly to
java scripts and lines that can be compile with pandoc

1.   use this to output plots to html

| plotly.offline.plot(fig, filename=f'{k}.html')

2.   and then convert it

| plotlyhtml2markdown(f'{k}.html')

3.   paste the output lines to markdown file and compile with pandoc

| pandoc --from markdown-markdown_in_html_blocks+raw_html --wrap=auto --mathjax --mathml -c pandoc.css -s --toc --toc-depth=3 $prefix.md -o $prefix.html

"""

from bs4 import BeautifulSoup
from bs4.element import Tag
from os import getcwd

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

def plotlyhtml2markdown(basename):
    with open(f"{basename}.html") as fin:
        soup = BeautifulSoup(fin, features="lxml")

    iterprint(soup)

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
