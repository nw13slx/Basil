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
from os import getcwd

def plotlyhtml2markdown(basename):
    with open(f"{basename}.html") as fin:
        soup = BeautifulSoup(fin, features="lxml")

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
                print(item)
                # print(len(item.contents))
            else:
                print(item)
        elif (item.name is not None):
            print(item)
    print("</html>")
