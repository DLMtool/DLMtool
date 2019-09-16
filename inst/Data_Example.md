title: Example Data Documentation File
author: A. Hordyk <a.hordyk@oceans.ubc.ca>
date: 16 September 2019

## Metadata
This is an example Data Documentation file. 

The file can be edited with any normal text editor, and must be saved with a '.md' file extension.

The documentation file use Markdown syntax. This is mainly just regular text, but 
also allows easy inclusion of lists and hyper-links. For example, to create an
unordered list (i.e., bullet points) begin each new line with the * symbol:

* List item 1
* List item 2
  * List item 2a
  
Numbered lists can be created in a similar manner:

1. Item 1
2. Item 2


The ## Symbol indicates a second-level heading. Do not delete the second level headings
in the Data Documentation template file (e.g., ## Metadata etc). You can add sub-headings
in each section using the ### symbol. For example:

### A sub-section in the Metadata Section
Any text below this heading will appear in the sub-section. Do not add any additional first
or second level headings (i.e., don't use # or ##, only ###)

Links can be added in the following manner: [More information on Markdown Syntax](https://guides.github.com/pdfs/markdown-cheatsheet-online.pdf)

## Biology
The example Cobia data file has estimates of natural mortality, the von Bertalanffy
growth parameters, and the length-at-maturity parameters.

### Natural Mortality
Information relating to the estimated natural mortality, including links to relevant literature,
could be added here. 

### Growth Parameters
Text describing the data and methods used to obtain the estimates of the von Bertalanffy growth parameters
would be added here. Include figures, links, and/or references to data or supporting literature.

### Other Parameters
No information is available for length-weight parameters, or steepness and sigmaR of the 
stock-recruitment relationship.

## Selectivity
In this example we only have estimates of the length-at-first-capture (LFC). Selectivity 
is assumed to be asymptotic. 

## Time-Series
Information describing the time-series data should be added here. Use sub-sections (e.g., ### Heading)
to describe each time-series data set that is included in the data object.


## Catch-at-Age
This is example text for the catch-at-age data.


## Catch-at-Length
This is example text for the catch-at-length data.

## Reference
This is an example text for reference section. In this case we do not have any estimates 
for the Reference values for the Cobia example.


## Reference List
References can be added here, e.g.,:

Carruthers, T.R., and Hordyk, A.R. 2018. The data-limited methods toolkit (DLMtool): an R package for informing management of data-limited populations. Methods Ecol. Evol. 9(12): 2388–2395. doi:10.1111/2041-210X.13081.

