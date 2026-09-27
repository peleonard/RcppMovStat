# RcppMovStat
## Background

Developed by Pengfei Li as a side product of his Ph.D. research in the Department of Mathematics
at Tulane University, advised by Dr. Michelle Lacey. The dissertation, *Statistical Models for the
Analysis of Whole Genome Bisulfite Sequencing Datasets* (2020), developed statistical methods for
DNA methylation data.

## Intro
This package provides several efficient functions to 
    calculate common moving (or rolling, running) statistics for both evenly and unevenly spaced time 
    series: moving average, moving median, moving maximum (minimum), and so on. Built on 'C++', these 
    functions are apparently more efficient than those written in a traditional 'R' way and also faster
    than others using package 'Rcpp'. 
    
Install me with devtools:

    install_github("peleonard/RcppMovStat")
