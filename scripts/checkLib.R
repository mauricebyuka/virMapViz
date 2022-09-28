#! /usr/bin/Rscript

programs <- c('Gviz', 'data.table', 'GenomicFeatures', 'GenomicRanges', 'tools')

ChekLib <- function (libraries) 
{
  #notInstalled <- c()
  for (prog in libraries)
  {
    all_installed <- rownames(installed.packages())
    if (prog %in% all_installed) 
    {
      cat(paste(' R library ', prog, ': OK ---', '\n', sep = ''))
    } else 
    {
      cat(paste(' R library ', prog, ': Library not found', '\n', sep = ''))
      #notInstalled <- c(notInstalled, prog)
    }
  }
  #cat(notInstalled)
}

ChekLib(programs)
