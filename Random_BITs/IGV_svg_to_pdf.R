install.packages('rsvg')

library(rsvg)
svg = './Destkop/Figures/CoCLIP/GeneTrack'

rsvg_pdf('SLK.svg', file = 'SLK.pdf', width = 1920, height = 1280, css = NULL)
