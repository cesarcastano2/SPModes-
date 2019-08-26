function CV = getCV(x)


CV = nanstd(x)/nanmean(x);
     