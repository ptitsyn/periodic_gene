# periodic_gene
Scripts for analysis of periodic patterns in gene expression, presumably circadian 

Previous publications based on these algorithms: 
True or false: all genes are rhythmic. Ptitsyn AA, Gimble JM. Ann Med. 2011 Feb;43(1):1-12. doi: 10.3109/07853890.2010.538078. Epub 2010 Dec 8. PMID: 21142579 
Ptitsyn AA, Zvonic S, Gimble JM. Permutation test for periodicity in short time series data. BMC Bioinformatics. 2006 Sep 6;7 Suppl 2(Suppl 2):S10. doi: 10.1186/1471-2105-7-S2-S10. PMID: 17118131; PMCID: PMC1683571.
Ptitsyn AA, Zvonic S, Gimble JM. Digital signal processing reveals circadian baseline oscillation in the majority of mammalian genes. PLoS Comput Biol. 2007 Jun;3(6):e120. doi: 10.1371/journal.pcbi.0030120. PMID: 17571920; PMCID: PMC1892608.
Ptitsyn AA, Zvonic S, Conrad SA, Scott LK, Mynatt RL, Gimble JM. Circadian clocks are resounding in peripheral tissues. PLoS Comput Biol. 2006 Mar;2(3):e16. doi: 10.1371/journal.pcbi.0020016. Epub 2006 Mar 10. PMID: 16532060; PMCID: PMC1391918.

I ported the old C/C++ code to Python. The lazy fat snake makes everything go painfully slow. Nevertheless, the overall time for the largest available data sets is still OK even for a cheap laptop. What's a couple of days compared to months if not longer to conduct the experiments in the lab to produce those data?

The code is commented to make it easy to adopt and modify, no manual is necessary.

