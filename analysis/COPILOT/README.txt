mt.threshold = 20
filtering.ratio = 1
top.percent = different cutoff corresponding to censoring UMI>38041 (which corresponds to 1% cells having the most UMI combining all libs)
min.UMI.low.quality = 10
min.UMI.high.quality = 25
Dot in gene_ids AND mt_genes

Note:
1) removing dot from gene_ids AND mt_genes does NOT change the result
2) removing dot from ONLY mt_genes also does NOT change the result (maybe COPILOT does "partial" str_detect)

UseIt is FilRatio1TopPercVariesMinUmi1025MtThresh20
