#	Custom graphics for Rank-based Gene Ontology Analysis with Adaptive Clustering (Wright et al. 2015)

Custom by Pierre-Louis STENGER - 06-02-2019
Pierrelouis.stenger@gmail.com
phd student  

All code is from https://github.com/z0on/GO_MWU for Rank-based Gene Ontology Analysis with Adaptive Clustering.

Here I just custom the gomwuPlot.

In order to: i) don't show the genes ("18/45" for example); ii) reverse the graph in other side (cluster tree at the right) in order to arrange later many graphics, iii) combine the both.

So there 3 customized function from the incredible work of Mikhail V Matz (https://github.com/z0on/GO_MWU):
```
gomwuPlot_reverse_without_genes()
gomwuPlot_reverse()
gomwuPlot_without_genes()
```

Normal plot from the gomwuPlot()
-----------
![alt tag](https://github.com/PLStenger/RBGOA_custom_graphics/blob/master/normal.png)

Custom plot gomwuPlot_reverse()
-----------
![alt tag](https://github.com/PLStenger/RBGOA_custom_graphics/blob/master/gomwuPlot_reverse.png)


Custom plot gomwuPlot_without_genes()
-----------
![alt tag](https://github.com/PLStenger/RBGOA_custom_graphics/blob/master/gomwuPlot_reverse_without_genes.png)


Custom plot gomwuPlot_reverse_without_genes()
-----------
![alt tag](https://github.com/PLStenger/RBGOA_custom_graphics/blob/master/gomwuPlot_without_genes.png)


Source:
- https://github.com/z0on/GO_MWU
- Wright, R. M., Aglyamova, G. V., Meyer, E.  and Matz, M. V. Gene expression associated with white syndromes in a reef-building coral, Acropora hyacinthus. BMC Genomics 2015, 16: 371. 
( http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1540-2 )
