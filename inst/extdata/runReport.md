Run Report
==========

Analysis Parameters
-------------------


```r
alpha = 0.05
```





Global Statistics
-----------------



```r
cuff
```

```
## CuffSet instance with:
## 	 3 samples
## 	 400 genes
## 	 1203 isoforms
## 	 662 TSS
## 	 906 CDS
## 	 1062 promoters
## 	 1986 splicing
## 	 990 relCDS
```




### FPKM distributions


```r
csDensity(genes(cuff))
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-31.png) 

```r
csDensity(isoforms(cuff))
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-32.png) 


Significant features
--------------------

### Genes


```r
sigMatrix(cuff, level = "genes")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 




```r
sigGeneIDs <- getSig(cuff, alpha = alpha)
sigGenes <- getGenes(cuff, sigGeneIDs)
```



There are `207` significantly different genes at a `5`$\%$ FDR. 

#### Visualizations


```r
csHeatmap(sigGenes, cluster = "both", labRow = F)
```

```
## Using tracking_id, sample_name as id variables
```

```
## Using as id variables
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 



### Isoforms


```r
sigMatrix(cuff, level = "isoforms")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 

