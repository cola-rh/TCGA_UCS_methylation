cola Report for Hierarchical Partitioning - 'TCGA_UCS_methylation'
==================

**Date**: 2021-07-22 16:09:49 CEST, **cola version**: 1.9.4

----------------------------------------------------------------

<style type='text/css'>

body, td, th {
   font-family: Arial,Helvetica,sans-serif;
   background-color: white;
   font-size: 13px;
  max-width: 800px;
  margin: auto;
  margin-left:210px;
  padding: 0px 10px 0px 10px;
  border-left: 1px solid #EEEEEE;
  line-height: 150%;
}

tt, code, pre {
   font-family: 'DejaVu Sans Mono', 'Droid Sans Mono', 'Lucida Console', Consolas, Monaco, 

monospace;
}

h1 {
   font-size:2.2em;
}

h2 {
   font-size:1.8em;
}

h3 {
   font-size:1.4em;
}

h4 {
   font-size:1.0em;
}

h5 {
   font-size:0.9em;
}

h6 {
   font-size:0.8em;
}

a {
  text-decoration: none;
  color: #0366d6;
}

a:hover {
  text-decoration: underline;
}

a:visited {
   color: #0366d6;
}

pre, img {
  max-width: 100%;
}
pre {
  overflow-x: auto;
}
pre code {
   display: block; padding: 0.5em;
}

code {
  font-size: 92%;
  border: 1px solid #ccc;
}

code[class] {
  background-color: #F8F8F8;
}

table, td, th {
  border: 1px solid #ccc;
}

blockquote {
   color:#666666;
   margin:0;
   padding-left: 1em;
   border-left: 0.5em #EEE solid;
}

hr {
   height: 0px;
   border-bottom: none;
   border-top-width: thin;
   border-top-style: dotted;
   border-top-color: #999999;
}

@media print {
   * {
      background: transparent !important;
      color: black !important;
      filter:none !important;
      -ms-filter: none !important;
   }

   body {
      font-size:12pt;
      max-width:100%;
   }

   a, a:visited {
      text-decoration: underline;
   }

   hr {
      visibility: hidden;
      page-break-before: always;
   }

   pre, blockquote {
      padding-right: 1em;
      page-break-inside: avoid;
   }

   tr, img {
      page-break-inside: avoid;
   }

   img {
      max-width: 100% !important;
   }

   @page :left {
      margin: 15mm 20mm 15mm 10mm;
   }

   @page :right {
      margin: 15mm 10mm 15mm 20mm;
   }

   p, h2, h3 {
      orphans: 3; widows: 3;
   }

   h2, h3 {
      page-break-after: avoid;
   }
}
</style>




## Summary



First the variable is renamed to `res_rh`.


```r
res_rh = rh
```



The partition hierarchy and all available functions which can be applied to `res_rh` object.


```r
res_rh
```

```
#> A 'HierarchicalPartition' object with 4 combinations of top-value methods and partitioning methods.
#>   On a matrix with 375356 rows and 57 columns.
#>   Performed in total 14000 partitions.
#>   There are 7 groups under the following parameters:
#>     - min_samples: 6
#>     - mean_silhouette_cutoff: 0.9
#>     - min_n_signatures: 1000 (signatures are selected based on:)
#>       - fdr_cutoff: 0.05
#>       - group_diff: 0.25
#> 
#> Hierarchy of the partition:
#>   0, 57 cols
#>   |-- 01, 20 cols, 4690 signatures
#>   |   |-- 011, 7 cols (b)
#>   |   |-- 012, 6 cols (b)
#>   |   `-- 013, 7 cols (b)
#>   |-- 02, 18 cols, 5323 signatures
#>   |   |-- 021, 6 cols (b)
#>   |   |-- 022, 9 cols (b)
#>   |   `-- 023, 3 cols (b)
#>   `-- 03, 19 cols, 856 signatures (c)
#> Stop reason:
#>   b) Subgroup had too few columns.
#>   c) There were too few signatures.
#> 
#> Following methods can be applied to this 'HierarchicalPartition' object:
#>  [1] "all_leaves"            "all_nodes"             "cola_report"           "collect_classes"      
#>  [5] "colnames"              "compare_signatures"    "dimension_reduction"   "functional_enrichment"
#>  [9] "get_anno_col"          "get_anno"              "get_children_nodes"    "get_classes"          
#> [13] "get_matrix"            "get_signatures"        "is_leaf_node"          "max_depth"            
#> [17] "merge_node"            "ncol"                  "node_info"             "node_level"           
#> [21] "nrow"                  "rownames"              "show"                  "split_node"           
#> [25] "suggest_best_k"        "test_to_known_factors" "top_rows_heatmap"      "top_rows_overlap"     
#> 
#> You can get result for a single node by e.g. object["01"]
```

The call of `hierarchical_partition()` was:


```
#> hierarchical_partition(data = mat, top_n = 1000, top_value_method = c("SD", "ATC"), 
#>     partition_method = c("kmeans", "skmeans"), subset = 500, group_diff = 0.25, min_n_signatures = 1000, 
#>     filter_fun = function(mat) {
#>         s = rowSds(mat)
#>         order(-s)[1:30000]
#>     }, max_k = 8, scale_rows = FALSE, cores = 4)
```

Dimension of the input matrix:


```r
mat = get_matrix(res_rh)
dim(mat)
```

```
#> [1] 375356     57
```

All the methods that were tried:


```r
res_rh@param$combination_method
```

```
#> [[1]]
#> [1] "SD"     "kmeans"
#> 
#> [[2]]
#> [1] "ATC"    "kmeans"
#> 
#> [[3]]
#> [1] "SD"      "skmeans"
#> 
#> [[4]]
#> [1] "ATC"     "skmeans"
```

### Density distribution

The density distribution for each sample is visualized as one column in the following heatmap.
The clustering is based on the distance which is the Kolmogorov-Smirnov statistic between two distributions.




```r
library(ComplexHeatmap)
densityHeatmap(mat, ylab = "value", cluster_columns = TRUE, show_column_names = FALSE,
    mc.cores = 1)
```

![plot of chunk density-heatmap](figure_cola/density-heatmap-1.png)



Some values about the hierarchy:


```r
all_nodes(res_rh)
```

```
#>  [1] "0"   "01"  "011" "012" "013" "02"  "021" "022" "023" "03"
```

```r
all_leaves(res_rh)
```

```
#> [1] "011" "012" "013" "021" "022" "023" "03"
```

```r
node_info(res_rh)
```

```
#>     id best_method depth best_k n_columns n_signatures p_signatures is_leaf
#> 1    0  ATC:kmeans     1      3        57        17040      0.04540   FALSE
#> 2   01  ATC:kmeans     2      3        20         4690      0.01249   FALSE
#> 3  011 not applied     3     NA         7           NA           NA    TRUE
#> 4  012 not applied     3     NA         6           NA           NA    TRUE
#> 5  013 not applied     3     NA         7           NA           NA    TRUE
#> 6   02 ATC:skmeans     2      3        18         5323      0.01418   FALSE
#> 7  021 not applied     3     NA         6           NA           NA    TRUE
#> 8  022 not applied     3     NA         9           NA           NA    TRUE
#> 9  023 not applied     3     NA         3           NA           NA    TRUE
#> 10  03  ATC:kmeans     2      2        19          856      0.00228    TRUE
```

In the output from `node_info()`, there are the following columns:

- `id`: The node id.
- `best_method`: The best method selected.
- `depth`: Depth of the node in the hierarchy.
- `best_k`: Best number of groups of the partition on that node.
- `n_columns`: Number of columns in the submatrix.
- `n_signatures`: Number of signatures with the `best_k`.
- `p_signatures`: Proportion of hte signatures in total number of rows in the matrix.
- `is_leaf`: Whether the node is a leaf.

Labels of nodes are encoded in a special way. The number of digits
correspond to the depth of the node in the hierarchy and the value of the
digits correspond to the index of the subgroup in the current node, E.g. a label
of “012” means the node is the second subgroup of the partition which is the
first subgroup of the root node.

### Suggest the best k



Following table shows the best `k` (number of partitions) for each node in the
partition hierarchy. Clicking on the node name in the table goes to the
corresponding section for the partitioning on that node.

[The cola vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13)
explains the definition of the metrics used for determining the best
number of partitions.



```r
suggest_best_k(res_rh)
```


|Node              |Best method                                         |Is leaf   |Best k |1-PAC |Mean silhouette |Concordance | #samples|   |
|:-----------------|:---------------------------------------------------|:---------|:------|:-----|:---------------|:-----------|--------:|:--|
|[Node0](#Node0)   |ATC:kmeans                                          |          |3      |1.00  |0.98            |0.99        |       57|** |
|[Node01](#Node01) |ATC:kmeans                                          |          |3      |1.00  |0.98            |0.98        |       20|** |
|Node011-leaf      |<span style='color:grey;'><i>not applied</i></span> |✓ (b)     |       |      |                |            |        7|   |
|Node012-leaf      |<span style='color:grey;'><i>not applied</i></span> |✓ (b)     |       |      |                |            |        6|   |
|Node013-leaf      |<span style='color:grey;'><i>not applied</i></span> |✓ (b)     |       |      |                |            |        7|   |
|[Node02](#Node02) |ATC:skmeans                                         |          |7      |0.95  |0.72            |0.95        |       18|** |
|Node021-leaf      |<span style='color:grey;'><i>not applied</i></span> |✓ (b)     |       |      |                |            |        6|   |
|Node022-leaf      |<span style='color:grey;'><i>not applied</i></span> |✓ (b)     |       |      |                |            |        9|   |
|Node023-leaf      |<span style='color:grey;'><i>not applied</i></span> |✓ (b)     |       |      |                |            |        3|   |
|Node03-leaf       |ATC:kmeans                                          |✓ (&#99;) |2      |1.00  |1.00            |1.00        |       19|** |


Stop reason: b) Subgroup had too few columns. c) There were too few signatures. 

\*\*: 1-PAC > 0.95, \*: 1-PAC > 0.9


### Partition hierarchy

The nodes of the hierarchy can be merged by setting the `merge_node` parameters. Here we 
control the hierarchy with the `min_n_signatures` parameter. The value of `min_n_signatures` is
from `node_info()`.





<style type='text/css'>



.ui-helper-hidden {
	display: none;
}
.ui-helper-hidden-accessible {
	border: 0;
	clip: rect(0 0 0 0);
	height: 1px;
	margin: -1px;
	overflow: hidden;
	padding: 0;
	position: absolute;
	width: 1px;
}
.ui-helper-reset {
	margin: 0;
	padding: 0;
	border: 0;
	outline: 0;
	line-height: 1.3;
	text-decoration: none;
	font-size: 100%;
	list-style: none;
}
.ui-helper-clearfix:before,
.ui-helper-clearfix:after {
	content: "";
	display: table;
	border-collapse: collapse;
}
.ui-helper-clearfix:after {
	clear: both;
}
.ui-helper-zfix {
	width: 100%;
	height: 100%;
	top: 0;
	left: 0;
	position: absolute;
	opacity: 0;
	filter:Alpha(Opacity=0); 
}

.ui-front {
	z-index: 100;
}



.ui-state-disabled {
	cursor: default !important;
	pointer-events: none;
}



.ui-icon {
	display: inline-block;
	vertical-align: middle;
	margin-top: -.25em;
	position: relative;
	text-indent: -99999px;
	overflow: hidden;
	background-repeat: no-repeat;
}

.ui-widget-icon-block {
	left: 50%;
	margin-left: -8px;
	display: block;
}




.ui-widget-overlay {
	position: fixed;
	top: 0;
	left: 0;
	width: 100%;
	height: 100%;
}
.ui-accordion .ui-accordion-header {
	display: block;
	cursor: pointer;
	position: relative;
	margin: 2px 0 0 0;
	padding: .5em .5em .5em .7em;
	font-size: 100%;
}
.ui-accordion .ui-accordion-content {
	padding: 1em 2.2em;
	border-top: 0;
	overflow: auto;
}
.ui-autocomplete {
	position: absolute;
	top: 0;
	left: 0;
	cursor: default;
}
.ui-menu {
	list-style: none;
	padding: 0;
	margin: 0;
	display: block;
	outline: 0;
}
.ui-menu .ui-menu {
	position: absolute;
}
.ui-menu .ui-menu-item {
	margin: 0;
	cursor: pointer;
	
	list-style-image: url("data:image/gif;base64,R0lGODlhAQABAIAAAAAAAP///yH5BAEAAAAALAAAAAABAAEAAAIBRAA7");
}
.ui-menu .ui-menu-item-wrapper {
	position: relative;
	padding: 3px 1em 3px .4em;
}
.ui-menu .ui-menu-divider {
	margin: 5px 0;
	height: 0;
	font-size: 0;
	line-height: 0;
	border-width: 1px 0 0 0;
}
.ui-menu .ui-state-focus,
.ui-menu .ui-state-active {
	margin: -1px;
}


.ui-menu-icons {
	position: relative;
}
.ui-menu-icons .ui-menu-item-wrapper {
	padding-left: 2em;
}


.ui-menu .ui-icon {
	position: absolute;
	top: 0;
	bottom: 0;
	left: .2em;
	margin: auto 0;
}


.ui-menu .ui-menu-icon {
	left: auto;
	right: 0;
}
.ui-button {
	padding: .4em 1em;
	display: inline-block;
	position: relative;
	line-height: normal;
	margin-right: .1em;
	cursor: pointer;
	vertical-align: middle;
	text-align: center;
	-webkit-user-select: none;
	-moz-user-select: none;
	-ms-user-select: none;
	user-select: none;

	
	overflow: visible;
}

.ui-button,
.ui-button:link,
.ui-button:visited,
.ui-button:hover,
.ui-button:active {
	text-decoration: none;
}


.ui-button-icon-only {
	width: 2em;
	box-sizing: border-box;
	text-indent: -9999px;
	white-space: nowrap;
}


input.ui-button.ui-button-icon-only {
	text-indent: 0;
}


.ui-button-icon-only .ui-icon {
	position: absolute;
	top: 50%;
	left: 50%;
	margin-top: -8px;
	margin-left: -8px;
}

.ui-button.ui-icon-notext .ui-icon {
	padding: 0;
	width: 2.1em;
	height: 2.1em;
	text-indent: -9999px;
	white-space: nowrap;

}

input.ui-button.ui-icon-notext .ui-icon {
	width: auto;
	height: auto;
	text-indent: 0;
	white-space: normal;
	padding: .4em 1em;
}



input.ui-button::-moz-focus-inner,
button.ui-button::-moz-focus-inner {
	border: 0;
	padding: 0;
}
.ui-controlgroup {
	vertical-align: middle;
	display: inline-block;
}
.ui-controlgroup > .ui-controlgroup-item {
	float: left;
	margin-left: 0;
	margin-right: 0;
}
.ui-controlgroup > .ui-controlgroup-item:focus,
.ui-controlgroup > .ui-controlgroup-item.ui-visual-focus {
	z-index: 9999;
}
.ui-controlgroup-vertical > .ui-controlgroup-item {
	display: block;
	float: none;
	width: 100%;
	margin-top: 0;
	margin-bottom: 0;
	text-align: left;
}
.ui-controlgroup-vertical .ui-controlgroup-item {
	box-sizing: border-box;
}
.ui-controlgroup .ui-controlgroup-label {
	padding: .4em 1em;
}
.ui-controlgroup .ui-controlgroup-label span {
	font-size: 80%;
}
.ui-controlgroup-horizontal .ui-controlgroup-label + .ui-controlgroup-item {
	border-left: none;
}
.ui-controlgroup-vertical .ui-controlgroup-label + .ui-controlgroup-item {
	border-top: none;
}
.ui-controlgroup-horizontal .ui-controlgroup-label.ui-widget-content {
	border-right: none;
}
.ui-controlgroup-vertical .ui-controlgroup-label.ui-widget-content {
	border-bottom: none;
}


.ui-controlgroup-vertical .ui-spinner-input {

	
	width: 75%;
	width: calc( 100% - 2.4em );
}
.ui-controlgroup-vertical .ui-spinner .ui-spinner-up {
	border-top-style: solid;
}

.ui-checkboxradio-label .ui-icon-background {
	box-shadow: inset 1px 1px 1px #ccc;
	border-radius: .12em;
	border: none;
}
.ui-checkboxradio-radio-label .ui-icon-background {
	width: 16px;
	height: 16px;
	border-radius: 1em;
	overflow: visible;
	border: none;
}
.ui-checkboxradio-radio-label.ui-checkboxradio-checked .ui-icon,
.ui-checkboxradio-radio-label.ui-checkboxradio-checked:hover .ui-icon {
	background-image: none;
	width: 8px;
	height: 8px;
	border-width: 4px;
	border-style: solid;
}
.ui-checkboxradio-disabled {
	pointer-events: none;
}
.ui-datepicker {
	width: 17em;
	padding: .2em .2em 0;
	display: none;
}
.ui-datepicker .ui-datepicker-header {
	position: relative;
	padding: .2em 0;
}
.ui-datepicker .ui-datepicker-prev,
.ui-datepicker .ui-datepicker-next {
	position: absolute;
	top: 2px;
	width: 1.8em;
	height: 1.8em;
}
.ui-datepicker .ui-datepicker-prev-hover,
.ui-datepicker .ui-datepicker-next-hover {
	top: 1px;
}
.ui-datepicker .ui-datepicker-prev {
	left: 2px;
}
.ui-datepicker .ui-datepicker-next {
	right: 2px;
}
.ui-datepicker .ui-datepicker-prev-hover {
	left: 1px;
}
.ui-datepicker .ui-datepicker-next-hover {
	right: 1px;
}
.ui-datepicker .ui-datepicker-prev span,
.ui-datepicker .ui-datepicker-next span {
	display: block;
	position: absolute;
	left: 50%;
	margin-left: -8px;
	top: 50%;
	margin-top: -8px;
}
.ui-datepicker .ui-datepicker-title {
	margin: 0 2.3em;
	line-height: 1.8em;
	text-align: center;
}
.ui-datepicker .ui-datepicker-title select {
	font-size: 1em;
	margin: 1px 0;
}
.ui-datepicker select.ui-datepicker-month,
.ui-datepicker select.ui-datepicker-year {
	width: 45%;
}
.ui-datepicker table {
	width: 100%;
	font-size: .9em;
	border-collapse: collapse;
	margin: 0 0 .4em;
}
.ui-datepicker th {
	padding: .7em .3em;
	text-align: center;
	font-weight: bold;
	border: 0;
}
.ui-datepicker td {
	border: 0;
	padding: 1px;
}
.ui-datepicker td span,
.ui-datepicker td a {
	display: block;
	padding: .2em;
	text-align: right;
	text-decoration: none;
}
.ui-datepicker .ui-datepicker-buttonpane {
	background-image: none;
	margin: .7em 0 0 0;
	padding: 0 .2em;
	border-left: 0;
	border-right: 0;
	border-bottom: 0;
}
.ui-datepicker .ui-datepicker-buttonpane button {
	float: right;
	margin: .5em .2em .4em;
	cursor: pointer;
	padding: .2em .6em .3em .6em;
	width: auto;
	overflow: visible;
}
.ui-datepicker .ui-datepicker-buttonpane button.ui-datepicker-current {
	float: left;
}


.ui-datepicker.ui-datepicker-multi {
	width: auto;
}
.ui-datepicker-multi .ui-datepicker-group {
	float: left;
}
.ui-datepicker-multi .ui-datepicker-group table {
	width: 95%;
	margin: 0 auto .4em;
}
.ui-datepicker-multi-2 .ui-datepicker-group {
	width: 50%;
}
.ui-datepicker-multi-3 .ui-datepicker-group {
	width: 33.3%;
}
.ui-datepicker-multi-4 .ui-datepicker-group {
	width: 25%;
}
.ui-datepicker-multi .ui-datepicker-group-last .ui-datepicker-header,
.ui-datepicker-multi .ui-datepicker-group-middle .ui-datepicker-header {
	border-left-width: 0;
}
.ui-datepicker-multi .ui-datepicker-buttonpane {
	clear: left;
}
.ui-datepicker-row-break {
	clear: both;
	width: 100%;
	font-size: 0;
}


.ui-datepicker-rtl {
	direction: rtl;
}
.ui-datepicker-rtl .ui-datepicker-prev {
	right: 2px;
	left: auto;
}
.ui-datepicker-rtl .ui-datepicker-next {
	left: 2px;
	right: auto;
}
.ui-datepicker-rtl .ui-datepicker-prev:hover {
	right: 1px;
	left: auto;
}
.ui-datepicker-rtl .ui-datepicker-next:hover {
	left: 1px;
	right: auto;
}
.ui-datepicker-rtl .ui-datepicker-buttonpane {
	clear: right;
}
.ui-datepicker-rtl .ui-datepicker-buttonpane button {
	float: left;
}
.ui-datepicker-rtl .ui-datepicker-buttonpane button.ui-datepicker-current,
.ui-datepicker-rtl .ui-datepicker-group {
	float: right;
}
.ui-datepicker-rtl .ui-datepicker-group-last .ui-datepicker-header,
.ui-datepicker-rtl .ui-datepicker-group-middle .ui-datepicker-header {
	border-right-width: 0;
	border-left-width: 1px;
}


.ui-datepicker .ui-icon {
	display: block;
	text-indent: -99999px;
	overflow: hidden;
	background-repeat: no-repeat;
	left: .5em;
	top: .3em;
}
.ui-dialog {
	position: absolute;
	top: 0;
	left: 0;
	padding: .2em;
	outline: 0;
}
.ui-dialog .ui-dialog-titlebar {
	padding: .4em 1em;
	position: relative;
}
.ui-dialog .ui-dialog-title {
	float: left;
	margin: .1em 0;
	white-space: nowrap;
	width: 90%;
	overflow: hidden;
	text-overflow: ellipsis;
}
.ui-dialog .ui-dialog-titlebar-close {
	position: absolute;
	right: .3em;
	top: 50%;
	width: 20px;
	margin: -10px 0 0 0;
	padding: 1px;
	height: 20px;
}
.ui-dialog .ui-dialog-content {
	position: relative;
	border: 0;
	padding: .5em 1em;
	background: none;
	overflow: auto;
}
.ui-dialog .ui-dialog-buttonpane {
	text-align: left;
	border-width: 1px 0 0 0;
	background-image: none;
	margin-top: .5em;
	padding: .3em 1em .5em .4em;
}
.ui-dialog .ui-dialog-buttonpane .ui-dialog-buttonset {
	float: right;
}
.ui-dialog .ui-dialog-buttonpane button {
	margin: .5em .4em .5em 0;
	cursor: pointer;
}
.ui-dialog .ui-resizable-n {
	height: 2px;
	top: 0;
}
.ui-dialog .ui-resizable-e {
	width: 2px;
	right: 0;
}
.ui-dialog .ui-resizable-s {
	height: 2px;
	bottom: 0;
}
.ui-dialog .ui-resizable-w {
	width: 2px;
	left: 0;
}
.ui-dialog .ui-resizable-se,
.ui-dialog .ui-resizable-sw,
.ui-dialog .ui-resizable-ne,
.ui-dialog .ui-resizable-nw {
	width: 7px;
	height: 7px;
}
.ui-dialog .ui-resizable-se {
	right: 0;
	bottom: 0;
}
.ui-dialog .ui-resizable-sw {
	left: 0;
	bottom: 0;
}
.ui-dialog .ui-resizable-ne {
	right: 0;
	top: 0;
}
.ui-dialog .ui-resizable-nw {
	left: 0;
	top: 0;
}
.ui-draggable .ui-dialog-titlebar {
	cursor: move;
}
.ui-draggable-handle {
	-ms-touch-action: none;
	touch-action: none;
}
.ui-resizable {
	position: relative;
}
.ui-resizable-handle {
	position: absolute;
	font-size: 0.1px;
	display: block;
	-ms-touch-action: none;
	touch-action: none;
}
.ui-resizable-disabled .ui-resizable-handle,
.ui-resizable-autohide .ui-resizable-handle {
	display: none;
}
.ui-resizable-n {
	cursor: n-resize;
	height: 7px;
	width: 100%;
	top: -5px;
	left: 0;
}
.ui-resizable-s {
	cursor: s-resize;
	height: 7px;
	width: 100%;
	bottom: -5px;
	left: 0;
}
.ui-resizable-e {
	cursor: e-resize;
	width: 7px;
	right: -5px;
	top: 0;
	height: 100%;
}
.ui-resizable-w {
	cursor: w-resize;
	width: 7px;
	left: -5px;
	top: 0;
	height: 100%;
}
.ui-resizable-se {
	cursor: se-resize;
	width: 12px;
	height: 12px;
	right: 1px;
	bottom: 1px;
}
.ui-resizable-sw {
	cursor: sw-resize;
	width: 9px;
	height: 9px;
	left: -5px;
	bottom: -5px;
}
.ui-resizable-nw {
	cursor: nw-resize;
	width: 9px;
	height: 9px;
	left: -5px;
	top: -5px;
}
.ui-resizable-ne {
	cursor: ne-resize;
	width: 9px;
	height: 9px;
	right: -5px;
	top: -5px;
}
.ui-progressbar {
	height: 2em;
	text-align: left;
	overflow: hidden;
}
.ui-progressbar .ui-progressbar-value {
	margin: -1px;
	height: 100%;
}
.ui-progressbar .ui-progressbar-overlay {
	background: url("data:image/gif;base64,R0lGODlhKAAoAIABAAAAAP///yH/C05FVFNDQVBFMi4wAwEAAAAh+QQJAQABACwAAAAAKAAoAAACkYwNqXrdC52DS06a7MFZI+4FHBCKoDeWKXqymPqGqxvJrXZbMx7Ttc+w9XgU2FB3lOyQRWET2IFGiU9m1frDVpxZZc6bfHwv4c1YXP6k1Vdy292Fb6UkuvFtXpvWSzA+HycXJHUXiGYIiMg2R6W459gnWGfHNdjIqDWVqemH2ekpObkpOlppWUqZiqr6edqqWQAAIfkECQEAAQAsAAAAACgAKAAAApSMgZnGfaqcg1E2uuzDmmHUBR8Qil95hiPKqWn3aqtLsS18y7G1SzNeowWBENtQd+T1JktP05nzPTdJZlR6vUxNWWjV+vUWhWNkWFwxl9VpZRedYcflIOLafaa28XdsH/ynlcc1uPVDZxQIR0K25+cICCmoqCe5mGhZOfeYSUh5yJcJyrkZWWpaR8doJ2o4NYq62lAAACH5BAkBAAEALAAAAAAoACgAAAKVDI4Yy22ZnINRNqosw0Bv7i1gyHUkFj7oSaWlu3ovC8GxNso5fluz3qLVhBVeT/Lz7ZTHyxL5dDalQWPVOsQWtRnuwXaFTj9jVVh8pma9JjZ4zYSj5ZOyma7uuolffh+IR5aW97cHuBUXKGKXlKjn+DiHWMcYJah4N0lYCMlJOXipGRr5qdgoSTrqWSq6WFl2ypoaUAAAIfkECQEAAQAsAAAAACgAKAAAApaEb6HLgd/iO7FNWtcFWe+ufODGjRfoiJ2akShbueb0wtI50zm02pbvwfWEMWBQ1zKGlLIhskiEPm9R6vRXxV4ZzWT2yHOGpWMyorblKlNp8HmHEb/lCXjcW7bmtXP8Xt229OVWR1fod2eWqNfHuMjXCPkIGNileOiImVmCOEmoSfn3yXlJWmoHGhqp6ilYuWYpmTqKUgAAIfkECQEAAQAsAAAAACgAKAAAApiEH6kb58biQ3FNWtMFWW3eNVcojuFGfqnZqSebuS06w5V80/X02pKe8zFwP6EFWOT1lDFk8rGERh1TTNOocQ61Hm4Xm2VexUHpzjymViHrFbiELsefVrn6XKfnt2Q9G/+Xdie499XHd2g4h7ioOGhXGJboGAnXSBnoBwKYyfioubZJ2Hn0RuRZaflZOil56Zp6iioKSXpUAAAh+QQJAQABACwAAAAAKAAoAAACkoQRqRvnxuI7kU1a1UU5bd5tnSeOZXhmn5lWK3qNTWvRdQxP8qvaC+/yaYQzXO7BMvaUEmJRd3TsiMAgswmNYrSgZdYrTX6tSHGZO73ezuAw2uxuQ+BbeZfMxsexY35+/Qe4J1inV0g4x3WHuMhIl2jXOKT2Q+VU5fgoSUI52VfZyfkJGkha6jmY+aaYdirq+lQAACH5BAkBAAEALAAAAAAoACgAAAKWBIKpYe0L3YNKToqswUlvznigd4wiR4KhZrKt9Upqip61i9E3vMvxRdHlbEFiEXfk9YARYxOZZD6VQ2pUunBmtRXo1Lf8hMVVcNl8JafV38aM2/Fu5V16Bn63r6xt97j09+MXSFi4BniGFae3hzbH9+hYBzkpuUh5aZmHuanZOZgIuvbGiNeomCnaxxap2upaCZsq+1kAACH5BAkBAAEALAAAAAAoACgAAAKXjI8By5zf4kOxTVrXNVlv1X0d8IGZGKLnNpYtm8Lr9cqVeuOSvfOW79D9aDHizNhDJidFZhNydEahOaDH6nomtJjp1tutKoNWkvA6JqfRVLHU/QUfau9l2x7G54d1fl995xcIGAdXqMfBNadoYrhH+Mg2KBlpVpbluCiXmMnZ2Sh4GBqJ+ckIOqqJ6LmKSllZmsoq6wpQAAAh+QQJAQABACwAAAAAKAAoAAAClYx/oLvoxuJDkU1a1YUZbJ59nSd2ZXhWqbRa2/gF8Gu2DY3iqs7yrq+xBYEkYvFSM8aSSObE+ZgRl1BHFZNr7pRCavZ5BW2142hY3AN/zWtsmf12p9XxxFl2lpLn1rseztfXZjdIWIf2s5dItwjYKBgo9yg5pHgzJXTEeGlZuenpyPmpGQoKOWkYmSpaSnqKileI2FAAACH5BAkBAAEALAAAAAAoACgAAAKVjB+gu+jG4kORTVrVhRlsnn2dJ3ZleFaptFrb+CXmO9OozeL5VfP99HvAWhpiUdcwkpBH3825AwYdU8xTqlLGhtCosArKMpvfa1mMRae9VvWZfeB2XfPkeLmm18lUcBj+p5dnN8jXZ3YIGEhYuOUn45aoCDkp16hl5IjYJvjWKcnoGQpqyPlpOhr3aElaqrq56Bq7VAAAOw==");
	height: 100%;
	filter: alpha(opacity=25); 
	opacity: 0.25;
}
.ui-progressbar-indeterminate .ui-progressbar-value {
	background-image: none;
}
.ui-selectable {
	-ms-touch-action: none;
	touch-action: none;
}
.ui-selectable-helper {
	position: absolute;
	z-index: 100;
	border: 1px dotted black;
}
.ui-selectmenu-menu {
	padding: 0;
	margin: 0;
	position: absolute;
	top: 0;
	left: 0;
	display: none;
}
.ui-selectmenu-menu .ui-menu {
	overflow: auto;
	overflow-x: hidden;
	padding-bottom: 1px;
}
.ui-selectmenu-menu .ui-menu .ui-selectmenu-optgroup {
	font-size: 1em;
	font-weight: bold;
	line-height: 1.5;
	padding: 2px 0.4em;
	margin: 0.5em 0 0 0;
	height: auto;
	border: 0;
}
.ui-selectmenu-open {
	display: block;
}
.ui-selectmenu-text {
	display: block;
	margin-right: 20px;
	overflow: hidden;
	text-overflow: ellipsis;
}
.ui-selectmenu-button.ui-button {
	text-align: left;
	white-space: nowrap;
	width: 14em;
}
.ui-selectmenu-icon.ui-icon {
	float: right;
	margin-top: 0;
}
.ui-slider {
	position: relative;
	text-align: left;
}
.ui-slider .ui-slider-handle {
	position: absolute;
	z-index: 2;
	width: 1.2em;
	height: 1.2em;
	cursor: default;
	-ms-touch-action: none;
	touch-action: none;
}
.ui-slider .ui-slider-range {
	position: absolute;
	z-index: 1;
	font-size: .7em;
	display: block;
	border: 0;
	background-position: 0 0;
}


.ui-slider.ui-state-disabled .ui-slider-handle,
.ui-slider.ui-state-disabled .ui-slider-range {
	filter: inherit;
}

.ui-slider-horizontal {
	height: .8em;
}
.ui-slider-horizontal .ui-slider-handle {
	top: -.3em;
	margin-left: -.6em;
}
.ui-slider-horizontal .ui-slider-range {
	top: 0;
	height: 100%;
}
.ui-slider-horizontal .ui-slider-range-min {
	left: 0;
}
.ui-slider-horizontal .ui-slider-range-max {
	right: 0;
}

.ui-slider-vertical {
	width: .8em;
	height: 100px;
}
.ui-slider-vertical .ui-slider-handle {
	left: -.3em;
	margin-left: 0;
	margin-bottom: -.6em;
}
.ui-slider-vertical .ui-slider-range {
	left: 0;
	width: 100%;
}
.ui-slider-vertical .ui-slider-range-min {
	bottom: 0;
}
.ui-slider-vertical .ui-slider-range-max {
	top: 0;
}
.ui-sortable-handle {
	-ms-touch-action: none;
	touch-action: none;
}
.ui-spinner {
	position: relative;
	display: inline-block;
	overflow: hidden;
	padding: 0;
	vertical-align: middle;
}
.ui-spinner-input {
	border: none;
	background: none;
	color: inherit;
	padding: .222em 0;
	margin: .2em 0;
	vertical-align: middle;
	margin-left: .4em;
	margin-right: 2em;
}
.ui-spinner-button {
	width: 1.6em;
	height: 50%;
	font-size: .5em;
	padding: 0;
	margin: 0;
	text-align: center;
	position: absolute;
	cursor: default;
	display: block;
	overflow: hidden;
	right: 0;
}

.ui-spinner a.ui-spinner-button {
	border-top-style: none;
	border-bottom-style: none;
	border-right-style: none;
}
.ui-spinner-up {
	top: 0;
}
.ui-spinner-down {
	bottom: 0;
}
.ui-tabs {
	position: relative;
	padding: .2em;
}
.ui-tabs .ui-tabs-nav {
	margin: 0;
	padding: .2em .2em 0;
}
.ui-tabs .ui-tabs-nav li {
	list-style: none;
	float: left;
	position: relative;
	top: 0;
	margin: 1px .2em 0 0;
	border-bottom-width: 0;
	padding: 0;
	white-space: nowrap;
}
.ui-tabs .ui-tabs-nav .ui-tabs-anchor {
	float: left;
	padding: .5em 1em;
	text-decoration: none;
}
.ui-tabs .ui-tabs-nav li.ui-tabs-active {
	margin-bottom: -1px;
	padding-bottom: 1px;
}
.ui-tabs .ui-tabs-nav li.ui-tabs-active .ui-tabs-anchor,
.ui-tabs .ui-tabs-nav li.ui-state-disabled .ui-tabs-anchor,
.ui-tabs .ui-tabs-nav li.ui-tabs-loading .ui-tabs-anchor {
	cursor: text;
}
.ui-tabs-collapsible .ui-tabs-nav li.ui-tabs-active .ui-tabs-anchor {
	cursor: pointer;
}
.ui-tabs .ui-tabs-panel {
	display: block;
	border-width: 0;
	padding: 1em 1.4em;
	background: none;
}
.ui-tooltip {
	padding: 8px;
	position: absolute;
	z-index: 9999;
	max-width: 300px;
}
body .ui-tooltip {
	border-width: 2px;
}

.ui-widget {
	font-family: Arial,Helvetica,sans-serif;
	font-size: 1em;
}
.ui-widget .ui-widget {
	font-size: 1em;
}
.ui-widget input,
.ui-widget select,
.ui-widget textarea,
.ui-widget button {
	font-family: Arial,Helvetica,sans-serif;
	font-size: 1em;
}
.ui-widget.ui-widget-content {
	border: 1px solid #c5c5c5;
}
.ui-widget-content {
	border: 1px solid #dddddd;
	background: #ffffff;
	color: #333333;
}
.ui-widget-content a {
	color: #333333;
}
.ui-widget-header {
	border: 1px solid #dddddd;
	background: #e9e9e9;
	color: #333333;
	font-weight: bold;
}
.ui-widget-header a {
	color: #333333;
}


.ui-state-default,
.ui-widget-content .ui-state-default,
.ui-widget-header .ui-state-default,
.ui-button,


html .ui-button.ui-state-disabled:hover,
html .ui-button.ui-state-disabled:active {
	border: 1px solid #c5c5c5;
	background: #f6f6f6;
	font-weight: normal;
	color: #454545;
}
.ui-state-default a,
.ui-state-default a:link,
.ui-state-default a:visited,
a.ui-button,
a:link.ui-button,
a:visited.ui-button,
.ui-button {
	color: #454545;
	text-decoration: none;
}
.ui-state-hover,
.ui-widget-content .ui-state-hover,
.ui-widget-header .ui-state-hover,
.ui-state-focus,
.ui-widget-content .ui-state-focus,
.ui-widget-header .ui-state-focus,
.ui-button:hover,
.ui-button:focus {
	border: 1px solid #cccccc;
	background: #ededed;
	font-weight: normal;
	color: #2b2b2b;
}
.ui-state-hover a,
.ui-state-hover a:hover,
.ui-state-hover a:link,
.ui-state-hover a:visited,
.ui-state-focus a,
.ui-state-focus a:hover,
.ui-state-focus a:link,
.ui-state-focus a:visited,
a.ui-button:hover,
a.ui-button:focus {
	color: #2b2b2b;
	text-decoration: none;
}

.ui-visual-focus {
	box-shadow: 0 0 3px 1px rgb(94, 158, 214);
}
.ui-state-active,
.ui-widget-content .ui-state-active,
.ui-widget-header .ui-state-active,
a.ui-button:active,
.ui-button:active,
.ui-button.ui-state-active:hover {
	border: 1px solid #003eff;
	background: #007fff;
	font-weight: normal;
	color: #ffffff;
}
.ui-icon-background,
.ui-state-active .ui-icon-background {
	border: #003eff;
	background-color: #ffffff;
}
.ui-state-active a,
.ui-state-active a:link,
.ui-state-active a:visited {
	color: #ffffff;
	text-decoration: none;
}


.ui-state-highlight,
.ui-widget-content .ui-state-highlight,
.ui-widget-header .ui-state-highlight {
	border: 1px solid #dad55e;
	background: #fffa90;
	color: #777620;
}
.ui-state-checked {
	border: 1px solid #dad55e;
	background: #fffa90;
}
.ui-state-highlight a,
.ui-widget-content .ui-state-highlight a,
.ui-widget-header .ui-state-highlight a {
	color: #777620;
}
.ui-state-error,
.ui-widget-content .ui-state-error,
.ui-widget-header .ui-state-error {
	border: 1px solid #f1a899;
	background: #fddfdf;
	color: #5f3f3f;
}
.ui-state-error a,
.ui-widget-content .ui-state-error a,
.ui-widget-header .ui-state-error a {
	color: #5f3f3f;
}
.ui-state-error-text,
.ui-widget-content .ui-state-error-text,
.ui-widget-header .ui-state-error-text {
	color: #5f3f3f;
}
.ui-priority-primary,
.ui-widget-content .ui-priority-primary,
.ui-widget-header .ui-priority-primary {
	font-weight: bold;
}
.ui-priority-secondary,
.ui-widget-content .ui-priority-secondary,
.ui-widget-header .ui-priority-secondary {
	opacity: .7;
	filter:Alpha(Opacity=70); 
	font-weight: normal;
}
.ui-state-disabled,
.ui-widget-content .ui-state-disabled,
.ui-widget-header .ui-state-disabled {
	opacity: .35;
	filter:Alpha(Opacity=35); 
	background-image: none;
}
.ui-state-disabled .ui-icon {
	filter:Alpha(Opacity=35); 
}




.ui-icon {
	width: 16px;
	height: 16px;
}
.ui-icon,
.ui-widget-content .ui-icon {
	background-image: url("images/ui-icons_444444_256x240.png");
}
.ui-widget-header .ui-icon {
	background-image: url("images/ui-icons_444444_256x240.png");
}
.ui-state-hover .ui-icon,
.ui-state-focus .ui-icon,
.ui-button:hover .ui-icon,
.ui-button:focus .ui-icon {
	background-image: url("images/ui-icons_555555_256x240.png");
}
.ui-state-active .ui-icon,
.ui-button:active .ui-icon {
	background-image: url("images/ui-icons_ffffff_256x240.png");
}
.ui-state-highlight .ui-icon,
.ui-button .ui-state-highlight.ui-icon {
	background-image: url("images/ui-icons_777620_256x240.png");
}
.ui-state-error .ui-icon,
.ui-state-error-text .ui-icon {
	background-image: url("images/ui-icons_cc0000_256x240.png");
}
.ui-button .ui-icon {
	background-image: url("images/ui-icons_777777_256x240.png");
}


.ui-icon-blank { background-position: 16px 16px; }
.ui-icon-caret-1-n { background-position: 0 0; }
.ui-icon-caret-1-ne { background-position: -16px 0; }
.ui-icon-caret-1-e { background-position: -32px 0; }
.ui-icon-caret-1-se { background-position: -48px 0; }
.ui-icon-caret-1-s { background-position: -65px 0; }
.ui-icon-caret-1-sw { background-position: -80px 0; }
.ui-icon-caret-1-w { background-position: -96px 0; }
.ui-icon-caret-1-nw { background-position: -112px 0; }
.ui-icon-caret-2-n-s { background-position: -128px 0; }
.ui-icon-caret-2-e-w { background-position: -144px 0; }
.ui-icon-triangle-1-n { background-position: 0 -16px; }
.ui-icon-triangle-1-ne { background-position: -16px -16px; }
.ui-icon-triangle-1-e { background-position: -32px -16px; }
.ui-icon-triangle-1-se { background-position: -48px -16px; }
.ui-icon-triangle-1-s { background-position: -65px -16px; }
.ui-icon-triangle-1-sw { background-position: -80px -16px; }
.ui-icon-triangle-1-w { background-position: -96px -16px; }
.ui-icon-triangle-1-nw { background-position: -112px -16px; }
.ui-icon-triangle-2-n-s { background-position: -128px -16px; }
.ui-icon-triangle-2-e-w { background-position: -144px -16px; }
.ui-icon-arrow-1-n { background-position: 0 -32px; }
.ui-icon-arrow-1-ne { background-position: -16px -32px; }
.ui-icon-arrow-1-e { background-position: -32px -32px; }
.ui-icon-arrow-1-se { background-position: -48px -32px; }
.ui-icon-arrow-1-s { background-position: -65px -32px; }
.ui-icon-arrow-1-sw { background-position: -80px -32px; }
.ui-icon-arrow-1-w { background-position: -96px -32px; }
.ui-icon-arrow-1-nw { background-position: -112px -32px; }
.ui-icon-arrow-2-n-s { background-position: -128px -32px; }
.ui-icon-arrow-2-ne-sw { background-position: -144px -32px; }
.ui-icon-arrow-2-e-w { background-position: -160px -32px; }
.ui-icon-arrow-2-se-nw { background-position: -176px -32px; }
.ui-icon-arrowstop-1-n { background-position: -192px -32px; }
.ui-icon-arrowstop-1-e { background-position: -208px -32px; }
.ui-icon-arrowstop-1-s { background-position: -224px -32px; }
.ui-icon-arrowstop-1-w { background-position: -240px -32px; }
.ui-icon-arrowthick-1-n { background-position: 1px -48px; }
.ui-icon-arrowthick-1-ne { background-position: -16px -48px; }
.ui-icon-arrowthick-1-e { background-position: -32px -48px; }
.ui-icon-arrowthick-1-se { background-position: -48px -48px; }
.ui-icon-arrowthick-1-s { background-position: -64px -48px; }
.ui-icon-arrowthick-1-sw { background-position: -80px -48px; }
.ui-icon-arrowthick-1-w { background-position: -96px -48px; }
.ui-icon-arrowthick-1-nw { background-position: -112px -48px; }
.ui-icon-arrowthick-2-n-s { background-position: -128px -48px; }
.ui-icon-arrowthick-2-ne-sw { background-position: -144px -48px; }
.ui-icon-arrowthick-2-e-w { background-position: -160px -48px; }
.ui-icon-arrowthick-2-se-nw { background-position: -176px -48px; }
.ui-icon-arrowthickstop-1-n { background-position: -192px -48px; }
.ui-icon-arrowthickstop-1-e { background-position: -208px -48px; }
.ui-icon-arrowthickstop-1-s { background-position: -224px -48px; }
.ui-icon-arrowthickstop-1-w { background-position: -240px -48px; }
.ui-icon-arrowreturnthick-1-w { background-position: 0 -64px; }
.ui-icon-arrowreturnthick-1-n { background-position: -16px -64px; }
.ui-icon-arrowreturnthick-1-e { background-position: -32px -64px; }
.ui-icon-arrowreturnthick-1-s { background-position: -48px -64px; }
.ui-icon-arrowreturn-1-w { background-position: -64px -64px; }
.ui-icon-arrowreturn-1-n { background-position: -80px -64px; }
.ui-icon-arrowreturn-1-e { background-position: -96px -64px; }
.ui-icon-arrowreturn-1-s { background-position: -112px -64px; }
.ui-icon-arrowrefresh-1-w { background-position: -128px -64px; }
.ui-icon-arrowrefresh-1-n { background-position: -144px -64px; }
.ui-icon-arrowrefresh-1-e { background-position: -160px -64px; }
.ui-icon-arrowrefresh-1-s { background-position: -176px -64px; }
.ui-icon-arrow-4 { background-position: 0 -80px; }
.ui-icon-arrow-4-diag { background-position: -16px -80px; }
.ui-icon-extlink { background-position: -32px -80px; }
.ui-icon-newwin { background-position: -48px -80px; }
.ui-icon-refresh { background-position: -64px -80px; }
.ui-icon-shuffle { background-position: -80px -80px; }
.ui-icon-transfer-e-w { background-position: -96px -80px; }
.ui-icon-transferthick-e-w { background-position: -112px -80px; }
.ui-icon-folder-collapsed { background-position: 0 -96px; }
.ui-icon-folder-open { background-position: -16px -96px; }
.ui-icon-document { background-position: -32px -96px; }
.ui-icon-document-b { background-position: -48px -96px; }
.ui-icon-note { background-position: -64px -96px; }
.ui-icon-mail-closed { background-position: -80px -96px; }
.ui-icon-mail-open { background-position: -96px -96px; }
.ui-icon-suitcase { background-position: -112px -96px; }
.ui-icon-comment { background-position: -128px -96px; }
.ui-icon-person { background-position: -144px -96px; }
.ui-icon-print { background-position: -160px -96px; }
.ui-icon-trash { background-position: -176px -96px; }
.ui-icon-locked { background-position: -192px -96px; }
.ui-icon-unlocked { background-position: -208px -96px; }
.ui-icon-bookmark { background-position: -224px -96px; }
.ui-icon-tag { background-position: -240px -96px; }
.ui-icon-home { background-position: 0 -112px; }
.ui-icon-flag { background-position: -16px -112px; }
.ui-icon-calendar { background-position: -32px -112px; }
.ui-icon-cart { background-position: -48px -112px; }
.ui-icon-pencil { background-position: -64px -112px; }
.ui-icon-clock { background-position: -80px -112px; }
.ui-icon-disk { background-position: -96px -112px; }
.ui-icon-calculator { background-position: -112px -112px; }
.ui-icon-zoomin { background-position: -128px -112px; }
.ui-icon-zoomout { background-position: -144px -112px; }
.ui-icon-search { background-position: -160px -112px; }
.ui-icon-wrench { background-position: -176px -112px; }
.ui-icon-gear { background-position: -192px -112px; }
.ui-icon-heart { background-position: -208px -112px; }
.ui-icon-star { background-position: -224px -112px; }
.ui-icon-link { background-position: -240px -112px; }
.ui-icon-cancel { background-position: 0 -128px; }
.ui-icon-plus { background-position: -16px -128px; }
.ui-icon-plusthick { background-position: -32px -128px; }
.ui-icon-minus { background-position: -48px -128px; }
.ui-icon-minusthick { background-position: -64px -128px; }
.ui-icon-close { background-position: -80px -128px; }
.ui-icon-closethick { background-position: -96px -128px; }
.ui-icon-key { background-position: -112px -128px; }
.ui-icon-lightbulb { background-position: -128px -128px; }
.ui-icon-scissors { background-position: -144px -128px; }
.ui-icon-clipboard { background-position: -160px -128px; }
.ui-icon-copy { background-position: -176px -128px; }
.ui-icon-contact { background-position: -192px -128px; }
.ui-icon-image { background-position: -208px -128px; }
.ui-icon-video { background-position: -224px -128px; }
.ui-icon-script { background-position: -240px -128px; }
.ui-icon-alert { background-position: 0 -144px; }
.ui-icon-info { background-position: -16px -144px; }
.ui-icon-notice { background-position: -32px -144px; }
.ui-icon-help { background-position: -48px -144px; }
.ui-icon-check { background-position: -64px -144px; }
.ui-icon-bullet { background-position: -80px -144px; }
.ui-icon-radio-on { background-position: -96px -144px; }
.ui-icon-radio-off { background-position: -112px -144px; }
.ui-icon-pin-w { background-position: -128px -144px; }
.ui-icon-pin-s { background-position: -144px -144px; }
.ui-icon-play { background-position: 0 -160px; }
.ui-icon-pause { background-position: -16px -160px; }
.ui-icon-seek-next { background-position: -32px -160px; }
.ui-icon-seek-prev { background-position: -48px -160px; }
.ui-icon-seek-end { background-position: -64px -160px; }
.ui-icon-seek-start { background-position: -80px -160px; }

.ui-icon-seek-first { background-position: -80px -160px; }
.ui-icon-stop { background-position: -96px -160px; }
.ui-icon-eject { background-position: -112px -160px; }
.ui-icon-volume-off { background-position: -128px -160px; }
.ui-icon-volume-on { background-position: -144px -160px; }
.ui-icon-power { background-position: 0 -176px; }
.ui-icon-signal-diag { background-position: -16px -176px; }
.ui-icon-signal { background-position: -32px -176px; }
.ui-icon-battery-0 { background-position: -48px -176px; }
.ui-icon-battery-1 { background-position: -64px -176px; }
.ui-icon-battery-2 { background-position: -80px -176px; }
.ui-icon-battery-3 { background-position: -96px -176px; }
.ui-icon-circle-plus { background-position: 0 -192px; }
.ui-icon-circle-minus { background-position: -16px -192px; }
.ui-icon-circle-close { background-position: -32px -192px; }
.ui-icon-circle-triangle-e { background-position: -48px -192px; }
.ui-icon-circle-triangle-s { background-position: -64px -192px; }
.ui-icon-circle-triangle-w { background-position: -80px -192px; }
.ui-icon-circle-triangle-n { background-position: -96px -192px; }
.ui-icon-circle-arrow-e { background-position: -112px -192px; }
.ui-icon-circle-arrow-s { background-position: -128px -192px; }
.ui-icon-circle-arrow-w { background-position: -144px -192px; }
.ui-icon-circle-arrow-n { background-position: -160px -192px; }
.ui-icon-circle-zoomin { background-position: -176px -192px; }
.ui-icon-circle-zoomout { background-position: -192px -192px; }
.ui-icon-circle-check { background-position: -208px -192px; }
.ui-icon-circlesmall-plus { background-position: 0 -208px; }
.ui-icon-circlesmall-minus { background-position: -16px -208px; }
.ui-icon-circlesmall-close { background-position: -32px -208px; }
.ui-icon-squaresmall-plus { background-position: -48px -208px; }
.ui-icon-squaresmall-minus { background-position: -64px -208px; }
.ui-icon-squaresmall-close { background-position: -80px -208px; }
.ui-icon-grip-dotted-vertical { background-position: 0 -224px; }
.ui-icon-grip-dotted-horizontal { background-position: -16px -224px; }
.ui-icon-grip-solid-vertical { background-position: -32px -224px; }
.ui-icon-grip-solid-horizontal { background-position: -48px -224px; }
.ui-icon-gripsmall-diagonal-se { background-position: -64px -224px; }
.ui-icon-grip-diagonal-se { background-position: -80px -224px; }





.ui-corner-all,
.ui-corner-top,
.ui-corner-left,
.ui-corner-tl {
	border-top-left-radius: 3px;
}
.ui-corner-all,
.ui-corner-top,
.ui-corner-right,
.ui-corner-tr {
	border-top-right-radius: 3px;
}
.ui-corner-all,
.ui-corner-bottom,
.ui-corner-left,
.ui-corner-bl {
	border-bottom-left-radius: 3px;
}
.ui-corner-all,
.ui-corner-bottom,
.ui-corner-right,
.ui-corner-br {
	border-bottom-right-radius: 3px;
}


.ui-widget-overlay {
	background: #aaaaaa;
	opacity: .3;
	filter: Alpha(Opacity=30); 
}
.ui-widget-shadow {
	-webkit-box-shadow: 0px 0px 5px #666666;
	box-shadow: 0px 0px 5px #666666;
} 
</style>
<script src='js/jquery-1.12.4.js'></script>
<script src='js/jquery-ui.js'></script>

<script>
$( function() {
	$( '#tabs-collect-classes-from-hierarchical-partition' ).tabs();
} );
</script>
<div id='tabs-collect-classes-from-hierarchical-partition'>
<ul>
<li><a href='#tab-collect-classes-from-hierarchical-partition-1'>n_signatures ≥ 4690</a></li>
<li><a href='#tab-collect-classes-from-hierarchical-partition-2'>n_signatures ≥ 5323</a></li>
<li><a href='#tab-collect-classes-from-hierarchical-partition-3'>n_signatures ≥ 17040</a></li>
</ul>
<div id='tab-collect-classes-from-hierarchical-partition-1'>
<pre><code class="r">collect_classes(res_rh, merge_node = merge_node_param(min_n_signatures = 4690))
</code></pre>

<p><img src="figure_cola/tab-collect-classes-from-hierarchical-partition-1-1.png" alt="plot of chunk tab-collect-classes-from-hierarchical-partition-1"/></p>

</div>
<div id='tab-collect-classes-from-hierarchical-partition-2'>
<pre><code class="r">collect_classes(res_rh, merge_node = merge_node_param(min_n_signatures = 5323))
</code></pre>

<p><img src="figure_cola/tab-collect-classes-from-hierarchical-partition-2-1.png" alt="plot of chunk tab-collect-classes-from-hierarchical-partition-2"/></p>

</div>
<div id='tab-collect-classes-from-hierarchical-partition-3'>
<pre><code class="r">collect_classes(res_rh, merge_node = merge_node_param(min_n_signatures = 17040))
</code></pre>

<p><img src="figure_cola/tab-collect-classes-from-hierarchical-partition-3-1.png" alt="plot of chunk tab-collect-classes-from-hierarchical-partition-3"/></p>

</div>
</div>

Following shows the table of the partitions (You need to click the **show/hide
code output** link to see it).


<script>
$( function() {
	$( '#tabs-get-classes-from-hierarchical-partition' ).tabs();
} );
</script>
<div id='tabs-get-classes-from-hierarchical-partition'>
<ul>
<li><a href='#tab-get-classes-from-hierarchical-partition-1'>n_signatures ≥ 4690</a></li>
<li><a href='#tab-get-classes-from-hierarchical-partition-2'>n_signatures ≥ 5323</a></li>
<li><a href='#tab-get-classes-from-hierarchical-partition-3'>n_signatures ≥ 17040</a></li>
</ul>

<div id='tab-get-classes-from-hierarchical-partition-1'>
<p><a id='tab-get-classes-from-hierarchical-partition-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">get_classes(res_rh, merge_node = merge_node_param(min_n_signatures = 4690))
</code></pre>

<pre><code>#&gt; TCGA.N6.A4VF.01 TCGA.N7.A59B.01 TCGA.N7.A4Y0.01 TCGA.QN.A5NN.01 TCGA.N6.A4VD.01 TCGA.N8.A4PN.01 
#&gt;            &quot;03&quot;           &quot;012&quot;            &quot;03&quot;            &quot;03&quot;            &quot;03&quot;           &quot;022&quot; 
#&gt; TCGA.N5.A59E.01 TCGA.NG.A4VW.01 TCGA.N8.A4PO.01 TCGA.N7.A4Y5.01 TCGA.ND.A4WA.01 TCGA.N6.A4VC.01 
#&gt;           &quot;021&quot;           &quot;012&quot;           &quot;013&quot;           &quot;022&quot;           &quot;011&quot;           &quot;013&quot; 
#&gt; TCGA.N6.A4VG.01 TCGA.N5.A59F.01 TCGA.N5.A4RD.01 TCGA.N9.A4Q3.01 TCGA.ND.A4WC.01 TCGA.N8.A4PI.01 
#&gt;           &quot;013&quot;            &quot;03&quot;           &quot;021&quot;           &quot;011&quot;           &quot;012&quot;            &quot;03&quot; 
#&gt; TCGA.NA.A5I1.01 TCGA.N8.A56S.01 TCGA.N8.A4PL.01 TCGA.NA.A4R0.01 TCGA.N8.A4PQ.01 TCGA.NA.A4QV.01 
#&gt;           &quot;011&quot;           &quot;011&quot;           &quot;011&quot;           &quot;013&quot;           &quot;022&quot;            &quot;03&quot; 
#&gt; TCGA.NF.A5CP.01 TCGA.ND.A4WF.01 TCGA.N5.A4R8.01 TCGA.NA.A4QY.01 TCGA.NA.A4QW.01 TCGA.N9.A4Q1.01 
#&gt;            &quot;03&quot;           &quot;021&quot;           &quot;023&quot;           &quot;013&quot;           &quot;022&quot;           &quot;022&quot; 
#&gt; TCGA.N5.A4RV.01 TCGA.NF.A4WX.01 TCGA.N8.A4PP.01 TCGA.N5.A4RS.01 TCGA.N5.A4RN.01 TCGA.N5.A4RT.01 
#&gt;           &quot;023&quot;           &quot;011&quot;           &quot;012&quot;           &quot;011&quot;            &quot;03&quot;            &quot;03&quot; 
#&gt; TCGA.N5.A4RF.01 TCGA.QM.A5NM.01 TCGA.N9.A4PZ.01 TCGA.N5.A4RO.01 TCGA.NA.A4QX.01 TCGA.N5.A4RA.01 
#&gt;            &quot;03&quot;            &quot;03&quot;            &quot;03&quot;           &quot;021&quot;           &quot;022&quot;           &quot;022&quot; 
#&gt; TCGA.N9.A4Q8.01 TCGA.NG.A4VU.01 TCGA.N9.A4Q4.01 TCGA.NF.A4WU.01 TCGA.N8.A4PM.01 TCGA.N9.A4Q7.01 
#&gt;           &quot;012&quot;           &quot;022&quot;           &quot;022&quot;           &quot;013&quot;            &quot;03&quot;           &quot;012&quot; 
#&gt; TCGA.NA.A4R1.01 TCGA.N5.A4RJ.01 TCGA.N5.A4RU.01 TCGA.N6.A4V9.01 TCGA.ND.A4W6.01 TCGA.N6.A4VE.01 
#&gt;            &quot;03&quot;            &quot;03&quot;           &quot;023&quot;            &quot;03&quot;           &quot;013&quot;            &quot;03&quot; 
#&gt; TCGA.N7.A4Y8.01 TCGA.N5.A4RM.01 TCGA.NF.A4X2.01 
#&gt;           &quot;021&quot;           &quot;021&quot;            &quot;03&quot;
</code></pre>

<script>
$('#tab-get-classes-from-hierarchical-partition-1-a').parent().next().next().hide();
$('#tab-get-classes-from-hierarchical-partition-1-a').click(function(){
  $('#tab-get-classes-from-hierarchical-partition-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-get-classes-from-hierarchical-partition-2'>
<p><a id='tab-get-classes-from-hierarchical-partition-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">get_classes(res_rh, merge_node = merge_node_param(min_n_signatures = 5323))
</code></pre>

<pre><code>#&gt; TCGA.N6.A4VF.01 TCGA.N7.A59B.01 TCGA.N7.A4Y0.01 TCGA.QN.A5NN.01 TCGA.N6.A4VD.01 TCGA.N8.A4PN.01 
#&gt;            &quot;03&quot;            &quot;01&quot;            &quot;03&quot;            &quot;03&quot;            &quot;03&quot;           &quot;022&quot; 
#&gt; TCGA.N5.A59E.01 TCGA.NG.A4VW.01 TCGA.N8.A4PO.01 TCGA.N7.A4Y5.01 TCGA.ND.A4WA.01 TCGA.N6.A4VC.01 
#&gt;           &quot;021&quot;            &quot;01&quot;            &quot;01&quot;           &quot;022&quot;            &quot;01&quot;            &quot;01&quot; 
#&gt; TCGA.N6.A4VG.01 TCGA.N5.A59F.01 TCGA.N5.A4RD.01 TCGA.N9.A4Q3.01 TCGA.ND.A4WC.01 TCGA.N8.A4PI.01 
#&gt;            &quot;01&quot;            &quot;03&quot;           &quot;021&quot;            &quot;01&quot;            &quot;01&quot;            &quot;03&quot; 
#&gt; TCGA.NA.A5I1.01 TCGA.N8.A56S.01 TCGA.N8.A4PL.01 TCGA.NA.A4R0.01 TCGA.N8.A4PQ.01 TCGA.NA.A4QV.01 
#&gt;            &quot;01&quot;            &quot;01&quot;            &quot;01&quot;            &quot;01&quot;           &quot;022&quot;            &quot;03&quot; 
#&gt; TCGA.NF.A5CP.01 TCGA.ND.A4WF.01 TCGA.N5.A4R8.01 TCGA.NA.A4QY.01 TCGA.NA.A4QW.01 TCGA.N9.A4Q1.01 
#&gt;            &quot;03&quot;           &quot;021&quot;           &quot;023&quot;            &quot;01&quot;           &quot;022&quot;           &quot;022&quot; 
#&gt; TCGA.N5.A4RV.01 TCGA.NF.A4WX.01 TCGA.N8.A4PP.01 TCGA.N5.A4RS.01 TCGA.N5.A4RN.01 TCGA.N5.A4RT.01 
#&gt;           &quot;023&quot;            &quot;01&quot;            &quot;01&quot;            &quot;01&quot;            &quot;03&quot;            &quot;03&quot; 
#&gt; TCGA.N5.A4RF.01 TCGA.QM.A5NM.01 TCGA.N9.A4PZ.01 TCGA.N5.A4RO.01 TCGA.NA.A4QX.01 TCGA.N5.A4RA.01 
#&gt;            &quot;03&quot;            &quot;03&quot;            &quot;03&quot;           &quot;021&quot;           &quot;022&quot;           &quot;022&quot; 
#&gt; TCGA.N9.A4Q8.01 TCGA.NG.A4VU.01 TCGA.N9.A4Q4.01 TCGA.NF.A4WU.01 TCGA.N8.A4PM.01 TCGA.N9.A4Q7.01 
#&gt;            &quot;01&quot;           &quot;022&quot;           &quot;022&quot;            &quot;01&quot;            &quot;03&quot;            &quot;01&quot; 
#&gt; TCGA.NA.A4R1.01 TCGA.N5.A4RJ.01 TCGA.N5.A4RU.01 TCGA.N6.A4V9.01 TCGA.ND.A4W6.01 TCGA.N6.A4VE.01 
#&gt;            &quot;03&quot;            &quot;03&quot;           &quot;023&quot;            &quot;03&quot;            &quot;01&quot;            &quot;03&quot; 
#&gt; TCGA.N7.A4Y8.01 TCGA.N5.A4RM.01 TCGA.NF.A4X2.01 
#&gt;           &quot;021&quot;           &quot;021&quot;            &quot;03&quot;
</code></pre>

<script>
$('#tab-get-classes-from-hierarchical-partition-2-a').parent().next().next().hide();
$('#tab-get-classes-from-hierarchical-partition-2-a').click(function(){
  $('#tab-get-classes-from-hierarchical-partition-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-get-classes-from-hierarchical-partition-3'>
<p><a id='tab-get-classes-from-hierarchical-partition-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">get_classes(res_rh, merge_node = merge_node_param(min_n_signatures = 17040))
</code></pre>

<pre><code>#&gt; TCGA.N6.A4VF.01 TCGA.N7.A59B.01 TCGA.N7.A4Y0.01 TCGA.QN.A5NN.01 TCGA.N6.A4VD.01 TCGA.N8.A4PN.01 
#&gt;            &quot;03&quot;            &quot;01&quot;            &quot;03&quot;            &quot;03&quot;            &quot;03&quot;            &quot;02&quot; 
#&gt; TCGA.N5.A59E.01 TCGA.NG.A4VW.01 TCGA.N8.A4PO.01 TCGA.N7.A4Y5.01 TCGA.ND.A4WA.01 TCGA.N6.A4VC.01 
#&gt;            &quot;02&quot;            &quot;01&quot;            &quot;01&quot;            &quot;02&quot;            &quot;01&quot;            &quot;01&quot; 
#&gt; TCGA.N6.A4VG.01 TCGA.N5.A59F.01 TCGA.N5.A4RD.01 TCGA.N9.A4Q3.01 TCGA.ND.A4WC.01 TCGA.N8.A4PI.01 
#&gt;            &quot;01&quot;            &quot;03&quot;            &quot;02&quot;            &quot;01&quot;            &quot;01&quot;            &quot;03&quot; 
#&gt; TCGA.NA.A5I1.01 TCGA.N8.A56S.01 TCGA.N8.A4PL.01 TCGA.NA.A4R0.01 TCGA.N8.A4PQ.01 TCGA.NA.A4QV.01 
#&gt;            &quot;01&quot;            &quot;01&quot;            &quot;01&quot;            &quot;01&quot;            &quot;02&quot;            &quot;03&quot; 
#&gt; TCGA.NF.A5CP.01 TCGA.ND.A4WF.01 TCGA.N5.A4R8.01 TCGA.NA.A4QY.01 TCGA.NA.A4QW.01 TCGA.N9.A4Q1.01 
#&gt;            &quot;03&quot;            &quot;02&quot;            &quot;02&quot;            &quot;01&quot;            &quot;02&quot;            &quot;02&quot; 
#&gt; TCGA.N5.A4RV.01 TCGA.NF.A4WX.01 TCGA.N8.A4PP.01 TCGA.N5.A4RS.01 TCGA.N5.A4RN.01 TCGA.N5.A4RT.01 
#&gt;            &quot;02&quot;            &quot;01&quot;            &quot;01&quot;            &quot;01&quot;            &quot;03&quot;            &quot;03&quot; 
#&gt; TCGA.N5.A4RF.01 TCGA.QM.A5NM.01 TCGA.N9.A4PZ.01 TCGA.N5.A4RO.01 TCGA.NA.A4QX.01 TCGA.N5.A4RA.01 
#&gt;            &quot;03&quot;            &quot;03&quot;            &quot;03&quot;            &quot;02&quot;            &quot;02&quot;            &quot;02&quot; 
#&gt; TCGA.N9.A4Q8.01 TCGA.NG.A4VU.01 TCGA.N9.A4Q4.01 TCGA.NF.A4WU.01 TCGA.N8.A4PM.01 TCGA.N9.A4Q7.01 
#&gt;            &quot;01&quot;            &quot;02&quot;            &quot;02&quot;            &quot;01&quot;            &quot;03&quot;            &quot;01&quot; 
#&gt; TCGA.NA.A4R1.01 TCGA.N5.A4RJ.01 TCGA.N5.A4RU.01 TCGA.N6.A4V9.01 TCGA.ND.A4W6.01 TCGA.N6.A4VE.01 
#&gt;            &quot;03&quot;            &quot;03&quot;            &quot;02&quot;            &quot;03&quot;            &quot;01&quot;            &quot;03&quot; 
#&gt; TCGA.N7.A4Y8.01 TCGA.N5.A4RM.01 TCGA.NF.A4X2.01 
#&gt;            &quot;02&quot;            &quot;02&quot;            &quot;03&quot;
</code></pre>

<script>
$('#tab-get-classes-from-hierarchical-partition-3-a').parent().next().next().hide();
$('#tab-get-classes-from-hierarchical-partition-3-a').click(function(){
  $('#tab-get-classes-from-hierarchical-partition-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>



### Top rows heatmap

Heatmaps of the top rows:





```r
top_rows_heatmap(res_rh)
```

![plot of chunk top-rows-heatmap](figure_cola/top-rows-heatmap-1.png)

```
#> Error in h(simpleError(msg, call)) : 
#>   error in evaluating the argument 'object' in selecting a method for function 'draw': no applicable method for 'height' applied to an object of class "list"
```

Top rows on each node:


```r
top_rows_overlap(res_rh, method = "upset")
```

![plot of chunk top-rows-overlap](figure_cola/top-rows-overlap-1.png)


### UMAP plot

UMAP plot which shows how samples are separated.




<script>
$( function() {
	$( '#tabs-dimension-reduction-by-depth' ).tabs();
} );
</script>
<div id='tabs-dimension-reduction-by-depth'>
<ul>
<li><a href='#tab-dimension-reduction-by-depth-1'>n_signatures ≥ 4690</a></li>
<li><a href='#tab-dimension-reduction-by-depth-2'>n_signatures ≥ 5323</a></li>
<li><a href='#tab-dimension-reduction-by-depth-3'>n_signatures ≥ 17040</a></li>
</ul>
<div id='tab-dimension-reduction-by-depth-1'>
<pre><code class="r">par(mfrow = c(1, 2))
dimension_reduction(res_rh, merge_node = merge_node_param(min_n_signatures = 4690),
    method = &quot;UMAP&quot;, top_value_method = &quot;SD&quot;, top_n = 40000, scale_rows = FALSE)
dimension_reduction(res_rh, merge_node = merge_node_param(min_n_signatures = 4690),
    method = &quot;UMAP&quot;, top_value_method = &quot;ATC&quot;, top_n = 40000, scale_rows = TRUE)
</code></pre>

<p><img src="figure_cola/tab-dimension-reduction-by-depth-1-1.png" title="plot of chunk tab-dimension-reduction-by-depth-1" alt="plot of chunk tab-dimension-reduction-by-depth-1" width="100%" /></p>

</div>
<div id='tab-dimension-reduction-by-depth-2'>
<pre><code class="r">par(mfrow = c(1, 2))
dimension_reduction(res_rh, merge_node = merge_node_param(min_n_signatures = 5323),
    method = &quot;UMAP&quot;, top_value_method = &quot;SD&quot;, top_n = 40000, scale_rows = FALSE)
dimension_reduction(res_rh, merge_node = merge_node_param(min_n_signatures = 5323),
    method = &quot;UMAP&quot;, top_value_method = &quot;ATC&quot;, top_n = 40000, scale_rows = TRUE)
</code></pre>

<p><img src="figure_cola/tab-dimension-reduction-by-depth-2-1.png" title="plot of chunk tab-dimension-reduction-by-depth-2" alt="plot of chunk tab-dimension-reduction-by-depth-2" width="100%" /></p>

</div>
<div id='tab-dimension-reduction-by-depth-3'>
<pre><code class="r">par(mfrow = c(1, 2))
dimension_reduction(res_rh, merge_node = merge_node_param(min_n_signatures = 17040),
    method = &quot;UMAP&quot;, top_value_method = &quot;SD&quot;, top_n = 40000, scale_rows = FALSE)
dimension_reduction(res_rh, merge_node = merge_node_param(min_n_signatures = 17040),
    method = &quot;UMAP&quot;, top_value_method = &quot;ATC&quot;, top_n = 40000, scale_rows = TRUE)
</code></pre>

<p><img src="figure_cola/tab-dimension-reduction-by-depth-3-1.png" title="plot of chunk tab-dimension-reduction-by-depth-3" alt="plot of chunk tab-dimension-reduction-by-depth-3" width="100%" /></p>

</div>
</div>




### Signature heatmap

Signatures on the heatmap are the union of all signatures found on every node
on the hierarchy. The number of k-means on rows are automatically selected by the function.




<script>
$( function() {
	$( '#tabs-get-signatures-from-hierarchical-partition' ).tabs();
} );
</script>
<div id='tabs-get-signatures-from-hierarchical-partition'>
<ul>
<li><a href='#tab-get-signatures-from-hierarchical-partition-1'>n_signatures ≥ 4690</a></li>
<li><a href='#tab-get-signatures-from-hierarchical-partition-2'>n_signatures ≥ 5323</a></li>
<li><a href='#tab-get-signatures-from-hierarchical-partition-3'>n_signatures ≥ 17040</a></li>
</ul>
<div id='tab-get-signatures-from-hierarchical-partition-1'>
<pre><code class="r">get_signatures(res_rh, merge_node = merge_node_param(min_n_signatures = 4690))
</code></pre>

<p><img src="figure_cola/tab-get-signatures-from-hierarchical-partition-1-1.png" alt="plot of chunk tab-get-signatures-from-hierarchical-partition-1"/></p>

</div>
<div id='tab-get-signatures-from-hierarchical-partition-2'>
<pre><code class="r">get_signatures(res_rh, merge_node = merge_node_param(min_n_signatures = 5323))
</code></pre>

<p><img src="figure_cola/tab-get-signatures-from-hierarchical-partition-2-1.png" alt="plot of chunk tab-get-signatures-from-hierarchical-partition-2"/></p>

</div>
<div id='tab-get-signatures-from-hierarchical-partition-3'>
<pre><code class="r">get_signatures(res_rh, merge_node = merge_node_param(min_n_signatures = 17040))
</code></pre>

<p><img src="figure_cola/tab-get-signatures-from-hierarchical-partition-3-1.png" alt="plot of chunk tab-get-signatures-from-hierarchical-partition-3"/></p>

</div>
</div>




Compare signatures from different nodes:


```r
compare_signatures(res_rh, verbose = FALSE)
```

![plot of chunk unnamed-chunk-24](figure_cola/unnamed-chunk-24-1.png)

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs. Note it only works on every node and the final signatures
are the union of all signatures of all nodes.


```r
# code only for demonstration
# e.g. to show the top 500 most significant rows on each node.
tb = get_signature(res_rh, top_signatures = 500)
```


## Results for each node


---------------------------------------------------




### Node0


Child nodes: 
                [Node01](#Node01)
        ,
                [Node02](#Node02)
        ,
                Node03-leaf
        .







The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

```r
res = res_rh["0"]
```

A summary of `res` and all the functions that can be applied to it:

```r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4, 5, 6, 7, 8.
#>   On a matrix with 30000 rows and 57 columns.
#>   Top rows (1000) are extracted by 'ATC' method.
#>   Subgroups are detected by 'kmeans' method.
#>   Performed in total 350 partitions by row resampling.
#>   Best k for subgroups seems to be 3.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_partitions"     
#>  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
#> [10] "functional_enrichment"   "get_anno_col"            "get_anno"               
#> [13] "get_classes"             "get_consensus"           "get_matrix"             
#> [16] "get_membership"          "get_param"               "get_signatures"         
#> [19] "get_stats"               "is_best_k"               "is_stable_k"            
#> [22] "membership_heatmap"      "ncol"                    "nrow"                   
#> [25] "plot_ecdf"               "predict_classes"         "rownames"               
#> [28] "select_partition_number" "show"                    "suggest_best_k"         
#> [31] "test_to_known_factors"   "top_rows_heatmap"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

```r
collect_plots(res)
```

![plot of chunk node-0-collect-plots](figure_cola/node-0-collect-plots-1.png)

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

```r
select_partition_number(res)
```

![plot of chunk node-0-select-partition-number](figure_cola/node-0-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

```r
get_stats(res)
```

```
#>   k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2 2 1.000           0.982       0.992         0.4993 0.499   0.499
#> 3 3 1.000           0.976       0.990         0.3576 0.718   0.489
#> 4 4 0.784           0.695       0.841         0.0896 0.977   0.930
#> 5 5 0.789           0.785       0.831         0.0733 0.836   0.499
#> 6 6 0.760           0.768       0.837         0.0363 0.930   0.672
#> 7 7 0.808           0.720       0.814         0.0337 0.962   0.791
#> 8 8 0.829           0.628       0.795         0.0212 0.966   0.811
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

```r
suggest_best_k(res)
```

```
#> [1] 3
#> attr(,"optional")
#> [1] 2
```

There is also optional best $k$ = 2 that is worth to check.

Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-node-0-get-classes' ).tabs();
} );
</script>
<div id='tabs-node-0-get-classes'>
<ul>
<li><a href='#tab-node-0-get-classes-1'>k = 2</a></li>
<li><a href='#tab-node-0-get-classes-2'>k = 3</a></li>
<li><a href='#tab-node-0-get-classes-3'>k = 4</a></li>
<li><a href='#tab-node-0-get-classes-4'>k = 5</a></li>
<li><a href='#tab-node-0-get-classes-5'>k = 6</a></li>
<li><a href='#tab-node-0-get-classes-6'>k = 7</a></li>
<li><a href='#tab-node-0-get-classes-7'>k = 8</a></li>
</ul>

<div id='tab-node-0-get-classes-1'>
<p><a id='tab-node-0-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1   p2
#&gt; TCGA.N6.A4VF.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N7.A59B.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N7.A4Y0.01     2   0.141      0.968 0.02 0.98
#&gt; TCGA.QN.A5NN.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N6.A4VD.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N8.A4PN.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.N5.A59E.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.NG.A4VW.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N8.A4PO.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N7.A4Y5.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.ND.A4WA.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N6.A4VC.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N6.A4VG.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N5.A59F.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N5.A4RD.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.N9.A4Q3.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.ND.A4WC.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N8.A4PI.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.NA.A5I1.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N8.A56S.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N8.A4PL.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.NA.A4R0.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N8.A4PQ.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.NA.A4QV.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.NF.A5CP.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.ND.A4WF.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.N5.A4R8.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.NA.A4QY.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.NA.A4QW.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.N9.A4Q1.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.N5.A4RV.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.NF.A4WX.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N8.A4PP.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N5.A4RS.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N5.A4RN.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.N5.A4RT.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.N5.A4RF.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.QM.A5NM.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N9.A4PZ.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N5.A4RO.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.NA.A4QX.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.N5.A4RA.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.N9.A4Q8.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.NG.A4VU.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.N9.A4Q4.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.NF.A4WU.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N8.A4PM.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N9.A4Q7.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.NA.A4R1.01     2   0.141      0.968 0.02 0.98
#&gt; TCGA.N5.A4RJ.01     2   0.943      0.446 0.36 0.64
#&gt; TCGA.N5.A4RU.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.N6.A4V9.01     2   0.141      0.968 0.02 0.98
#&gt; TCGA.ND.A4W6.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N6.A4VE.01     1   0.000      1.000 1.00 0.00
#&gt; TCGA.N7.A4Y8.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.N5.A4RM.01     2   0.000      0.982 0.00 1.00
#&gt; TCGA.NF.A4X2.01     2   0.141      0.968 0.02 0.98
</code></pre>

<script>
$('#tab-node-0-get-classes-1-a').parent().next().next().hide();
$('#tab-node-0-get-classes-1-a').click(function(){
  $('#tab-node-0-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-0-get-classes-2'>
<p><a id='tab-node-0-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1   p2   p3
#&gt; TCGA.N6.A4VF.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.N7.A59B.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.N7.A4Y0.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.QN.A5NN.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.N6.A4VD.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.N8.A4PN.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.N5.A59E.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.NG.A4VW.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.N8.A4PO.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.N7.A4Y5.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.ND.A4WA.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.N6.A4VC.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.N6.A4VG.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.N5.A59F.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.N5.A4RD.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.N9.A4Q3.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.ND.A4WC.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.N8.A4PI.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.NA.A5I1.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.N8.A56S.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.N8.A4PL.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.NA.A4R0.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.N8.A4PQ.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.NA.A4QV.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.NF.A5CP.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.ND.A4WF.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.N5.A4R8.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.NA.A4QY.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.NA.A4QW.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.N9.A4Q1.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.N5.A4RV.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.NF.A4WX.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.N8.A4PP.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.N5.A4RS.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.N5.A4RN.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.N5.A4RT.01     3  0.5016      0.673 0.00 0.24 0.76
#&gt; TCGA.N5.A4RF.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.QM.A5NM.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.N9.A4PZ.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.N5.A4RO.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.NA.A4QX.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.N5.A4RA.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.N9.A4Q8.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.NG.A4VU.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.N9.A4Q4.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.NF.A4WU.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.N8.A4PM.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.N9.A4Q7.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.NA.A4R1.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.N5.A4RJ.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.N5.A4RU.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.N6.A4V9.01     3  0.0000      0.985 0.00 0.00 1.00
#&gt; TCGA.ND.A4W6.01     1  0.0000      1.000 1.00 0.00 0.00
#&gt; TCGA.N6.A4VE.01     3  0.0892      0.966 0.02 0.00 0.98
#&gt; TCGA.N7.A4Y8.01     2  0.0000      0.982 0.00 1.00 0.00
#&gt; TCGA.N5.A4RM.01     2  0.5560      0.558 0.00 0.70 0.30
#&gt; TCGA.NF.A4X2.01     3  0.0000      0.985 0.00 0.00 1.00
</code></pre>

<script>
$('#tab-node-0-get-classes-2-a').parent().next().next().hide();
$('#tab-node-0-get-classes-2-a').click(function(){
  $('#tab-node-0-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-0-get-classes-3'>
<p><a id='tab-node-0-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1   p2   p3   p4
#&gt; TCGA.N6.A4VF.01     3  0.0000     0.8165 0.00 0.00 1.00 0.00
#&gt; TCGA.N7.A59B.01     1  0.0000     0.8411 1.00 0.00 0.00 0.00
#&gt; TCGA.N7.A4Y0.01     3  0.4624     0.4742 0.00 0.00 0.66 0.34
#&gt; TCGA.QN.A5NN.01     3  0.0707     0.8161 0.00 0.00 0.98 0.02
#&gt; TCGA.N6.A4VD.01     3  0.0707     0.8161 0.00 0.00 0.98 0.02
#&gt; TCGA.N8.A4PN.01     2  0.3172     0.7126 0.00 0.84 0.00 0.16
#&gt; TCGA.N5.A59E.01     2  0.0707     0.7999 0.00 0.98 0.00 0.02
#&gt; TCGA.NG.A4VW.01     1  0.0000     0.8411 1.00 0.00 0.00 0.00
#&gt; TCGA.N8.A4PO.01     1  0.0000     0.8411 1.00 0.00 0.00 0.00
#&gt; TCGA.N7.A4Y5.01     2  0.0000     0.8099 0.00 1.00 0.00 0.00
#&gt; TCGA.ND.A4WA.01     1  0.0000     0.8411 1.00 0.00 0.00 0.00
#&gt; TCGA.N6.A4VC.01     1  0.0000     0.8411 1.00 0.00 0.00 0.00
#&gt; TCGA.N6.A4VG.01     1  0.2647     0.8321 0.88 0.00 0.00 0.12
#&gt; TCGA.N5.A59F.01     3  0.0000     0.8165 0.00 0.00 1.00 0.00
#&gt; TCGA.N5.A4RD.01     2  0.0707     0.7999 0.00 0.98 0.00 0.02
#&gt; TCGA.N9.A4Q3.01     1  0.4790     0.7802 0.62 0.00 0.00 0.38
#&gt; TCGA.ND.A4WC.01     1  0.0000     0.8411 1.00 0.00 0.00 0.00
#&gt; TCGA.N8.A4PI.01     3  0.0707     0.8161 0.00 0.00 0.98 0.02
#&gt; TCGA.NA.A5I1.01     1  0.4790     0.7802 0.62 0.00 0.00 0.38
#&gt; TCGA.N8.A56S.01     1  0.4855     0.7689 0.60 0.00 0.00 0.40
#&gt; TCGA.N8.A4PL.01     1  0.4790     0.7802 0.62 0.00 0.00 0.38
#&gt; TCGA.NA.A4R0.01     1  0.0000     0.8411 1.00 0.00 0.00 0.00
#&gt; TCGA.N8.A4PQ.01     2  0.5428     0.4100 0.00 0.60 0.02 0.38
#&gt; TCGA.NA.A4QV.01     3  0.0000     0.8165 0.00 0.00 1.00 0.00
#&gt; TCGA.NF.A5CP.01     3  0.0707     0.8161 0.00 0.00 0.98 0.02
#&gt; TCGA.ND.A4WF.01     2  0.0000     0.8099 0.00 1.00 0.00 0.00
#&gt; TCGA.N5.A4R8.01     2  0.0000     0.8099 0.00 1.00 0.00 0.00
#&gt; TCGA.NA.A4QY.01     1  0.4790     0.7802 0.62 0.00 0.00 0.38
#&gt; TCGA.NA.A4QW.01     2  0.3172     0.7126 0.00 0.84 0.00 0.16
#&gt; TCGA.N9.A4Q1.01     2  0.0000     0.8099 0.00 1.00 0.00 0.00
#&gt; TCGA.N5.A4RV.01     2  0.4907    -0.0762 0.00 0.58 0.00 0.42
#&gt; TCGA.NF.A4WX.01     1  0.4790     0.7802 0.62 0.00 0.00 0.38
#&gt; TCGA.N8.A4PP.01     1  0.0000     0.8411 1.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RS.01     1  0.4855     0.7689 0.60 0.00 0.00 0.40
#&gt; TCGA.N5.A4RN.01     3  0.4948     0.1965 0.00 0.00 0.56 0.44
#&gt; TCGA.N5.A4RT.01     4  0.7602     0.3923 0.00 0.20 0.38 0.42
#&gt; TCGA.N5.A4RF.01     3  0.0000     0.8165 0.00 0.00 1.00 0.00
#&gt; TCGA.QM.A5NM.01     3  0.0707     0.8161 0.00 0.00 0.98 0.02
#&gt; TCGA.N9.A4PZ.01     3  0.0000     0.8165 0.00 0.00 1.00 0.00
#&gt; TCGA.N5.A4RO.01     2  0.4907    -0.0762 0.00 0.58 0.00 0.42
#&gt; TCGA.NA.A4QX.01     2  0.0000     0.8099 0.00 1.00 0.00 0.00
#&gt; TCGA.N5.A4RA.01     2  0.3172     0.7126 0.00 0.84 0.00 0.16
#&gt; TCGA.N9.A4Q8.01     1  0.0000     0.8411 1.00 0.00 0.00 0.00
#&gt; TCGA.NG.A4VU.01     2  0.0000     0.8099 0.00 1.00 0.00 0.00
#&gt; TCGA.N9.A4Q4.01     2  0.0000     0.8099 0.00 1.00 0.00 0.00
#&gt; TCGA.NF.A4WU.01     1  0.2647     0.8321 0.88 0.00 0.00 0.12
#&gt; TCGA.N8.A4PM.01     3  0.0707     0.8051 0.00 0.00 0.98 0.02
#&gt; TCGA.N9.A4Q7.01     1  0.0000     0.8411 1.00 0.00 0.00 0.00
#&gt; TCGA.NA.A4R1.01     3  0.4624     0.4742 0.00 0.00 0.66 0.34
#&gt; TCGA.N5.A4RJ.01     3  0.4134     0.5835 0.00 0.00 0.74 0.26
#&gt; TCGA.N5.A4RU.01     2  0.4907    -0.0762 0.00 0.58 0.00 0.42
#&gt; TCGA.N6.A4V9.01     3  0.4948     0.1472 0.00 0.00 0.56 0.44
#&gt; TCGA.ND.A4W6.01     1  0.4855     0.7689 0.60 0.00 0.00 0.40
#&gt; TCGA.N6.A4VE.01     3  0.1637     0.7606 0.06 0.00 0.94 0.00
#&gt; TCGA.N7.A4Y8.01     2  0.0000     0.8099 0.00 1.00 0.00 0.00
#&gt; TCGA.N5.A4RM.01     4  0.7220     0.3681 0.00 0.42 0.14 0.44
#&gt; TCGA.NF.A4X2.01     3  0.4277     0.5748 0.00 0.00 0.72 0.28
</code></pre>

<script>
$('#tab-node-0-get-classes-3-a').parent().next().next().hide();
$('#tab-node-0-get-classes-3-a').click(function(){
  $('#tab-node-0-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-0-get-classes-4'>
<p><a id='tab-node-0-get-classes-4-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 5), get_membership(res, k = 5))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1   p2   p3   p4   p5
#&gt; TCGA.N6.A4VF.01     3  0.0000      0.916 0.00 0.00 1.00 0.00 0.00
#&gt; TCGA.N7.A59B.01     1  0.4748      0.774 0.66 0.00 0.04 0.00 0.30
#&gt; TCGA.N7.A4Y0.01     4  0.3561      0.569 0.00 0.00 0.26 0.74 0.00
#&gt; TCGA.QN.A5NN.01     3  0.3424      0.708 0.00 0.00 0.76 0.24 0.00
#&gt; TCGA.N6.A4VD.01     3  0.0609      0.916 0.00 0.00 0.98 0.02 0.00
#&gt; TCGA.N8.A4PN.01     2  0.2754      0.885 0.08 0.88 0.00 0.04 0.00
#&gt; TCGA.N5.A59E.01     2  0.4263      0.720 0.18 0.76 0.00 0.06 0.00
#&gt; TCGA.NG.A4VW.01     1  0.4613      0.781 0.62 0.00 0.00 0.02 0.36
#&gt; TCGA.N8.A4PO.01     1  0.4748      0.774 0.66 0.00 0.04 0.00 0.30
#&gt; TCGA.N7.A4Y5.01     2  0.0000      0.936 0.00 1.00 0.00 0.00 0.00
#&gt; TCGA.ND.A4WA.01     1  0.3983      0.810 0.66 0.00 0.00 0.00 0.34
#&gt; TCGA.N6.A4VC.01     1  0.3983      0.810 0.66 0.00 0.00 0.00 0.34
#&gt; TCGA.N6.A4VG.01     5  0.4252      0.357 0.28 0.00 0.00 0.02 0.70
#&gt; TCGA.N5.A59F.01     3  0.1648      0.899 0.04 0.00 0.94 0.02 0.00
#&gt; TCGA.N5.A4RD.01     2  0.2077      0.893 0.04 0.92 0.00 0.04 0.00
#&gt; TCGA.N9.A4Q3.01     5  0.0609      0.864 0.02 0.00 0.00 0.00 0.98
#&gt; TCGA.ND.A4WC.01     1  0.4540      0.805 0.64 0.00 0.00 0.02 0.34
#&gt; TCGA.N8.A4PI.01     3  0.0609      0.916 0.00 0.00 0.98 0.02 0.00
#&gt; TCGA.NA.A5I1.01     5  0.0000      0.872 0.00 0.00 0.00 0.00 1.00
#&gt; TCGA.N8.A56S.01     5  0.1043      0.861 0.00 0.00 0.00 0.04 0.96
#&gt; TCGA.N8.A4PL.01     5  0.0000      0.872 0.00 0.00 0.00 0.00 1.00
#&gt; TCGA.NA.A4R0.01     1  0.3983      0.810 0.66 0.00 0.00 0.00 0.34
#&gt; TCGA.N8.A4PQ.01     1  0.7761     -0.465 0.34 0.32 0.02 0.30 0.02
#&gt; TCGA.NA.A4QV.01     3  0.0609      0.914 0.02 0.00 0.98 0.00 0.00
#&gt; TCGA.NF.A5CP.01     3  0.0609      0.916 0.00 0.00 0.98 0.02 0.00
#&gt; TCGA.ND.A4WF.01     2  0.2077      0.893 0.04 0.92 0.00 0.04 0.00
#&gt; TCGA.N5.A4R8.01     2  0.0000      0.936 0.00 1.00 0.00 0.00 0.00
#&gt; TCGA.NA.A4QY.01     5  0.0609      0.864 0.02 0.00 0.00 0.00 0.98
#&gt; TCGA.NA.A4QW.01     2  0.2754      0.885 0.08 0.88 0.00 0.04 0.00
#&gt; TCGA.N9.A4Q1.01     2  0.0000      0.936 0.00 1.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RV.01     4  0.5597      0.707 0.20 0.16 0.00 0.64 0.00
#&gt; TCGA.NF.A4WX.01     5  0.0000      0.872 0.00 0.00 0.00 0.00 1.00
#&gt; TCGA.N8.A4PP.01     1  0.3983      0.810 0.66 0.00 0.00 0.00 0.34
#&gt; TCGA.N5.A4RS.01     5  0.1043      0.861 0.00 0.00 0.00 0.04 0.96
#&gt; TCGA.N5.A4RN.01     4  0.2020      0.705 0.00 0.00 0.10 0.90 0.00
#&gt; TCGA.N5.A4RT.01     4  0.5829      0.762 0.22 0.04 0.08 0.66 0.00
#&gt; TCGA.N5.A4RF.01     3  0.1648      0.899 0.04 0.00 0.94 0.02 0.00
#&gt; TCGA.QM.A5NM.01     3  0.1216      0.908 0.02 0.00 0.96 0.02 0.00
#&gt; TCGA.N9.A4PZ.01     3  0.0609      0.910 0.02 0.00 0.98 0.00 0.00
#&gt; TCGA.N5.A4RO.01     4  0.5759      0.688 0.20 0.18 0.00 0.62 0.00
#&gt; TCGA.NA.A4QX.01     2  0.0000      0.936 0.00 1.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RA.01     2  0.2754      0.885 0.08 0.88 0.00 0.04 0.00
#&gt; TCGA.N9.A4Q8.01     1  0.4540      0.805 0.64 0.00 0.00 0.02 0.34
#&gt; TCGA.NG.A4VU.01     2  0.0000      0.936 0.00 1.00 0.00 0.00 0.00
#&gt; TCGA.N9.A4Q4.01     2  0.0000      0.936 0.00 1.00 0.00 0.00 0.00
#&gt; TCGA.NF.A4WU.01     5  0.4252      0.357 0.28 0.00 0.00 0.02 0.70
#&gt; TCGA.N8.A4PM.01     3  0.1043      0.908 0.04 0.00 0.96 0.00 0.00
#&gt; TCGA.N9.A4Q7.01     1  0.4540      0.805 0.64 0.00 0.00 0.02 0.34
#&gt; TCGA.NA.A4R1.01     4  0.3561      0.569 0.00 0.00 0.26 0.74 0.00
#&gt; TCGA.N5.A4RJ.01     3  0.5013      0.475 0.08 0.00 0.68 0.24 0.00
#&gt; TCGA.N5.A4RU.01     4  0.5597      0.707 0.20 0.16 0.00 0.64 0.00
#&gt; TCGA.N6.A4V9.01     4  0.5599      0.740 0.26 0.00 0.12 0.62 0.00
#&gt; TCGA.ND.A4W6.01     5  0.1648      0.841 0.02 0.00 0.00 0.04 0.94
#&gt; TCGA.N6.A4VE.01     3  0.1043      0.894 0.04 0.00 0.96 0.00 0.00
#&gt; TCGA.N7.A4Y8.01     2  0.0000      0.936 0.00 1.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RM.01     4  0.5958      0.761 0.24 0.04 0.08 0.64 0.00
#&gt; TCGA.NF.A4X2.01     4  0.3684      0.535 0.00 0.00 0.28 0.72 0.00
</code></pre>

<script>
$('#tab-node-0-get-classes-4-a').parent().next().next().hide();
$('#tab-node-0-get-classes-4-a').click(function(){
  $('#tab-node-0-get-classes-4-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-0-get-classes-5'>
<p><a id='tab-node-0-get-classes-5-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 6), get_membership(res, k = 6))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1   p2   p3   p4   p5   p6
#&gt; TCGA.N6.A4VF.01     3  0.0000     0.9100 0.00 0.00 1.00 0.00 0.00 0.00
#&gt; TCGA.N7.A59B.01     1  0.2725     0.7701 0.88 0.00 0.02 0.00 0.04 0.06
#&gt; TCGA.N7.A4Y0.01     6  0.5422     0.8226 0.00 0.00 0.12 0.24 0.02 0.62
#&gt; TCGA.QN.A5NN.01     6  0.4282     0.3900 0.00 0.00 0.42 0.00 0.02 0.56
#&gt; TCGA.N6.A4VD.01     3  0.1092     0.8996 0.00 0.00 0.96 0.00 0.02 0.02
#&gt; TCGA.N8.A4PN.01     2  0.3970     0.7797 0.00 0.80 0.00 0.04 0.08 0.08
#&gt; TCGA.N5.A59E.01     4  0.5424    -0.0705 0.00 0.46 0.00 0.46 0.04 0.04
#&gt; TCGA.NG.A4VW.01     1  0.1092     0.8226 0.96 0.00 0.00 0.00 0.02 0.02
#&gt; TCGA.N8.A4PO.01     1  0.2725     0.7701 0.88 0.00 0.02 0.00 0.04 0.06
#&gt; TCGA.N7.A4Y5.01     2  0.0000     0.8762 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; TCGA.ND.A4WA.01     1  0.2581     0.7752 0.86 0.00 0.00 0.02 0.00 0.12
#&gt; TCGA.N6.A4VC.01     1  0.2581     0.7752 0.86 0.00 0.00 0.02 0.00 0.12
#&gt; TCGA.N6.A4VG.01     1  0.4078     0.2129 0.64 0.00 0.00 0.00 0.34 0.02
#&gt; TCGA.N5.A59F.01     3  0.0547     0.9081 0.00 0.00 0.98 0.00 0.02 0.00
#&gt; TCGA.N5.A4RD.01     2  0.3409     0.5075 0.00 0.70 0.00 0.30 0.00 0.00
#&gt; TCGA.N9.A4Q3.01     5  0.4067     0.9070 0.26 0.00 0.00 0.00 0.70 0.04
#&gt; TCGA.ND.A4WC.01     1  0.0000     0.8434 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N8.A4PI.01     3  0.1092     0.8996 0.00 0.00 0.96 0.00 0.02 0.02
#&gt; TCGA.NA.A5I1.01     5  0.3711     0.9110 0.26 0.00 0.00 0.00 0.72 0.02
#&gt; TCGA.N8.A56S.01     5  0.5007     0.8847 0.24 0.00 0.00 0.02 0.66 0.08
#&gt; TCGA.N8.A4PL.01     5  0.3711     0.9110 0.26 0.00 0.00 0.00 0.72 0.02
#&gt; TCGA.NA.A4R0.01     1  0.0000     0.8434 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N8.A4PQ.01     4  0.6168     0.5035 0.00 0.10 0.00 0.60 0.16 0.14
#&gt; TCGA.NA.A4QV.01     3  0.0547     0.9081 0.00 0.00 0.98 0.00 0.02 0.00
#&gt; TCGA.NF.A5CP.01     3  0.0000     0.9100 0.00 0.00 1.00 0.00 0.00 0.00
#&gt; TCGA.ND.A4WF.01     2  0.3409     0.5075 0.00 0.70 0.00 0.30 0.00 0.00
#&gt; TCGA.N5.A4R8.01     2  0.0000     0.8762 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; TCGA.NA.A4QY.01     5  0.4067     0.9070 0.26 0.00 0.00 0.00 0.70 0.04
#&gt; TCGA.NA.A4QW.01     2  0.3697     0.7895 0.00 0.82 0.00 0.04 0.08 0.06
#&gt; TCGA.N9.A4Q1.01     2  0.0000     0.8762 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RV.01     4  0.1556     0.7387 0.00 0.08 0.00 0.92 0.00 0.00
#&gt; TCGA.NF.A4WX.01     5  0.3711     0.9110 0.26 0.00 0.00 0.00 0.72 0.02
#&gt; TCGA.N8.A4PP.01     1  0.0000     0.8434 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RS.01     5  0.5007     0.8847 0.24 0.00 0.00 0.02 0.66 0.08
#&gt; TCGA.N5.A4RN.01     6  0.4552     0.7579 0.00 0.00 0.06 0.30 0.00 0.64
#&gt; TCGA.N5.A4RT.01     4  0.1807     0.7016 0.00 0.02 0.06 0.92 0.00 0.00
#&gt; TCGA.N5.A4RF.01     3  0.0547     0.9081 0.00 0.00 0.98 0.00 0.02 0.00
#&gt; TCGA.QM.A5NM.01     3  0.2794     0.8425 0.00 0.00 0.86 0.00 0.06 0.08
#&gt; TCGA.N9.A4PZ.01     3  0.1807     0.8743 0.00 0.00 0.92 0.00 0.02 0.06
#&gt; TCGA.N5.A4RO.01     4  0.1814     0.7318 0.00 0.10 0.00 0.90 0.00 0.00
#&gt; TCGA.NA.A4QX.01     2  0.0000     0.8762 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RA.01     2  0.3970     0.7797 0.00 0.80 0.00 0.04 0.08 0.08
#&gt; TCGA.N9.A4Q8.01     1  0.0000     0.8434 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.NG.A4VU.01     2  0.0000     0.8762 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N9.A4Q4.01     2  0.0000     0.8762 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; TCGA.NF.A4WU.01     1  0.4078     0.2129 0.64 0.00 0.00 0.00 0.34 0.02
#&gt; TCGA.N8.A4PM.01     3  0.1092     0.9032 0.00 0.00 0.96 0.00 0.02 0.02
#&gt; TCGA.N9.A4Q7.01     1  0.0000     0.8434 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.NA.A4R1.01     6  0.4834     0.8151 0.00 0.00 0.10 0.26 0.00 0.64
#&gt; TCGA.N5.A4RJ.01     3  0.4094     0.5994 0.00 0.00 0.74 0.18 0.08 0.00
#&gt; TCGA.N5.A4RU.01     4  0.1556     0.7387 0.00 0.08 0.00 0.92 0.00 0.00
#&gt; TCGA.N6.A4V9.01     4  0.3475     0.5631 0.00 0.00 0.14 0.80 0.06 0.00
#&gt; TCGA.ND.A4W6.01     5  0.4887     0.8664 0.22 0.00 0.00 0.02 0.68 0.08
#&gt; TCGA.N6.A4VE.01     3  0.4008     0.7698 0.04 0.00 0.80 0.02 0.02 0.12
#&gt; TCGA.N7.A4Y8.01     2  0.0000     0.8762 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RM.01     4  0.3045     0.6819 0.00 0.02 0.06 0.86 0.06 0.00
#&gt; TCGA.NF.A4X2.01     6  0.4989     0.8196 0.00 0.00 0.14 0.22 0.00 0.64
</code></pre>

<script>
$('#tab-node-0-get-classes-5-a').parent().next().next().hide();
$('#tab-node-0-get-classes-5-a').click(function(){
  $('#tab-node-0-get-classes-5-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-0-get-classes-6'>
<p><a id='tab-node-0-get-classes-6-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 7), get_membership(res, k = 7))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1   p2   p3   p4   p5   p6   p7
#&gt; TCGA.N6.A4VF.01     3  0.0000      0.827 0.00 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N7.A59B.01     1  0.3546      0.840 0.82 0.00 0.02 0.00 0.08 0.02 0.06
#&gt; TCGA.N7.A4Y0.01     6  0.3911      0.828 0.00 0.00 0.06 0.24 0.00 0.70 0.00
#&gt; TCGA.QN.A5NN.01     6  0.3943      0.177 0.00 0.00 0.42 0.00 0.00 0.56 0.02
#&gt; TCGA.N6.A4VD.01     3  0.2376      0.793 0.02 0.00 0.86 0.00 0.00 0.00 0.12
#&gt; TCGA.N8.A4PN.01     2  0.3221      0.614 0.00 0.68 0.00 0.00 0.00 0.00 0.32
#&gt; TCGA.N5.A59E.01     2  0.4820      0.241 0.00 0.52 0.00 0.38 0.00 0.00 0.10
#&gt; TCGA.NG.A4VW.01     1  0.1671      0.841 0.90 0.00 0.00 0.00 0.10 0.00 0.00
#&gt; TCGA.N8.A4PO.01     1  0.3546      0.840 0.82 0.00 0.02 0.00 0.08 0.02 0.06
#&gt; TCGA.N7.A4Y5.01     2  0.0000      0.819 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.ND.A4WA.01     1  0.6405      0.533 0.46 0.00 0.00 0.00 0.10 0.26 0.18
#&gt; TCGA.N6.A4VC.01     1  0.6405      0.533 0.46 0.00 0.00 0.00 0.10 0.26 0.18
#&gt; TCGA.N6.A4VG.01     5  0.3755      0.431 0.34 0.00 0.00 0.00 0.64 0.00 0.02
#&gt; TCGA.N5.A59F.01     3  0.3000      0.800 0.00 0.00 0.84 0.02 0.00 0.04 0.10
#&gt; TCGA.N5.A4RD.01     2  0.2945      0.617 0.00 0.74 0.00 0.26 0.00 0.00 0.00
#&gt; TCGA.N9.A4Q3.01     5  0.4177      0.698 0.00 0.00 0.00 0.00 0.70 0.12 0.18
#&gt; TCGA.ND.A4WC.01     1  0.1433      0.854 0.92 0.00 0.00 0.00 0.08 0.00 0.00
#&gt; TCGA.N8.A4PI.01     3  0.2376      0.793 0.02 0.00 0.86 0.00 0.00 0.00 0.12
#&gt; TCGA.NA.A5I1.01     5  0.0504      0.762 0.02 0.00 0.00 0.00 0.98 0.00 0.00
#&gt; TCGA.N8.A56S.01     5  0.2572      0.740 0.00 0.00 0.00 0.00 0.80 0.00 0.20
#&gt; TCGA.N8.A4PL.01     5  0.0863      0.759 0.04 0.00 0.00 0.00 0.96 0.00 0.00
#&gt; TCGA.NA.A4R0.01     1  0.2803      0.849 0.84 0.00 0.00 0.00 0.10 0.06 0.00
#&gt; TCGA.N8.A4PQ.01     7  0.4408      0.000 0.00 0.02 0.00 0.44 0.02 0.00 0.52
#&gt; TCGA.NA.A4QV.01     3  0.2509      0.813 0.00 0.00 0.88 0.02 0.00 0.04 0.06
#&gt; TCGA.NF.A5CP.01     3  0.0504      0.826 0.02 0.00 0.98 0.00 0.00 0.00 0.00
#&gt; TCGA.ND.A4WF.01     2  0.2945      0.617 0.00 0.74 0.00 0.26 0.00 0.00 0.00
#&gt; TCGA.N5.A4R8.01     2  0.0000      0.819 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.NA.A4QY.01     5  0.4177      0.698 0.00 0.00 0.00 0.00 0.70 0.12 0.18
#&gt; TCGA.NA.A4QW.01     2  0.3047      0.646 0.00 0.72 0.00 0.00 0.00 0.00 0.28
#&gt; TCGA.N9.A4Q1.01     2  0.0000      0.819 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RV.01     4  0.0504      0.827 0.00 0.02 0.00 0.98 0.00 0.00 0.00
#&gt; TCGA.NF.A4WX.01     5  0.0504      0.762 0.02 0.00 0.00 0.00 0.98 0.00 0.00
#&gt; TCGA.N8.A4PP.01     1  0.2803      0.849 0.84 0.00 0.00 0.00 0.10 0.06 0.00
#&gt; TCGA.N5.A4RS.01     5  0.2572      0.740 0.00 0.00 0.00 0.00 0.80 0.00 0.20
#&gt; TCGA.N5.A4RN.01     6  0.3745      0.798 0.00 0.00 0.04 0.26 0.00 0.70 0.00
#&gt; TCGA.N5.A4RT.01     4  0.0863      0.817 0.04 0.00 0.00 0.96 0.00 0.00 0.00
#&gt; TCGA.N5.A4RF.01     3  0.3691      0.790 0.02 0.00 0.80 0.02 0.00 0.04 0.12
#&gt; TCGA.QM.A5NM.01     3  0.2906      0.762 0.02 0.00 0.80 0.00 0.00 0.00 0.18
#&gt; TCGA.N9.A4PZ.01     3  0.2213      0.812 0.02 0.00 0.90 0.00 0.00 0.04 0.04
#&gt; TCGA.N5.A4RO.01     4  0.0863      0.794 0.00 0.04 0.00 0.96 0.00 0.00 0.00
#&gt; TCGA.NA.A4QX.01     2  0.0000      0.819 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RA.01     2  0.3221      0.614 0.00 0.68 0.00 0.00 0.00 0.00 0.32
#&gt; TCGA.N9.A4Q8.01     1  0.1433      0.854 0.92 0.00 0.00 0.00 0.08 0.00 0.00
#&gt; TCGA.NG.A4VU.01     2  0.0000      0.819 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N9.A4Q4.01     2  0.0000      0.819 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.NF.A4WU.01     5  0.3755      0.431 0.34 0.00 0.00 0.00 0.64 0.00 0.02
#&gt; TCGA.N8.A4PM.01     3  0.4440      0.785 0.02 0.00 0.74 0.02 0.00 0.10 0.12
#&gt; TCGA.N9.A4Q7.01     1  0.1433      0.854 0.92 0.00 0.00 0.00 0.08 0.00 0.00
#&gt; TCGA.NA.A4R1.01     6  0.3911      0.828 0.00 0.00 0.06 0.24 0.00 0.70 0.00
#&gt; TCGA.N5.A4RJ.01     3  0.4992      0.613 0.00 0.00 0.66 0.16 0.00 0.04 0.14
#&gt; TCGA.N5.A4RU.01     4  0.0504      0.827 0.00 0.02 0.00 0.98 0.00 0.00 0.00
#&gt; TCGA.N6.A4V9.01     4  0.3879      0.507 0.04 0.00 0.10 0.78 0.00 0.00 0.08
#&gt; TCGA.ND.A4W6.01     5  0.2572      0.740 0.00 0.00 0.00 0.00 0.80 0.00 0.20
#&gt; TCGA.N6.A4VE.01     3  0.3765      0.738 0.02 0.00 0.78 0.00 0.00 0.10 0.10
#&gt; TCGA.N7.A4Y8.01     2  0.0000      0.819 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RM.01     4  0.1718      0.781 0.04 0.00 0.00 0.92 0.00 0.00 0.04
#&gt; TCGA.NF.A4X2.01     6  0.3911      0.828 0.00 0.00 0.06 0.24 0.00 0.70 0.00
</code></pre>

<script>
$('#tab-node-0-get-classes-6-a').parent().next().next().hide();
$('#tab-node-0-get-classes-6-a').click(function(){
  $('#tab-node-0-get-classes-6-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-0-get-classes-7'>
<p><a id='tab-node-0-get-classes-7-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 8), get_membership(res, k = 8))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1   p2   p3   p4   p5   p6   p7   p8
#&gt; TCGA.N6.A4VF.01     3  0.0000      0.722 0.00 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N7.A59B.01     1  0.3036      0.641 0.78 0.00 0.00 0.04 0.00 0.00 0.18 0.00
#&gt; TCGA.N7.A4Y0.01     6  0.2674      0.859 0.00 0.00 0.04 0.00 0.00 0.86 0.04 0.06
#&gt; TCGA.QN.A5NN.01     6  0.3845      0.524 0.00 0.00 0.28 0.00 0.00 0.66 0.06 0.00
#&gt; TCGA.N6.A4VD.01     3  0.4057      0.619 0.00 0.00 0.72 0.20 0.02 0.02 0.04 0.00
#&gt; TCGA.N8.A4PN.01     2  0.3303      0.690 0.00 0.76 0.00 0.20 0.00 0.02 0.02 0.00
#&gt; TCGA.N5.A59E.01     8  0.4672      0.132 0.00 0.24 0.00 0.12 0.00 0.02 0.00 0.62
#&gt; TCGA.NG.A4VW.01     1  0.0941      0.742 0.96 0.00 0.00 0.00 0.02 0.02 0.00 0.00
#&gt; TCGA.N8.A4PO.01     1  0.3036      0.641 0.78 0.00 0.00 0.04 0.00 0.00 0.18 0.00
#&gt; TCGA.N7.A4Y5.01     2  0.0000      0.825 0.00 1.00 0.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.ND.A4WA.01     1  0.3737      0.558 0.50 0.00 0.00 0.48 0.02 0.00 0.00 0.00
#&gt; TCGA.N6.A4VC.01     1  0.3737      0.558 0.50 0.00 0.00 0.48 0.02 0.00 0.00 0.00
#&gt; TCGA.N6.A4VG.01     1  0.4607      0.250 0.50 0.00 0.00 0.06 0.42 0.02 0.00 0.00
#&gt; TCGA.N5.A59F.01     3  0.3198      0.661 0.00 0.00 0.72 0.00 0.00 0.00 0.26 0.02
#&gt; TCGA.N5.A4RD.01     2  0.3971      0.257 0.00 0.54 0.00 0.04 0.00 0.00 0.00 0.42
#&gt; TCGA.N9.A4Q3.01     5  0.4358      0.577 0.04 0.00 0.00 0.40 0.54 0.02 0.00 0.00
#&gt; TCGA.ND.A4WC.01     1  0.0471      0.741 0.98 0.00 0.00 0.00 0.02 0.00 0.00 0.00
#&gt; TCGA.N8.A4PI.01     3  0.4057      0.619 0.00 0.00 0.72 0.20 0.02 0.02 0.04 0.00
#&gt; TCGA.NA.A5I1.01     5  0.0808      0.794 0.04 0.00 0.00 0.00 0.96 0.00 0.00 0.00
#&gt; TCGA.N8.A56S.01     5  0.3262      0.775 0.02 0.00 0.00 0.02 0.82 0.04 0.10 0.00
#&gt; TCGA.N8.A4PL.01     5  0.0808      0.794 0.04 0.00 0.00 0.00 0.96 0.00 0.00 0.00
#&gt; TCGA.NA.A4R0.01     1  0.2859      0.728 0.82 0.00 0.00 0.14 0.02 0.00 0.02 0.00
#&gt; TCGA.N8.A4PQ.01     7  0.6225      0.000 0.00 0.00 0.00 0.20 0.02 0.06 0.42 0.30
#&gt; TCGA.NA.A4QV.01     3  0.2224      0.714 0.00 0.00 0.86 0.00 0.00 0.00 0.12 0.02
#&gt; TCGA.NF.A5CP.01     3  0.0471      0.720 0.00 0.00 0.98 0.02 0.00 0.00 0.00 0.00
#&gt; TCGA.ND.A4WF.01     2  0.3971      0.257 0.00 0.54 0.00 0.04 0.00 0.00 0.00 0.42
#&gt; TCGA.N5.A4R8.01     2  0.0471      0.819 0.00 0.98 0.00 0.02 0.00 0.00 0.00 0.00
#&gt; TCGA.NA.A4QY.01     5  0.4740      0.582 0.04 0.00 0.00 0.38 0.54 0.02 0.02 0.00
#&gt; TCGA.NA.A4QW.01     2  0.3169      0.703 0.00 0.78 0.00 0.18 0.00 0.02 0.02 0.00
#&gt; TCGA.N9.A4Q1.01     2  0.0000      0.825 0.00 1.00 0.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RV.01     8  0.0471      0.629 0.00 0.02 0.00 0.00 0.00 0.00 0.00 0.98
#&gt; TCGA.NF.A4WX.01     5  0.0808      0.794 0.04 0.00 0.00 0.00 0.96 0.00 0.00 0.00
#&gt; TCGA.N8.A4PP.01     1  0.3169      0.720 0.78 0.00 0.00 0.18 0.02 0.00 0.02 0.00
#&gt; TCGA.N5.A4RS.01     5  0.3262      0.775 0.02 0.00 0.00 0.02 0.82 0.04 0.10 0.00
#&gt; TCGA.N5.A4RN.01     6  0.3369      0.820 0.00 0.00 0.02 0.04 0.00 0.82 0.04 0.08
#&gt; TCGA.N5.A4RT.01     8  0.1341      0.573 0.00 0.00 0.00 0.00 0.00 0.00 0.08 0.92
#&gt; TCGA.N5.A4RF.01     3  0.3374      0.632 0.00 0.00 0.68 0.00 0.00 0.00 0.30 0.02
#&gt; TCGA.QM.A5NM.01     3  0.5914      0.471 0.02 0.00 0.52 0.20 0.02 0.02 0.22 0.00
#&gt; TCGA.N9.A4PZ.01     3  0.2719      0.653 0.00 0.00 0.80 0.02 0.00 0.00 0.18 0.00
#&gt; TCGA.N5.A4RO.01     8  0.0941      0.618 0.00 0.02 0.00 0.02 0.00 0.00 0.00 0.96
#&gt; TCGA.NA.A4QX.01     2  0.0000      0.825 0.00 1.00 0.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RA.01     2  0.3303      0.690 0.00 0.76 0.00 0.20 0.00 0.02 0.02 0.00
#&gt; TCGA.N9.A4Q8.01     1  0.0941      0.742 0.96 0.00 0.00 0.00 0.02 0.02 0.00 0.00
#&gt; TCGA.NG.A4VU.01     2  0.0000      0.825 0.00 1.00 0.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N9.A4Q4.01     2  0.0000      0.825 0.00 1.00 0.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.NF.A4WU.01     1  0.4607      0.250 0.50 0.00 0.00 0.06 0.42 0.02 0.00 0.00
#&gt; TCGA.N8.A4PM.01     3  0.3941      0.658 0.00 0.00 0.68 0.04 0.00 0.00 0.26 0.02
#&gt; TCGA.N9.A4Q7.01     1  0.0941      0.742 0.96 0.00 0.00 0.00 0.02 0.02 0.00 0.00
#&gt; TCGA.NA.A4R1.01     6  0.2348      0.864 0.00 0.00 0.04 0.00 0.00 0.88 0.02 0.06
#&gt; TCGA.N5.A4RJ.01     3  0.4302      0.476 0.00 0.00 0.56 0.00 0.00 0.00 0.36 0.08
#&gt; TCGA.N5.A4RU.01     8  0.0471      0.629 0.00 0.02 0.00 0.00 0.00 0.00 0.00 0.98
#&gt; TCGA.N6.A4V9.01     8  0.4141     -0.117 0.00 0.00 0.06 0.00 0.00 0.00 0.38 0.56
#&gt; TCGA.ND.A4W6.01     5  0.3454      0.769 0.02 0.00 0.00 0.02 0.80 0.04 0.12 0.00
#&gt; TCGA.N6.A4VE.01     3  0.3385      0.648 0.00 0.00 0.76 0.08 0.00 0.00 0.16 0.00
#&gt; TCGA.N7.A4Y8.01     2  0.0808      0.810 0.00 0.96 0.00 0.04 0.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RM.01     8  0.2852      0.266 0.00 0.00 0.00 0.00 0.00 0.00 0.28 0.72
#&gt; TCGA.NF.A4X2.01     6  0.1887      0.864 0.00 0.00 0.04 0.00 0.00 0.90 0.00 0.06
</code></pre>

<script>
$('#tab-node-0-get-classes-7-a').parent().next().next().hide();
$('#tab-node-0-get-classes-7-a').click(function(){
  $('#tab-node-0-get-classes-7-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-node-0-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-node-0-consensus-heatmap'>
<ul>
<li><a href='#tab-node-0-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-node-0-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-node-0-consensus-heatmap-3'>k = 4</a></li>
<li><a href='#tab-node-0-consensus-heatmap-4'>k = 5</a></li>
<li><a href='#tab-node-0-consensus-heatmap-5'>k = 6</a></li>
<li><a href='#tab-node-0-consensus-heatmap-6'>k = 7</a></li>
<li><a href='#tab-node-0-consensus-heatmap-7'>k = 8</a></li>
</ul>
<div id='tab-node-0-consensus-heatmap-1'>
<pre><code class="r">consensus_heatmap(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-0-consensus-heatmap-1-1.png" alt="plot of chunk tab-node-0-consensus-heatmap-1"/></p>

</div>
<div id='tab-node-0-consensus-heatmap-2'>
<pre><code class="r">consensus_heatmap(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-0-consensus-heatmap-2-1.png" alt="plot of chunk tab-node-0-consensus-heatmap-2"/></p>

</div>
<div id='tab-node-0-consensus-heatmap-3'>
<pre><code class="r">consensus_heatmap(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-0-consensus-heatmap-3-1.png" alt="plot of chunk tab-node-0-consensus-heatmap-3"/></p>

</div>
<div id='tab-node-0-consensus-heatmap-4'>
<pre><code class="r">consensus_heatmap(res, k = 5)
</code></pre>

<p><img src="figure_cola/tab-node-0-consensus-heatmap-4-1.png" alt="plot of chunk tab-node-0-consensus-heatmap-4"/></p>

</div>
<div id='tab-node-0-consensus-heatmap-5'>
<pre><code class="r">consensus_heatmap(res, k = 6)
</code></pre>

<p><img src="figure_cola/tab-node-0-consensus-heatmap-5-1.png" alt="plot of chunk tab-node-0-consensus-heatmap-5"/></p>

</div>
<div id='tab-node-0-consensus-heatmap-6'>
<pre><code class="r">consensus_heatmap(res, k = 7)
</code></pre>

<p><img src="figure_cola/tab-node-0-consensus-heatmap-6-1.png" alt="plot of chunk tab-node-0-consensus-heatmap-6"/></p>

</div>
<div id='tab-node-0-consensus-heatmap-7'>
<pre><code class="r">consensus_heatmap(res, k = 8)
</code></pre>

<p><img src="figure_cola/tab-node-0-consensus-heatmap-7-1.png" alt="plot of chunk tab-node-0-consensus-heatmap-7"/></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-node-0-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-node-0-membership-heatmap'>
<ul>
<li><a href='#tab-node-0-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-node-0-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-node-0-membership-heatmap-3'>k = 4</a></li>
<li><a href='#tab-node-0-membership-heatmap-4'>k = 5</a></li>
<li><a href='#tab-node-0-membership-heatmap-5'>k = 6</a></li>
<li><a href='#tab-node-0-membership-heatmap-6'>k = 7</a></li>
<li><a href='#tab-node-0-membership-heatmap-7'>k = 8</a></li>
</ul>
<div id='tab-node-0-membership-heatmap-1'>
<pre><code class="r">membership_heatmap(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-0-membership-heatmap-1-1.png" alt="plot of chunk tab-node-0-membership-heatmap-1"/></p>

</div>
<div id='tab-node-0-membership-heatmap-2'>
<pre><code class="r">membership_heatmap(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-0-membership-heatmap-2-1.png" alt="plot of chunk tab-node-0-membership-heatmap-2"/></p>

</div>
<div id='tab-node-0-membership-heatmap-3'>
<pre><code class="r">membership_heatmap(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-0-membership-heatmap-3-1.png" alt="plot of chunk tab-node-0-membership-heatmap-3"/></p>

</div>
<div id='tab-node-0-membership-heatmap-4'>
<pre><code class="r">membership_heatmap(res, k = 5)
</code></pre>

<p><img src="figure_cola/tab-node-0-membership-heatmap-4-1.png" alt="plot of chunk tab-node-0-membership-heatmap-4"/></p>

</div>
<div id='tab-node-0-membership-heatmap-5'>
<pre><code class="r">membership_heatmap(res, k = 6)
</code></pre>

<p><img src="figure_cola/tab-node-0-membership-heatmap-5-1.png" alt="plot of chunk tab-node-0-membership-heatmap-5"/></p>

</div>
<div id='tab-node-0-membership-heatmap-6'>
<pre><code class="r">membership_heatmap(res, k = 7)
</code></pre>

<p><img src="figure_cola/tab-node-0-membership-heatmap-6-1.png" alt="plot of chunk tab-node-0-membership-heatmap-6"/></p>

</div>
<div id='tab-node-0-membership-heatmap-7'>
<pre><code class="r">membership_heatmap(res, k = 8)
</code></pre>

<p><img src="figure_cola/tab-node-0-membership-heatmap-7-1.png" alt="plot of chunk tab-node-0-membership-heatmap-7"/></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.






<script>
$( function() {
	$( '#tabs-node-0-get-signatures' ).tabs();
} );
</script>
<div id='tabs-node-0-get-signatures'>
<ul>
<li><a href='#tab-node-0-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-node-0-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-node-0-get-signatures-3'>k = 4</a></li>
<li><a href='#tab-node-0-get-signatures-4'>k = 5</a></li>
<li><a href='#tab-node-0-get-signatures-5'>k = 6</a></li>
<li><a href='#tab-node-0-get-signatures-6'>k = 7</a></li>
<li><a href='#tab-node-0-get-signatures-7'>k = 8</a></li>
</ul>
<div id='tab-node-0-get-signatures-1'>
<pre><code class="r">get_signatures(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-0-get-signatures-1-1.png" alt="plot of chunk tab-node-0-get-signatures-1"/></p>

</div>
<div id='tab-node-0-get-signatures-2'>
<pre><code class="r">get_signatures(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-0-get-signatures-2-1.png" alt="plot of chunk tab-node-0-get-signatures-2"/></p>

</div>
<div id='tab-node-0-get-signatures-3'>
<pre><code class="r">get_signatures(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-0-get-signatures-3-1.png" alt="plot of chunk tab-node-0-get-signatures-3"/></p>

</div>
<div id='tab-node-0-get-signatures-4'>
<pre><code class="r">get_signatures(res, k = 5)
</code></pre>

<p><img src="figure_cola/tab-node-0-get-signatures-4-1.png" alt="plot of chunk tab-node-0-get-signatures-4"/></p>

</div>
<div id='tab-node-0-get-signatures-5'>
<pre><code class="r">get_signatures(res, k = 6)
</code></pre>

<p><img src="figure_cola/tab-node-0-get-signatures-5-1.png" alt="plot of chunk tab-node-0-get-signatures-5"/></p>

</div>
<div id='tab-node-0-get-signatures-6'>
<pre><code class="r">get_signatures(res, k = 7)
</code></pre>

<p><img src="figure_cola/tab-node-0-get-signatures-6-1.png" alt="plot of chunk tab-node-0-get-signatures-6"/></p>

</div>
<div id='tab-node-0-get-signatures-7'>
<pre><code class="r">get_signatures(res, k = 8)
</code></pre>

<p><img src="figure_cola/tab-node-0-get-signatures-7-1.png" alt="plot of chunk tab-node-0-get-signatures-7"/></p>

</div>
</div>



Compare the overlap of signatures from different k:

```r
compare_signatures(res)
```

![plot of chunk node-0-signature_compare](figure_cola/node-0-signature_compare-1.png)

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

```r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

```r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```

If the signatures are defined as these which are uniquely high in current group, `diff_method` argument
can be set to `"uniquely_high_in_one_group"`:

```r
# code only for demonstration
tb = get_signature(res, k = ..., diff_method = "uniquely_high_in_one_group")
```




UMAP plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-node-0-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-node-0-dimension-reduction'>
<ul>
<li><a href='#tab-node-0-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-node-0-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-node-0-dimension-reduction-3'>k = 4</a></li>
<li><a href='#tab-node-0-dimension-reduction-4'>k = 5</a></li>
<li><a href='#tab-node-0-dimension-reduction-5'>k = 6</a></li>
<li><a href='#tab-node-0-dimension-reduction-6'>k = 7</a></li>
<li><a href='#tab-node-0-dimension-reduction-7'>k = 8</a></li>
</ul>
<div id='tab-node-0-dimension-reduction-1'>
<pre><code class="r">dimension_reduction(res, k = 2, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-0-dimension-reduction-1-1.png" alt="plot of chunk tab-node-0-dimension-reduction-1"/></p>

</div>
<div id='tab-node-0-dimension-reduction-2'>
<pre><code class="r">dimension_reduction(res, k = 3, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-0-dimension-reduction-2-1.png" alt="plot of chunk tab-node-0-dimension-reduction-2"/></p>

</div>
<div id='tab-node-0-dimension-reduction-3'>
<pre><code class="r">dimension_reduction(res, k = 4, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-0-dimension-reduction-3-1.png" alt="plot of chunk tab-node-0-dimension-reduction-3"/></p>

</div>
<div id='tab-node-0-dimension-reduction-4'>
<pre><code class="r">dimension_reduction(res, k = 5, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-0-dimension-reduction-4-1.png" alt="plot of chunk tab-node-0-dimension-reduction-4"/></p>

</div>
<div id='tab-node-0-dimension-reduction-5'>
<pre><code class="r">dimension_reduction(res, k = 6, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-0-dimension-reduction-5-1.png" alt="plot of chunk tab-node-0-dimension-reduction-5"/></p>

</div>
<div id='tab-node-0-dimension-reduction-6'>
<pre><code class="r">dimension_reduction(res, k = 7, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-0-dimension-reduction-6-1.png" alt="plot of chunk tab-node-0-dimension-reduction-6"/></p>

</div>
<div id='tab-node-0-dimension-reduction-7'>
<pre><code class="r">dimension_reduction(res, k = 8, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-0-dimension-reduction-7-1.png" alt="plot of chunk tab-node-0-dimension-reduction-7"/></p>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

```r
collect_classes(res)
```

![plot of chunk node-0-collect-classes](figure_cola/node-0-collect-classes-1.png)



If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html) for more detailed explanations.


 

---------------------------------------------------




### Node01


Parent node: [Node0](#Node0).
Child nodes: 
                Node011-leaf
        ,
                Node012-leaf
        ,
                Node013-leaf
        ,
                Node021-leaf
        ,
                Node022-leaf
        ,
                Node023-leaf
        .







The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

```r
res = res_rh["01"]
```

A summary of `res` and all the functions that can be applied to it:

```r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4, 5, 6, 7, 8.
#>   On a matrix with 30000 rows and 20 columns.
#>   Top rows (1000) are extracted by 'ATC' method.
#>   Subgroups are detected by 'kmeans' method.
#>   Performed in total 350 partitions by row resampling.
#>   Best k for subgroups seems to be 3.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_partitions"     
#>  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
#> [10] "functional_enrichment"   "get_anno_col"            "get_anno"               
#> [13] "get_classes"             "get_consensus"           "get_matrix"             
#> [16] "get_membership"          "get_param"               "get_signatures"         
#> [19] "get_stats"               "is_best_k"               "is_stable_k"            
#> [22] "membership_heatmap"      "ncol"                    "nrow"                   
#> [25] "plot_ecdf"               "predict_classes"         "rownames"               
#> [28] "select_partition_number" "show"                    "suggest_best_k"         
#> [31] "test_to_known_factors"   "top_rows_heatmap"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

```r
collect_plots(res)
```

![plot of chunk node-01-collect-plots](figure_cola/node-01-collect-plots-1.png)

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

```r
select_partition_number(res)
```

![plot of chunk node-01-select-partition-number](figure_cola/node-01-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

```r
get_stats(res)
```

```
#>   k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2 2 1.000           1.000       1.000         0.5058 0.495   0.495
#> 3 3 1.000           0.979       0.978         0.3715 0.700   0.452
#> 4 4 0.884           0.831       0.925         0.0871 0.937   0.793
#> 5 5 0.805           0.793       0.823         0.0449 1.000   1.000
#> 6 6 0.800           0.476       0.718         0.0475 0.921   0.681
#> 7 7 0.816           0.678       0.807         0.0394 0.884   0.421
#> 8 8 0.832           0.709       0.693         0.0214 0.947   0.615
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

```r
suggest_best_k(res)
```

```
#> [1] 3
#> attr(,"optional")
#> [1] 2
```

There is also optional best $k$ = 2 that is worth to check.

Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-node-01-get-classes' ).tabs();
} );
</script>
<div id='tabs-node-01-get-classes'>
<ul>
<li><a href='#tab-node-01-get-classes-1'>k = 2</a></li>
<li><a href='#tab-node-01-get-classes-2'>k = 3</a></li>
<li><a href='#tab-node-01-get-classes-3'>k = 4</a></li>
<li><a href='#tab-node-01-get-classes-4'>k = 5</a></li>
<li><a href='#tab-node-01-get-classes-5'>k = 6</a></li>
<li><a href='#tab-node-01-get-classes-6'>k = 7</a></li>
<li><a href='#tab-node-01-get-classes-7'>k = 8</a></li>
</ul>

<div id='tab-node-01-get-classes-1'>
<p><a id='tab-node-01-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette p1 p2
#&gt; TCGA.N7.A59B.01     2       0          1  0  1
#&gt; TCGA.NG.A4VW.01     2       0          1  0  1
#&gt; TCGA.N8.A4PO.01     2       0          1  0  1
#&gt; TCGA.ND.A4WA.01     1       0          1  1  0
#&gt; TCGA.N6.A4VC.01     1       0          1  1  0
#&gt; TCGA.N6.A4VG.01     1       0          1  1  0
#&gt; TCGA.N9.A4Q3.01     1       0          1  1  0
#&gt; TCGA.ND.A4WC.01     2       0          1  0  1
#&gt; TCGA.NA.A5I1.01     1       0          1  1  0
#&gt; TCGA.N8.A56S.01     1       0          1  1  0
#&gt; TCGA.N8.A4PL.01     1       0          1  1  0
#&gt; TCGA.NA.A4R0.01     2       0          1  0  1
#&gt; TCGA.NA.A4QY.01     1       0          1  1  0
#&gt; TCGA.NF.A4WX.01     1       0          1  1  0
#&gt; TCGA.N8.A4PP.01     2       0          1  0  1
#&gt; TCGA.N5.A4RS.01     1       0          1  1  0
#&gt; TCGA.N9.A4Q8.01     2       0          1  0  1
#&gt; TCGA.NF.A4WU.01     1       0          1  1  0
#&gt; TCGA.N9.A4Q7.01     2       0          1  0  1
#&gt; TCGA.ND.A4W6.01     1       0          1  1  0
</code></pre>

<script>
$('#tab-node-01-get-classes-1-a').parent().next().next().hide();
$('#tab-node-01-get-classes-1-a').click(function(){
  $('#tab-node-01-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-01-get-classes-2'>
<p><a id='tab-node-01-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1   p2   p3
#&gt; TCGA.N7.A59B.01     2   0.000      0.981 0.00 1.00 0.00
#&gt; TCGA.NG.A4VW.01     2   0.207      0.961 0.00 0.94 0.06
#&gt; TCGA.N8.A4PO.01     3   0.000      0.940 0.00 0.00 1.00
#&gt; TCGA.ND.A4WA.01     1   0.000      1.000 1.00 0.00 0.00
#&gt; TCGA.N6.A4VC.01     3   0.207      0.974 0.06 0.00 0.94
#&gt; TCGA.N6.A4VG.01     3   0.207      0.974 0.06 0.00 0.94
#&gt; TCGA.N9.A4Q3.01     1   0.000      1.000 1.00 0.00 0.00
#&gt; TCGA.ND.A4WC.01     2   0.000      0.981 0.00 1.00 0.00
#&gt; TCGA.NA.A5I1.01     1   0.000      1.000 1.00 0.00 0.00
#&gt; TCGA.N8.A56S.01     1   0.000      1.000 1.00 0.00 0.00
#&gt; TCGA.N8.A4PL.01     1   0.000      1.000 1.00 0.00 0.00
#&gt; TCGA.NA.A4R0.01     3   0.000      0.940 0.00 0.00 1.00
#&gt; TCGA.NA.A4QY.01     3   0.254      0.959 0.08 0.00 0.92
#&gt; TCGA.NF.A4WX.01     1   0.000      1.000 1.00 0.00 0.00
#&gt; TCGA.N8.A4PP.01     2   0.207      0.961 0.00 0.94 0.06
#&gt; TCGA.N5.A4RS.01     1   0.000      1.000 1.00 0.00 0.00
#&gt; TCGA.N9.A4Q8.01     2   0.000      0.981 0.00 1.00 0.00
#&gt; TCGA.NF.A4WU.01     3   0.207      0.974 0.06 0.00 0.94
#&gt; TCGA.N9.A4Q7.01     2   0.000      0.981 0.00 1.00 0.00
#&gt; TCGA.ND.A4W6.01     3   0.207      0.974 0.06 0.00 0.94
</code></pre>

<script>
$('#tab-node-01-get-classes-2-a').parent().next().next().hide();
$('#tab-node-01-get-classes-2-a').click(function(){
  $('#tab-node-01-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-01-get-classes-3'>
<p><a id='tab-node-01-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1   p2   p3   p4
#&gt; TCGA.N7.A59B.01     2  0.4406      0.765 0.00 0.70 0.00 0.30
#&gt; TCGA.NG.A4VW.01     2  0.0000      0.647 0.00 1.00 0.00 0.00
#&gt; TCGA.N8.A4PO.01     4  0.5392      0.494 0.00 0.28 0.04 0.68
#&gt; TCGA.ND.A4WA.01     1  0.0707      0.988 0.98 0.00 0.00 0.02
#&gt; TCGA.N6.A4VC.01     3  0.0707      0.932 0.02 0.00 0.98 0.00
#&gt; TCGA.N6.A4VG.01     3  0.0707      0.932 0.02 0.00 0.98 0.00
#&gt; TCGA.N9.A4Q3.01     1  0.0000      0.991 1.00 0.00 0.00 0.00
#&gt; TCGA.ND.A4WC.01     2  0.4406      0.765 0.00 0.70 0.00 0.30
#&gt; TCGA.NA.A5I1.01     1  0.0000      0.991 1.00 0.00 0.00 0.00
#&gt; TCGA.N8.A56S.01     1  0.0000      0.991 1.00 0.00 0.00 0.00
#&gt; TCGA.N8.A4PL.01     1  0.0707      0.988 0.98 0.00 0.00 0.02
#&gt; TCGA.NA.A4R0.01     3  0.4406      0.570 0.00 0.30 0.70 0.00
#&gt; TCGA.NA.A4QY.01     3  0.0707      0.932 0.02 0.00 0.98 0.00
#&gt; TCGA.NF.A4WX.01     1  0.0707      0.988 0.98 0.00 0.00 0.02
#&gt; TCGA.N8.A4PP.01     2  0.0000      0.647 0.00 1.00 0.00 0.00
#&gt; TCGA.N5.A4RS.01     1  0.0000      0.991 1.00 0.00 0.00 0.00
#&gt; TCGA.N9.A4Q8.01     4  0.2647      0.368 0.00 0.12 0.00 0.88
#&gt; TCGA.NF.A4WU.01     3  0.0707      0.932 0.02 0.00 0.98 0.00
#&gt; TCGA.N9.A4Q7.01     2  0.4277      0.769 0.00 0.72 0.00 0.28
#&gt; TCGA.ND.A4W6.01     3  0.0707      0.932 0.02 0.00 0.98 0.00
</code></pre>

<script>
$('#tab-node-01-get-classes-3-a').parent().next().next().hide();
$('#tab-node-01-get-classes-3-a').click(function(){
  $('#tab-node-01-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-01-get-classes-4'>
<p><a id='tab-node-01-get-classes-4-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 5), get_membership(res, k = 5))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1   p2   p3   p4   p5
#&gt; TCGA.N7.A59B.01     2  0.3274      0.773 0.00 0.78 0.00 0.00 0.22
#&gt; TCGA.NG.A4VW.01     2  0.1732      0.816 0.00 0.92 0.00 0.08 0.00
#&gt; TCGA.N8.A4PO.01     4  0.0000      0.737 0.00 0.00 0.00 1.00 0.00
#&gt; TCGA.ND.A4WA.01     1  0.3796      0.800 0.70 0.00 0.00 0.00 0.30
#&gt; TCGA.N6.A4VC.01     3  0.4527      0.622 0.04 0.00 0.70 0.00 0.26
#&gt; TCGA.N6.A4VG.01     3  0.1043      0.812 0.04 0.00 0.96 0.00 0.00
#&gt; TCGA.N9.A4Q3.01     1  0.0000      0.850 1.00 0.00 0.00 0.00 0.00
#&gt; TCGA.ND.A4WC.01     2  0.3274      0.773 0.00 0.78 0.00 0.00 0.22
#&gt; TCGA.NA.A5I1.01     1  0.0609      0.850 0.98 0.00 0.00 0.00 0.02
#&gt; TCGA.N8.A56S.01     1  0.0000      0.850 1.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N8.A4PL.01     1  0.3796      0.800 0.70 0.00 0.00 0.00 0.30
#&gt; TCGA.NA.A4R0.01     3  0.5685      0.734 0.00 0.04 0.68 0.08 0.20
#&gt; TCGA.NA.A4QY.01     3  0.1043      0.812 0.04 0.00 0.96 0.00 0.00
#&gt; TCGA.NF.A4WX.01     1  0.3796      0.800 0.70 0.00 0.00 0.00 0.30
#&gt; TCGA.N8.A4PP.01     2  0.1732      0.816 0.00 0.92 0.00 0.08 0.00
#&gt; TCGA.N5.A4RS.01     1  0.0000      0.850 1.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N9.A4Q8.01     4  0.5604      0.703 0.00 0.16 0.04 0.70 0.10
#&gt; TCGA.NF.A4WU.01     3  0.4096      0.815 0.04 0.00 0.76 0.00 0.20
#&gt; TCGA.N9.A4Q7.01     2  0.0000      0.827 0.00 1.00 0.00 0.00 0.00
#&gt; TCGA.ND.A4W6.01     3  0.4096      0.815 0.04 0.00 0.76 0.00 0.20
</code></pre>

<script>
$('#tab-node-01-get-classes-4-a').parent().next().next().hide();
$('#tab-node-01-get-classes-4-a').click(function(){
  $('#tab-node-01-get-classes-4-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-01-get-classes-5'>
<p><a id='tab-node-01-get-classes-5-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 6), get_membership(res, k = 6))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1   p2   p3   p4   p5   p6
#&gt; TCGA.N7.A59B.01     2   0.000      0.638 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; TCGA.NG.A4VW.01     2   0.425      0.754 0.00 0.58 0.00 0.02 0.40 0.00
#&gt; TCGA.N8.A4PO.01     4   0.371      0.589 0.00 0.00 0.00 0.72 0.02 0.26
#&gt; TCGA.ND.A4WA.01     1   0.148      0.481 0.94 0.00 0.00 0.00 0.04 0.02
#&gt; TCGA.N6.A4VC.01     5   0.558      0.000 0.00 0.00 0.40 0.00 0.46 0.14
#&gt; TCGA.N6.A4VG.01     3   0.000      0.463 0.00 0.00 1.00 0.00 0.00 0.00
#&gt; TCGA.N9.A4Q3.01     1   0.387     -1.000 0.50 0.00 0.00 0.00 0.00 0.50
#&gt; TCGA.ND.A4WC.01     2   0.000      0.638 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; TCGA.NA.A5I1.01     1   0.376     -0.711 0.60 0.00 0.00 0.00 0.00 0.40
#&gt; TCGA.N8.A56S.01     6   0.387      1.000 0.50 0.00 0.00 0.00 0.00 0.50
#&gt; TCGA.N8.A4PL.01     1   0.000      0.510 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.NA.A4R0.01     3   0.520      0.669 0.00 0.00 0.58 0.30 0.12 0.00
#&gt; TCGA.NA.A4QY.01     3   0.000      0.463 0.00 0.00 1.00 0.00 0.00 0.00
#&gt; TCGA.NF.A4WX.01     1   0.000      0.510 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N8.A4PP.01     2   0.425      0.754 0.00 0.58 0.00 0.02 0.40 0.00
#&gt; TCGA.N5.A4RS.01     6   0.387      1.000 0.50 0.00 0.00 0.00 0.00 0.50
#&gt; TCGA.N9.A4Q8.01     4   0.649      0.601 0.00 0.28 0.00 0.40 0.02 0.30
#&gt; TCGA.NF.A4WU.01     3   0.493      0.700 0.00 0.00 0.62 0.28 0.10 0.00
#&gt; TCGA.N9.A4Q7.01     2   0.365      0.760 0.00 0.64 0.00 0.00 0.36 0.00
#&gt; TCGA.ND.A4W6.01     3   0.493      0.700 0.00 0.00 0.62 0.28 0.10 0.00
</code></pre>

<script>
$('#tab-node-01-get-classes-5-a').parent().next().next().hide();
$('#tab-node-01-get-classes-5-a').click(function(){
  $('#tab-node-01-get-classes-5-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-01-get-classes-6'>
<p><a id='tab-node-01-get-classes-6-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 7), get_membership(res, k = 7))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1   p2   p3   p4   p5   p6   p7
#&gt; TCGA.N7.A59B.01     2  0.3562    0.00743 0.00 0.50 0.00 0.50 0.00 0.00 0.00
#&gt; TCGA.NG.A4VW.01     2  0.0000    0.78001 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N8.A4PO.01     4  0.5282    0.09077 0.18 0.00 0.00 0.48 0.00 0.00 0.34
#&gt; TCGA.ND.A4WA.01     1  0.4461    0.90149 0.62 0.00 0.00 0.00 0.02 0.32 0.04
#&gt; TCGA.N6.A4VC.01     5  0.1664    0.00000 0.00 0.00 0.06 0.00 0.92 0.00 0.02
#&gt; TCGA.N6.A4VG.01     7  0.4418    1.00000 0.00 0.00 0.30 0.00 0.08 0.00 0.62
#&gt; TCGA.N9.A4Q3.01     6  0.0000    0.93947 0.00 0.00 0.00 0.00 0.00 1.00 0.00
#&gt; TCGA.ND.A4WC.01     4  0.3562   -0.50372 0.00 0.50 0.00 0.50 0.00 0.00 0.00
#&gt; TCGA.NA.A5I1.01     6  0.2081    0.79612 0.14 0.00 0.00 0.00 0.00 0.86 0.00
#&gt; TCGA.N8.A56S.01     6  0.0000    0.93947 0.00 0.00 0.00 0.00 0.00 1.00 0.00
#&gt; TCGA.N8.A4PL.01     1  0.3139    0.95155 0.70 0.00 0.00 0.00 0.00 0.30 0.00
#&gt; TCGA.NA.A4R0.01     3  0.0504    0.96202 0.00 0.02 0.98 0.00 0.00 0.00 0.00
#&gt; TCGA.NA.A4QY.01     7  0.4418    1.00000 0.00 0.00 0.30 0.00 0.08 0.00 0.62
#&gt; TCGA.NF.A4WX.01     1  0.3139    0.95155 0.70 0.00 0.00 0.00 0.00 0.30 0.00
#&gt; TCGA.N8.A4PP.01     2  0.0000    0.78001 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RS.01     6  0.0000    0.93947 0.00 0.00 0.00 0.00 0.00 1.00 0.00
#&gt; TCGA.N9.A4Q8.01     4  0.6534    0.27993 0.16 0.12 0.00 0.56 0.06 0.00 0.10
#&gt; TCGA.NF.A4WU.01     3  0.0000    0.98087 0.00 0.00 1.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N9.A4Q7.01     2  0.0000    0.78001 0.00 1.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.ND.A4W6.01     3  0.0000    0.98087 0.00 0.00 1.00 0.00 0.00 0.00 0.00
</code></pre>

<script>
$('#tab-node-01-get-classes-6-a').parent().next().next().hide();
$('#tab-node-01-get-classes-6-a').click(function(){
  $('#tab-node-01-get-classes-6-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-01-get-classes-7'>
<p><a id='tab-node-01-get-classes-7-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 8), get_membership(res, k = 8))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1   p2   p3   p4   p5   p6   p7   p8
#&gt; TCGA.N7.A59B.01     4  0.2406      1.000 0.00 0.20 0.00 0.80 0.00 0.00 0.00 0.00
#&gt; TCGA.NG.A4VW.01     2  0.0471      0.958 0.00 0.98 0.00 0.00 0.02 0.00 0.00 0.00
#&gt; TCGA.N8.A4PO.01     8  0.0000      0.000 0.00 0.00 0.00 0.00 0.00 0.00 0.00 1.00
#&gt; TCGA.ND.A4WA.01     1  0.3808      0.510 0.62 0.00 0.00 0.04 0.00 0.34 0.00 0.00
#&gt; TCGA.N6.A4VC.01     5  0.1563      0.000 0.00 0.00 0.00 0.00 0.90 0.00 0.10 0.00
#&gt; TCGA.N6.A4VG.01     7  0.0000      1.000 0.00 0.00 0.00 0.00 0.00 0.00 1.00 0.00
#&gt; TCGA.N9.A4Q3.01     6  0.0000      0.907 0.00 0.00 0.00 0.00 0.00 1.00 0.00 0.00
#&gt; TCGA.ND.A4WC.01     4  0.2406      1.000 0.00 0.20 0.00 0.80 0.00 0.00 0.00 0.00
#&gt; TCGA.NA.A5I1.01     6  0.3434      0.676 0.14 0.00 0.00 0.10 0.00 0.76 0.00 0.00
#&gt; TCGA.N8.A56S.01     6  0.0000      0.907 0.00 0.00 0.00 0.00 0.00 1.00 0.00 0.00
#&gt; TCGA.N8.A4PL.01     1  0.3083      0.518 0.66 0.00 0.00 0.00 0.00 0.34 0.00 0.00
#&gt; TCGA.NA.A4R0.01     3  0.3036      0.940 0.00 0.04 0.78 0.00 0.00 0.00 0.18 0.00
#&gt; TCGA.NA.A4QY.01     7  0.0000      1.000 0.00 0.00 0.00 0.00 0.00 0.00 1.00 0.00
#&gt; TCGA.NF.A4WX.01     1  0.3142      0.498 0.64 0.00 0.00 0.00 0.00 0.36 0.00 0.00
#&gt; TCGA.N8.A4PP.01     2  0.0000      0.962 0.00 1.00 0.00 0.00 0.00 0.00 0.00 0.00
#&gt; TCGA.N5.A4RS.01     6  0.0000      0.907 0.00 0.00 0.00 0.00 0.00 1.00 0.00 0.00
#&gt; TCGA.N9.A4Q8.01     1  0.7339     -0.472 0.32 0.00 0.22 0.16 0.08 0.00 0.00 0.22
#&gt; TCGA.NF.A4WU.01     3  0.2534      0.970 0.00 0.00 0.78 0.00 0.00 0.00 0.22 0.00
#&gt; TCGA.N9.A4Q7.01     2  0.0808      0.938 0.00 0.96 0.00 0.04 0.00 0.00 0.00 0.00
#&gt; TCGA.ND.A4W6.01     3  0.2534      0.970 0.00 0.00 0.78 0.00 0.00 0.00 0.22 0.00
</code></pre>

<script>
$('#tab-node-01-get-classes-7-a').parent().next().next().hide();
$('#tab-node-01-get-classes-7-a').click(function(){
  $('#tab-node-01-get-classes-7-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-node-01-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-node-01-consensus-heatmap'>
<ul>
<li><a href='#tab-node-01-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-node-01-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-node-01-consensus-heatmap-3'>k = 4</a></li>
<li><a href='#tab-node-01-consensus-heatmap-4'>k = 5</a></li>
<li><a href='#tab-node-01-consensus-heatmap-5'>k = 6</a></li>
<li><a href='#tab-node-01-consensus-heatmap-6'>k = 7</a></li>
<li><a href='#tab-node-01-consensus-heatmap-7'>k = 8</a></li>
</ul>
<div id='tab-node-01-consensus-heatmap-1'>
<pre><code class="r">consensus_heatmap(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-01-consensus-heatmap-1-1.png" alt="plot of chunk tab-node-01-consensus-heatmap-1"/></p>

</div>
<div id='tab-node-01-consensus-heatmap-2'>
<pre><code class="r">consensus_heatmap(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-01-consensus-heatmap-2-1.png" alt="plot of chunk tab-node-01-consensus-heatmap-2"/></p>

</div>
<div id='tab-node-01-consensus-heatmap-3'>
<pre><code class="r">consensus_heatmap(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-01-consensus-heatmap-3-1.png" alt="plot of chunk tab-node-01-consensus-heatmap-3"/></p>

</div>
<div id='tab-node-01-consensus-heatmap-4'>
<pre><code class="r">consensus_heatmap(res, k = 5)
</code></pre>

<p><img src="figure_cola/tab-node-01-consensus-heatmap-4-1.png" alt="plot of chunk tab-node-01-consensus-heatmap-4"/></p>

</div>
<div id='tab-node-01-consensus-heatmap-5'>
<pre><code class="r">consensus_heatmap(res, k = 6)
</code></pre>

<p><img src="figure_cola/tab-node-01-consensus-heatmap-5-1.png" alt="plot of chunk tab-node-01-consensus-heatmap-5"/></p>

</div>
<div id='tab-node-01-consensus-heatmap-6'>
<pre><code class="r">consensus_heatmap(res, k = 7)
</code></pre>

<p><img src="figure_cola/tab-node-01-consensus-heatmap-6-1.png" alt="plot of chunk tab-node-01-consensus-heatmap-6"/></p>

</div>
<div id='tab-node-01-consensus-heatmap-7'>
<pre><code class="r">consensus_heatmap(res, k = 8)
</code></pre>

<p><img src="figure_cola/tab-node-01-consensus-heatmap-7-1.png" alt="plot of chunk tab-node-01-consensus-heatmap-7"/></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-node-01-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-node-01-membership-heatmap'>
<ul>
<li><a href='#tab-node-01-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-node-01-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-node-01-membership-heatmap-3'>k = 4</a></li>
<li><a href='#tab-node-01-membership-heatmap-4'>k = 5</a></li>
<li><a href='#tab-node-01-membership-heatmap-5'>k = 6</a></li>
<li><a href='#tab-node-01-membership-heatmap-6'>k = 7</a></li>
<li><a href='#tab-node-01-membership-heatmap-7'>k = 8</a></li>
</ul>
<div id='tab-node-01-membership-heatmap-1'>
<pre><code class="r">membership_heatmap(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-01-membership-heatmap-1-1.png" alt="plot of chunk tab-node-01-membership-heatmap-1"/></p>

</div>
<div id='tab-node-01-membership-heatmap-2'>
<pre><code class="r">membership_heatmap(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-01-membership-heatmap-2-1.png" alt="plot of chunk tab-node-01-membership-heatmap-2"/></p>

</div>
<div id='tab-node-01-membership-heatmap-3'>
<pre><code class="r">membership_heatmap(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-01-membership-heatmap-3-1.png" alt="plot of chunk tab-node-01-membership-heatmap-3"/></p>

</div>
<div id='tab-node-01-membership-heatmap-4'>
<pre><code class="r">membership_heatmap(res, k = 5)
</code></pre>

<p><img src="figure_cola/tab-node-01-membership-heatmap-4-1.png" alt="plot of chunk tab-node-01-membership-heatmap-4"/></p>

</div>
<div id='tab-node-01-membership-heatmap-5'>
<pre><code class="r">membership_heatmap(res, k = 6)
</code></pre>

<p><img src="figure_cola/tab-node-01-membership-heatmap-5-1.png" alt="plot of chunk tab-node-01-membership-heatmap-5"/></p>

</div>
<div id='tab-node-01-membership-heatmap-6'>
<pre><code class="r">membership_heatmap(res, k = 7)
</code></pre>

<p><img src="figure_cola/tab-node-01-membership-heatmap-6-1.png" alt="plot of chunk tab-node-01-membership-heatmap-6"/></p>

</div>
<div id='tab-node-01-membership-heatmap-7'>
<pre><code class="r">membership_heatmap(res, k = 8)
</code></pre>

<p><img src="figure_cola/tab-node-01-membership-heatmap-7-1.png" alt="plot of chunk tab-node-01-membership-heatmap-7"/></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.






<script>
$( function() {
	$( '#tabs-node-01-get-signatures' ).tabs();
} );
</script>
<div id='tabs-node-01-get-signatures'>
<ul>
<li><a href='#tab-node-01-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-node-01-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-node-01-get-signatures-3'>k = 4</a></li>
<li><a href='#tab-node-01-get-signatures-4'>k = 5</a></li>
<li><a href='#tab-node-01-get-signatures-5'>k = 6</a></li>
<li><a href='#tab-node-01-get-signatures-6'>k = 7</a></li>
<li><a href='#tab-node-01-get-signatures-7'>k = 8</a></li>
</ul>
<div id='tab-node-01-get-signatures-1'>
<pre><code class="r">get_signatures(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-01-get-signatures-1-1.png" alt="plot of chunk tab-node-01-get-signatures-1"/></p>

</div>
<div id='tab-node-01-get-signatures-2'>
<pre><code class="r">get_signatures(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-01-get-signatures-2-1.png" alt="plot of chunk tab-node-01-get-signatures-2"/></p>

</div>
<div id='tab-node-01-get-signatures-3'>
<pre><code class="r">get_signatures(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-01-get-signatures-3-1.png" alt="plot of chunk tab-node-01-get-signatures-3"/></p>

</div>
<div id='tab-node-01-get-signatures-4'>
<pre><code class="r">get_signatures(res, k = 5)
</code></pre>

<p><img src="figure_cola/tab-node-01-get-signatures-4-1.png" alt="plot of chunk tab-node-01-get-signatures-4"/></p>

</div>
<div id='tab-node-01-get-signatures-5'>
<pre><code class="r">get_signatures(res, k = 6)
</code></pre>

<p><img src="figure_cola/tab-node-01-get-signatures-5-1.png" alt="plot of chunk tab-node-01-get-signatures-5"/></p>

</div>
<div id='tab-node-01-get-signatures-6'>
<pre><code class="r">get_signatures(res, k = 7)
</code></pre>

<p><img src="figure_cola/tab-node-01-get-signatures-6-1.png" alt="plot of chunk tab-node-01-get-signatures-6"/></p>

</div>
<div id='tab-node-01-get-signatures-7'>
<pre><code class="r">get_signatures(res, k = 8)
</code></pre>

<p><img src="figure_cola/tab-node-01-get-signatures-7-1.png" alt="plot of chunk tab-node-01-get-signatures-7"/></p>

</div>
</div>



Compare the overlap of signatures from different k:

```r
compare_signatures(res)
```

![plot of chunk node-01-signature_compare](figure_cola/node-01-signature_compare-1.png)

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

```r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

```r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```

If the signatures are defined as these which are uniquely high in current group, `diff_method` argument
can be set to `"uniquely_high_in_one_group"`:

```r
# code only for demonstration
tb = get_signature(res, k = ..., diff_method = "uniquely_high_in_one_group")
```




UMAP plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-node-01-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-node-01-dimension-reduction'>
<ul>
<li><a href='#tab-node-01-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-node-01-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-node-01-dimension-reduction-3'>k = 4</a></li>
<li><a href='#tab-node-01-dimension-reduction-4'>k = 5</a></li>
<li><a href='#tab-node-01-dimension-reduction-5'>k = 6</a></li>
<li><a href='#tab-node-01-dimension-reduction-6'>k = 7</a></li>
<li><a href='#tab-node-01-dimension-reduction-7'>k = 8</a></li>
</ul>
<div id='tab-node-01-dimension-reduction-1'>
<pre><code class="r">dimension_reduction(res, k = 2, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-01-dimension-reduction-1-1.png" alt="plot of chunk tab-node-01-dimension-reduction-1"/></p>

</div>
<div id='tab-node-01-dimension-reduction-2'>
<pre><code class="r">dimension_reduction(res, k = 3, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-01-dimension-reduction-2-1.png" alt="plot of chunk tab-node-01-dimension-reduction-2"/></p>

</div>
<div id='tab-node-01-dimension-reduction-3'>
<pre><code class="r">dimension_reduction(res, k = 4, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-01-dimension-reduction-3-1.png" alt="plot of chunk tab-node-01-dimension-reduction-3"/></p>

</div>
<div id='tab-node-01-dimension-reduction-4'>
<pre><code class="r">dimension_reduction(res, k = 5, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-01-dimension-reduction-4-1.png" alt="plot of chunk tab-node-01-dimension-reduction-4"/></p>

</div>
<div id='tab-node-01-dimension-reduction-5'>
<pre><code class="r">dimension_reduction(res, k = 6, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-01-dimension-reduction-5-1.png" alt="plot of chunk tab-node-01-dimension-reduction-5"/></p>

</div>
<div id='tab-node-01-dimension-reduction-6'>
<pre><code class="r">dimension_reduction(res, k = 7, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-01-dimension-reduction-6-1.png" alt="plot of chunk tab-node-01-dimension-reduction-6"/></p>

</div>
<div id='tab-node-01-dimension-reduction-7'>
<pre><code class="r">dimension_reduction(res, k = 8, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-01-dimension-reduction-7-1.png" alt="plot of chunk tab-node-01-dimension-reduction-7"/></p>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

```r
collect_classes(res)
```

![plot of chunk node-01-collect-classes](figure_cola/node-01-collect-classes-1.png)



If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html) for more detailed explanations.


 

---------------------------------------------------




### Node02


Parent node: [Node0](#Node0).
Child nodes: 
                Node011-leaf
        ,
                Node012-leaf
        ,
                Node013-leaf
        ,
                Node021-leaf
        ,
                Node022-leaf
        ,
                Node023-leaf
        .







The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

```r
res = res_rh["02"]
```

A summary of `res` and all the functions that can be applied to it:

```r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4, 5, 6, 7, 8.
#>   On a matrix with 30000 rows and 18 columns.
#>   Top rows (1000) are extracted by 'ATC' method.
#>   Subgroups are detected by 'skmeans' method.
#>   Performed in total 350 partitions by row resampling.
#>   Best k for subgroups seems to be 7.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_partitions"     
#>  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
#> [10] "functional_enrichment"   "get_anno_col"            "get_anno"               
#> [13] "get_classes"             "get_consensus"           "get_matrix"             
#> [16] "get_membership"          "get_param"               "get_signatures"         
#> [19] "get_stats"               "is_best_k"               "is_stable_k"            
#> [22] "membership_heatmap"      "ncol"                    "nrow"                   
#> [25] "plot_ecdf"               "predict_classes"         "rownames"               
#> [28] "select_partition_number" "show"                    "suggest_best_k"         
#> [31] "test_to_known_factors"   "top_rows_heatmap"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

```r
collect_plots(res)
```

![plot of chunk node-02-collect-plots](figure_cola/node-02-collect-plots-1.png)

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

```r
select_partition_number(res)
```

![plot of chunk node-02-select-partition-number](figure_cola/node-02-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

```r
get_stats(res)
```

```
#>   k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2 2 1.000           1.000       1.000         0.5038 0.497   0.497
#> 3 3 0.941           0.959       0.976         0.2851 0.765   0.566
#> 4 4 0.902           0.925       0.963         0.1977 0.784   0.450
#> 5 5 0.908           0.789       0.948         0.0364 0.974   0.879
#> 6 6 0.948           0.859       0.979         0.0287 0.987   0.931
#> 7 7 0.954           0.720       0.948         0.0368 0.974   0.852
#> 8 8 0.980           0.759       0.991         0.0155 0.980   0.870
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

```r
suggest_best_k(res)
```

```
#> [1] 7
#> attr(,"optional")
#> [1] 2 3 4
```

There is also optional best $k$ = 2 3 4 that is worth to check.

Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-node-02-get-classes' ).tabs();
} );
</script>
<div id='tabs-node-02-get-classes'>
<ul>
<li><a href='#tab-node-02-get-classes-1'>k = 2</a></li>
<li><a href='#tab-node-02-get-classes-2'>k = 3</a></li>
<li><a href='#tab-node-02-get-classes-3'>k = 4</a></li>
<li><a href='#tab-node-02-get-classes-4'>k = 5</a></li>
<li><a href='#tab-node-02-get-classes-5'>k = 6</a></li>
<li><a href='#tab-node-02-get-classes-6'>k = 7</a></li>
<li><a href='#tab-node-02-get-classes-7'>k = 8</a></li>
</ul>

<div id='tab-node-02-get-classes-1'>
<p><a id='tab-node-02-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette p1 p2
#&gt; TCGA.N8.A4PN.01     2       0          1  0  1
#&gt; TCGA.N5.A59E.01     1       0          1  1  0
#&gt; TCGA.N7.A4Y5.01     2       0          1  0  1
#&gt; TCGA.N5.A4RD.01     1       0          1  1  0
#&gt; TCGA.N8.A4PQ.01     2       0          1  0  1
#&gt; TCGA.ND.A4WF.01     1       0          1  1  0
#&gt; TCGA.N5.A4R8.01     2       0          1  0  1
#&gt; TCGA.NA.A4QW.01     2       0          1  0  1
#&gt; TCGA.N9.A4Q1.01     2       0          1  0  1
#&gt; TCGA.N5.A4RV.01     1       0          1  1  0
#&gt; TCGA.N5.A4RO.01     1       0          1  1  0
#&gt; TCGA.NA.A4QX.01     2       0          1  0  1
#&gt; TCGA.N5.A4RA.01     2       0          1  0  1
#&gt; TCGA.NG.A4VU.01     2       0          1  0  1
#&gt; TCGA.N9.A4Q4.01     2       0          1  0  1
#&gt; TCGA.N5.A4RU.01     1       0          1  1  0
#&gt; TCGA.N7.A4Y8.01     2       0          1  0  1
#&gt; TCGA.N5.A4RM.01     1       0          1  1  0
</code></pre>

<script>
$('#tab-node-02-get-classes-1-a').parent().next().next().hide();
$('#tab-node-02-get-classes-1-a').click(function(){
  $('#tab-node-02-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-02-get-classes-2'>
<p><a id='tab-node-02-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1   p2   p3
#&gt; TCGA.N8.A4PN.01     2   0.153      0.974 0.00 0.96 0.04
#&gt; TCGA.N5.A59E.01     1   0.000      0.974 1.00 0.00 0.00
#&gt; TCGA.N7.A4Y5.01     2   0.000      0.980 0.00 1.00 0.00
#&gt; TCGA.N5.A4RD.01     1   0.000      0.974 1.00 0.00 0.00
#&gt; TCGA.N8.A4PQ.01     2   0.153      0.974 0.00 0.96 0.04
#&gt; TCGA.ND.A4WF.01     1   0.254      0.912 0.92 0.00 0.08
#&gt; TCGA.N5.A4R8.01     3   0.254      0.861 0.00 0.08 0.92
#&gt; TCGA.NA.A4QW.01     2   0.000      0.980 0.00 1.00 0.00
#&gt; TCGA.N9.A4Q1.01     2   0.000      0.980 0.00 1.00 0.00
#&gt; TCGA.N5.A4RV.01     3   0.153      0.929 0.04 0.00 0.96
#&gt; TCGA.N5.A4RO.01     1   0.000      0.974 1.00 0.00 0.00
#&gt; TCGA.NA.A4QX.01     2   0.153      0.974 0.00 0.96 0.04
#&gt; TCGA.N5.A4RA.01     2   0.153      0.974 0.00 0.96 0.04
#&gt; TCGA.NG.A4VU.01     2   0.000      0.980 0.00 1.00 0.00
#&gt; TCGA.N9.A4Q4.01     2   0.000      0.980 0.00 1.00 0.00
#&gt; TCGA.N5.A4RU.01     3   0.153      0.929 0.04 0.00 0.96
#&gt; TCGA.N7.A4Y8.01     1   0.153      0.939 0.96 0.00 0.04
#&gt; TCGA.N5.A4RM.01     1   0.000      0.974 1.00 0.00 0.00
</code></pre>

<script>
$('#tab-node-02-get-classes-2-a').parent().next().next().hide();
$('#tab-node-02-get-classes-2-a').click(function(){
  $('#tab-node-02-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-02-get-classes-3'>
<p><a id='tab-node-02-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1   p2   p3   p4
#&gt; TCGA.N8.A4PN.01     4   0.121      0.988 0.00 0.04 0.00 0.96
#&gt; TCGA.N5.A59E.01     1   0.000      0.918 1.00 0.00 0.00 0.00
#&gt; TCGA.N7.A4Y5.01     2   0.000      1.000 0.00 1.00 0.00 0.00
#&gt; TCGA.N5.A4RD.01     1   0.000      0.918 1.00 0.00 0.00 0.00
#&gt; TCGA.N8.A4PQ.01     4   0.121      0.988 0.00 0.04 0.00 0.96
#&gt; TCGA.ND.A4WF.01     3   0.413      0.671 0.26 0.00 0.74 0.00
#&gt; TCGA.N5.A4R8.01     4   0.000      0.953 0.00 0.00 0.00 1.00
#&gt; TCGA.NA.A4QW.01     2   0.000      1.000 0.00 1.00 0.00 0.00
#&gt; TCGA.N9.A4Q1.01     2   0.000      1.000 0.00 1.00 0.00 0.00
#&gt; TCGA.N5.A4RV.01     3   0.000      0.867 0.00 0.00 1.00 0.00
#&gt; TCGA.N5.A4RO.01     1   0.000      0.918 1.00 0.00 0.00 0.00
#&gt; TCGA.NA.A4QX.01     4   0.121      0.988 0.00 0.04 0.00 0.96
#&gt; TCGA.N5.A4RA.01     4   0.121      0.988 0.00 0.04 0.00 0.96
#&gt; TCGA.NG.A4VU.01     2   0.000      1.000 0.00 1.00 0.00 0.00
#&gt; TCGA.N9.A4Q4.01     2   0.000      1.000 0.00 1.00 0.00 0.00
#&gt; TCGA.N5.A4RU.01     3   0.000      0.867 0.00 0.00 1.00 0.00
#&gt; TCGA.N7.A4Y8.01     1   0.398      0.660 0.76 0.00 0.00 0.24
#&gt; TCGA.N5.A4RM.01     1   0.000      0.918 1.00 0.00 0.00 0.00
</code></pre>

<script>
$('#tab-node-02-get-classes-3-a').parent().next().next().hide();
$('#tab-node-02-get-classes-3-a').click(function(){
  $('#tab-node-02-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-02-get-classes-4'>
<p><a id='tab-node-02-get-classes-4-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 5), get_membership(res, k = 5))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1 p2  p3   p4   p5
#&gt; TCGA.N8.A4PN.01     4   0.000      1.000 0.00  0 0.0 1.00 0.00
#&gt; TCGA.N5.A59E.01     1   0.000      0.856 1.00  0 0.0 0.00 0.00
#&gt; TCGA.N7.A4Y5.01     2   0.000      1.000 0.00  1 0.0 0.00 0.00
#&gt; TCGA.N5.A4RD.01     1   0.000      0.856 1.00  0 0.0 0.00 0.00
#&gt; TCGA.N8.A4PQ.01     4   0.000      1.000 0.00  0 0.0 1.00 0.00
#&gt; TCGA.ND.A4WF.01     3   0.517      0.175 0.46  0 0.5 0.00 0.04
#&gt; TCGA.N5.A4R8.01     5   0.104      0.000 0.00  0 0.0 0.04 0.96
#&gt; TCGA.NA.A4QW.01     2   0.000      1.000 0.00  1 0.0 0.00 0.00
#&gt; TCGA.N9.A4Q1.01     2   0.000      1.000 0.00  1 0.0 0.00 0.00
#&gt; TCGA.N5.A4RV.01     3   0.000      0.649 0.00  0 1.0 0.00 0.00
#&gt; TCGA.N5.A4RO.01     1   0.000      0.856 1.00  0 0.0 0.00 0.00
#&gt; TCGA.NA.A4QX.01     4   0.000      1.000 0.00  0 0.0 1.00 0.00
#&gt; TCGA.N5.A4RA.01     4   0.000      1.000 0.00  0 0.0 1.00 0.00
#&gt; TCGA.NG.A4VU.01     2   0.000      1.000 0.00  1 0.0 0.00 0.00
#&gt; TCGA.N9.A4Q4.01     2   0.000      1.000 0.00  1 0.0 0.00 0.00
#&gt; TCGA.N5.A4RU.01     3   0.000      0.649 0.00  0 1.0 0.00 0.00
#&gt; TCGA.N7.A4Y8.01     1   0.418      0.302 0.60  0 0.0 0.40 0.00
#&gt; TCGA.N5.A4RM.01     1   0.000      0.856 1.00  0 0.0 0.00 0.00
</code></pre>

<script>
$('#tab-node-02-get-classes-4-a').parent().next().next().hide();
$('#tab-node-02-get-classes-4-a').click(function(){
  $('#tab-node-02-get-classes-4-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-02-get-classes-5'>
<p><a id='tab-node-02-get-classes-5-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 6), get_membership(res, k = 6))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1 p2   p3   p4 p5   p6
#&gt; TCGA.N8.A4PN.01     4  0.0000      1.000 0.00  0 0.00 1.00  0 0.00
#&gt; TCGA.N5.A59E.01     1  0.0000      0.913 1.00  0 0.00 0.00  0 0.00
#&gt; TCGA.N7.A4Y5.01     2  0.0000      1.000 0.00  1 0.00 0.00  0 0.00
#&gt; TCGA.N5.A4RD.01     1  0.1556      0.903 0.92  0 0.00 0.00  0 0.08
#&gt; TCGA.N8.A4PQ.01     4  0.0000      1.000 0.00  0 0.00 1.00  0 0.00
#&gt; TCGA.ND.A4WF.01     6  0.1480      0.000 0.04  0 0.02 0.00  0 0.94
#&gt; TCGA.N5.A4R8.01     5  0.0000      0.000 0.00  0 0.00 0.00  1 0.00
#&gt; TCGA.NA.A4QW.01     2  0.0000      1.000 0.00  1 0.00 0.00  0 0.00
#&gt; TCGA.N9.A4Q1.01     2  0.0000      1.000 0.00  1 0.00 0.00  0 0.00
#&gt; TCGA.N5.A4RV.01     3  0.0000      1.000 0.00  0 1.00 0.00  0 0.00
#&gt; TCGA.N5.A4RO.01     1  0.1556      0.903 0.92  0 0.00 0.00  0 0.08
#&gt; TCGA.NA.A4QX.01     4  0.0000      1.000 0.00  0 0.00 1.00  0 0.00
#&gt; TCGA.N5.A4RA.01     4  0.0000      1.000 0.00  0 0.00 1.00  0 0.00
#&gt; TCGA.NG.A4VU.01     2  0.0000      1.000 0.00  1 0.00 0.00  0 0.00
#&gt; TCGA.N9.A4Q4.01     2  0.0000      1.000 0.00  1 0.00 0.00  0 0.00
#&gt; TCGA.N5.A4RU.01     3  0.0000      1.000 0.00  0 1.00 0.00  0 0.00
#&gt; TCGA.N7.A4Y8.01     1  0.2512      0.842 0.88  0 0.00 0.06  0 0.06
#&gt; TCGA.N5.A4RM.01     1  0.0937      0.905 0.96  0 0.00 0.00  0 0.04
</code></pre>

<script>
$('#tab-node-02-get-classes-5-a').parent().next().next().hide();
$('#tab-node-02-get-classes-5-a').click(function(){
  $('#tab-node-02-get-classes-5-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-02-get-classes-6'>
<p><a id='tab-node-02-get-classes-6-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 7), get_membership(res, k = 7))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1 p2 p3 p4 p5   p6   p7
#&gt; TCGA.N8.A4PN.01     4   0.000      1.000 0.00  0  0  1  0 0.00 0.00
#&gt; TCGA.N5.A59E.01     1   0.387      0.456 0.60  0  0  0  0 0.02 0.38
#&gt; TCGA.N7.A4Y5.01     2   0.000      1.000 0.00  1  0  0  0 0.00 0.00
#&gt; TCGA.N5.A4RD.01     1   0.000      0.630 1.00  0  0  0  0 0.00 0.00
#&gt; TCGA.N8.A4PQ.01     4   0.000      1.000 0.00  0  0  1  0 0.00 0.00
#&gt; TCGA.ND.A4WF.01     6   0.000      0.000 0.00  0  0  0  0 1.00 0.00
#&gt; TCGA.N5.A4R8.01     5   0.000      0.000 0.00  0  0  0  1 0.00 0.00
#&gt; TCGA.NA.A4QW.01     2   0.000      1.000 0.00  1  0  0  0 0.00 0.00
#&gt; TCGA.N9.A4Q1.01     2   0.000      1.000 0.00  1  0  0  0 0.00 0.00
#&gt; TCGA.N5.A4RV.01     3   0.000      1.000 0.00  0  1  0  0 0.00 0.00
#&gt; TCGA.N5.A4RO.01     1   0.000      0.630 1.00  0  0  0  0 0.00 0.00
#&gt; TCGA.NA.A4QX.01     4   0.000      1.000 0.00  0  0  1  0 0.00 0.00
#&gt; TCGA.N5.A4RA.01     4   0.000      1.000 0.00  0  0  1  0 0.00 0.00
#&gt; TCGA.NG.A4VU.01     2   0.000      1.000 0.00  1  0  0  0 0.00 0.00
#&gt; TCGA.N9.A4Q4.01     2   0.000      1.000 0.00  1  0  0  0 0.00 0.00
#&gt; TCGA.N5.A4RU.01     3   0.000      1.000 0.00  0  1  0  0 0.00 0.00
#&gt; TCGA.N7.A4Y8.01     7   0.000      0.000 0.00  0  0  0  0 0.00 1.00
#&gt; TCGA.N5.A4RM.01     1   0.471      0.249 0.46  0  0  0  0 0.08 0.46
</code></pre>

<script>
$('#tab-node-02-get-classes-6-a').parent().next().next().hide();
$('#tab-node-02-get-classes-6-a').click(function(){
  $('#tab-node-02-get-classes-6-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-02-get-classes-7'>
<p><a id='tab-node-02-get-classes-7-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 8), get_membership(res, k = 8))
</code></pre>

<pre><code>#&gt;                 class entropy silhouette   p1 p2 p3 p4 p5 p6 p7   p8
#&gt; TCGA.N8.A4PN.01     4  0.0000      1.000 0.00  0  0  1  0  0  0 0.00
#&gt; TCGA.N5.A59E.01     1  0.1947      0.834 0.86  0  0  0  0  0  0 0.14
#&gt; TCGA.N7.A4Y5.01     2  0.0000      1.000 0.00  1  0  0  0  0  0 0.00
#&gt; TCGA.N5.A4RD.01     1  0.0471      0.912 0.98  0  0  0  0  0  0 0.02
#&gt; TCGA.N8.A4PQ.01     4  0.0000      1.000 0.00  0  0  1  0  0  0 0.00
#&gt; TCGA.ND.A4WF.01     6  0.0000      0.000 0.00  0  0  0  0  1  0 0.00
#&gt; TCGA.N5.A4R8.01     5  0.0000      0.000 0.00  0  0  0  1  0  0 0.00
#&gt; TCGA.NA.A4QW.01     2  0.0000      1.000 0.00  1  0  0  0  0  0 0.00
#&gt; TCGA.N9.A4Q1.01     2  0.0000      1.000 0.00  1  0  0  0  0  0 0.00
#&gt; TCGA.N5.A4RV.01     3  0.0000      1.000 0.00  0  1  0  0  0  0 0.00
#&gt; TCGA.N5.A4RO.01     1  0.0000      0.914 1.00  0  0  0  0  0  0 0.00
#&gt; TCGA.NA.A4QX.01     4  0.0000      1.000 0.00  0  0  1  0  0  0 0.00
#&gt; TCGA.N5.A4RA.01     4  0.0000      1.000 0.00  0  0  1  0  0  0 0.00
#&gt; TCGA.NG.A4VU.01     2  0.0000      1.000 0.00  1  0  0  0  0  0 0.00
#&gt; TCGA.N9.A4Q4.01     2  0.0000      1.000 0.00  1  0  0  0  0  0 0.00
#&gt; TCGA.N5.A4RU.01     3  0.0000      1.000 0.00  0  1  0  0  0  0 0.00
#&gt; TCGA.N7.A4Y8.01     7  0.0000      0.000 0.00  0  0  0  0  0  1 0.00
#&gt; TCGA.N5.A4RM.01     8  0.0000      0.000 0.00  0  0  0  0  0  0 1.00
</code></pre>

<script>
$('#tab-node-02-get-classes-7-a').parent().next().next().hide();
$('#tab-node-02-get-classes-7-a').click(function(){
  $('#tab-node-02-get-classes-7-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-node-02-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-node-02-consensus-heatmap'>
<ul>
<li><a href='#tab-node-02-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-node-02-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-node-02-consensus-heatmap-3'>k = 4</a></li>
<li><a href='#tab-node-02-consensus-heatmap-4'>k = 5</a></li>
<li><a href='#tab-node-02-consensus-heatmap-5'>k = 6</a></li>
<li><a href='#tab-node-02-consensus-heatmap-6'>k = 7</a></li>
<li><a href='#tab-node-02-consensus-heatmap-7'>k = 8</a></li>
</ul>
<div id='tab-node-02-consensus-heatmap-1'>
<pre><code class="r">consensus_heatmap(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-02-consensus-heatmap-1-1.png" alt="plot of chunk tab-node-02-consensus-heatmap-1"/></p>

</div>
<div id='tab-node-02-consensus-heatmap-2'>
<pre><code class="r">consensus_heatmap(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-02-consensus-heatmap-2-1.png" alt="plot of chunk tab-node-02-consensus-heatmap-2"/></p>

</div>
<div id='tab-node-02-consensus-heatmap-3'>
<pre><code class="r">consensus_heatmap(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-02-consensus-heatmap-3-1.png" alt="plot of chunk tab-node-02-consensus-heatmap-3"/></p>

</div>
<div id='tab-node-02-consensus-heatmap-4'>
<pre><code class="r">consensus_heatmap(res, k = 5)
</code></pre>

<p><img src="figure_cola/tab-node-02-consensus-heatmap-4-1.png" alt="plot of chunk tab-node-02-consensus-heatmap-4"/></p>

</div>
<div id='tab-node-02-consensus-heatmap-5'>
<pre><code class="r">consensus_heatmap(res, k = 6)
</code></pre>

<p><img src="figure_cola/tab-node-02-consensus-heatmap-5-1.png" alt="plot of chunk tab-node-02-consensus-heatmap-5"/></p>

</div>
<div id='tab-node-02-consensus-heatmap-6'>
<pre><code class="r">consensus_heatmap(res, k = 7)
</code></pre>

<p><img src="figure_cola/tab-node-02-consensus-heatmap-6-1.png" alt="plot of chunk tab-node-02-consensus-heatmap-6"/></p>

</div>
<div id='tab-node-02-consensus-heatmap-7'>
<pre><code class="r">consensus_heatmap(res, k = 8)
</code></pre>

<p><img src="figure_cola/tab-node-02-consensus-heatmap-7-1.png" alt="plot of chunk tab-node-02-consensus-heatmap-7"/></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-node-02-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-node-02-membership-heatmap'>
<ul>
<li><a href='#tab-node-02-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-node-02-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-node-02-membership-heatmap-3'>k = 4</a></li>
<li><a href='#tab-node-02-membership-heatmap-4'>k = 5</a></li>
<li><a href='#tab-node-02-membership-heatmap-5'>k = 6</a></li>
<li><a href='#tab-node-02-membership-heatmap-6'>k = 7</a></li>
<li><a href='#tab-node-02-membership-heatmap-7'>k = 8</a></li>
</ul>
<div id='tab-node-02-membership-heatmap-1'>
<pre><code class="r">membership_heatmap(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-02-membership-heatmap-1-1.png" alt="plot of chunk tab-node-02-membership-heatmap-1"/></p>

</div>
<div id='tab-node-02-membership-heatmap-2'>
<pre><code class="r">membership_heatmap(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-02-membership-heatmap-2-1.png" alt="plot of chunk tab-node-02-membership-heatmap-2"/></p>

</div>
<div id='tab-node-02-membership-heatmap-3'>
<pre><code class="r">membership_heatmap(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-02-membership-heatmap-3-1.png" alt="plot of chunk tab-node-02-membership-heatmap-3"/></p>

</div>
<div id='tab-node-02-membership-heatmap-4'>
<pre><code class="r">membership_heatmap(res, k = 5)
</code></pre>

<p><img src="figure_cola/tab-node-02-membership-heatmap-4-1.png" alt="plot of chunk tab-node-02-membership-heatmap-4"/></p>

</div>
<div id='tab-node-02-membership-heatmap-5'>
<pre><code class="r">membership_heatmap(res, k = 6)
</code></pre>

<p><img src="figure_cola/tab-node-02-membership-heatmap-5-1.png" alt="plot of chunk tab-node-02-membership-heatmap-5"/></p>

</div>
<div id='tab-node-02-membership-heatmap-6'>
<pre><code class="r">membership_heatmap(res, k = 7)
</code></pre>

<p><img src="figure_cola/tab-node-02-membership-heatmap-6-1.png" alt="plot of chunk tab-node-02-membership-heatmap-6"/></p>

</div>
<div id='tab-node-02-membership-heatmap-7'>
<pre><code class="r">membership_heatmap(res, k = 8)
</code></pre>

<p><img src="figure_cola/tab-node-02-membership-heatmap-7-1.png" alt="plot of chunk tab-node-02-membership-heatmap-7"/></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.






<script>
$( function() {
	$( '#tabs-node-02-get-signatures' ).tabs();
} );
</script>
<div id='tabs-node-02-get-signatures'>
<ul>
<li><a href='#tab-node-02-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-node-02-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-node-02-get-signatures-3'>k = 4</a></li>
<li><a href='#tab-node-02-get-signatures-4'>k = 5</a></li>
<li><a href='#tab-node-02-get-signatures-5'>k = 6</a></li>
<li><a href='#tab-node-02-get-signatures-6'>k = 7</a></li>
<li><a href='#tab-node-02-get-signatures-7'>k = 8</a></li>
</ul>
<div id='tab-node-02-get-signatures-1'>
<pre><code class="r">get_signatures(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-02-get-signatures-1-1.png" alt="plot of chunk tab-node-02-get-signatures-1"/></p>

</div>
<div id='tab-node-02-get-signatures-2'>
<pre><code class="r">get_signatures(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-02-get-signatures-2-1.png" alt="plot of chunk tab-node-02-get-signatures-2"/></p>

</div>
<div id='tab-node-02-get-signatures-3'>
<pre><code class="r">get_signatures(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-02-get-signatures-3-1.png" alt="plot of chunk tab-node-02-get-signatures-3"/></p>

</div>
<div id='tab-node-02-get-signatures-4'>
<pre><code class="r">get_signatures(res, k = 5)
</code></pre>

<p><img src="figure_cola/tab-node-02-get-signatures-4-1.png" alt="plot of chunk tab-node-02-get-signatures-4"/></p>

</div>
<div id='tab-node-02-get-signatures-5'>
<pre><code class="r">get_signatures(res, k = 6)
</code></pre>

<p><img src="figure_cola/tab-node-02-get-signatures-5-1.png" alt="plot of chunk tab-node-02-get-signatures-5"/></p>

</div>
<div id='tab-node-02-get-signatures-6'>
<pre><code class="r">get_signatures(res, k = 7)
</code></pre>

<p><img src="figure_cola/tab-node-02-get-signatures-6-1.png" alt="plot of chunk tab-node-02-get-signatures-6"/></p>

</div>
<div id='tab-node-02-get-signatures-7'>
<pre><code class="r">get_signatures(res, k = 8)
</code></pre>

<p><img src="figure_cola/tab-node-02-get-signatures-7-1.png" alt="plot of chunk tab-node-02-get-signatures-7"/></p>

</div>
</div>



Compare the overlap of signatures from different k:

```r
compare_signatures(res)
```

![plot of chunk node-02-signature_compare](figure_cola/node-02-signature_compare-1.png)

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

```r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

```r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```

If the signatures are defined as these which are uniquely high in current group, `diff_method` argument
can be set to `"uniquely_high_in_one_group"`:

```r
# code only for demonstration
tb = get_signature(res, k = ..., diff_method = "uniquely_high_in_one_group")
```




UMAP plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-node-02-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-node-02-dimension-reduction'>
<ul>
<li><a href='#tab-node-02-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-node-02-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-node-02-dimension-reduction-3'>k = 4</a></li>
<li><a href='#tab-node-02-dimension-reduction-4'>k = 5</a></li>
<li><a href='#tab-node-02-dimension-reduction-5'>k = 6</a></li>
<li><a href='#tab-node-02-dimension-reduction-6'>k = 7</a></li>
<li><a href='#tab-node-02-dimension-reduction-7'>k = 8</a></li>
</ul>
<div id='tab-node-02-dimension-reduction-1'>
<pre><code class="r">dimension_reduction(res, k = 2, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-02-dimension-reduction-1-1.png" alt="plot of chunk tab-node-02-dimension-reduction-1"/></p>

</div>
<div id='tab-node-02-dimension-reduction-2'>
<pre><code class="r">dimension_reduction(res, k = 3, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-02-dimension-reduction-2-1.png" alt="plot of chunk tab-node-02-dimension-reduction-2"/></p>

</div>
<div id='tab-node-02-dimension-reduction-3'>
<pre><code class="r">dimension_reduction(res, k = 4, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-02-dimension-reduction-3-1.png" alt="plot of chunk tab-node-02-dimension-reduction-3"/></p>

</div>
<div id='tab-node-02-dimension-reduction-4'>
<pre><code class="r">dimension_reduction(res, k = 5, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-02-dimension-reduction-4-1.png" alt="plot of chunk tab-node-02-dimension-reduction-4"/></p>

</div>
<div id='tab-node-02-dimension-reduction-5'>
<pre><code class="r">dimension_reduction(res, k = 6, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-02-dimension-reduction-5-1.png" alt="plot of chunk tab-node-02-dimension-reduction-5"/></p>

</div>
<div id='tab-node-02-dimension-reduction-6'>
<pre><code class="r">dimension_reduction(res, k = 7, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-02-dimension-reduction-6-1.png" alt="plot of chunk tab-node-02-dimension-reduction-6"/></p>

</div>
<div id='tab-node-02-dimension-reduction-7'>
<pre><code class="r">dimension_reduction(res, k = 8, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-02-dimension-reduction-7-1.png" alt="plot of chunk tab-node-02-dimension-reduction-7"/></p>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

```r
collect_classes(res)
```

![plot of chunk node-02-collect-classes](figure_cola/node-02-collect-classes-1.png)



If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html) for more detailed explanations.


 

## Session info


```r
sessionInfo()
```

```
#> R version 4.1.0 (2021-05-18)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: CentOS Linux 7 (Core)
#> 
#> Matrix products: default
#> BLAS/LAPACK: /usr/lib64/libopenblas-r0.3.3.so
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
#>  [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
#> [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] genefilter_1.74.0    ComplexHeatmap_2.8.0 markdown_1.1         knitr_1.33          
#> [5] matrixStats_0.59.0   cola_1.9.4          
#> 
#> loaded via a namespace (and not attached):
#>   [1] bitops_1.0-7           bit64_4.0.5            doParallel_1.0.16      RColorBrewer_1.1-2    
#>   [5] httr_1.4.2             GenomeInfoDb_1.28.1    data.tree_1.0.0        tools_4.1.0           
#>   [9] utf8_1.2.1             R6_2.5.0               irlba_2.3.3            DBI_1.1.1             
#>  [13] BiocGenerics_0.38.0    colorspace_2.0-2       GetoptLong_1.0.5       gridExtra_2.3         
#>  [17] tidyselect_1.1.1       bit_4.0.4              compiler_4.1.0         Biobase_2.52.0        
#>  [21] Cairo_1.5-12.2         xml2_1.3.2             microbenchmark_1.4-7   slam_0.1-48           
#>  [25] scales_1.1.1           askpass_1.1            stringr_1.4.0          digest_0.6.27         
#>  [29] XVector_0.32.0         pkgconfig_2.0.3        umap_0.2.7.0           fastmap_1.1.0         
#>  [33] highr_0.9              rlang_0.4.11           GlobalOptions_0.1.2    rstudioapi_0.13       
#>  [37] RSQLite_2.2.7          impute_1.66.0          generics_0.1.0         shape_1.4.6           
#>  [41] jsonlite_1.7.2         mclust_5.4.7           dplyr_1.0.7            dendextend_1.15.1     
#>  [45] RCurl_1.98-1.3         magrittr_2.0.1         GenomeInfoDbData_1.2.6 Matrix_1.3-4          
#>  [49] fansi_0.5.0            Rcpp_1.0.7             munsell_0.5.0          S4Vectors_0.30.0      
#>  [53] viridis_0.6.1          reticulate_1.20        lifecycle_1.0.0        scatterplot3d_0.3-41  
#>  [57] stringi_1.7.3          zlibbioc_1.38.0        blob_1.2.1             parallel_4.1.0        
#>  [61] crayon_1.4.1           lattice_0.20-44        Biostrings_2.60.1      splines_4.1.0         
#>  [65] annotate_1.70.0        circlize_0.4.13        KEGGREST_1.32.0        polylabelr_0.2.0      
#>  [69] pillar_1.6.1           rjson_0.2.20           codetools_0.2-18       stats4_4.1.0          
#>  [73] XML_3.99-0.6           glue_1.4.2             evaluate_0.14          png_0.1-7             
#>  [77] vctrs_0.3.8            foreach_1.5.1          polyclip_1.10-0        purrr_0.3.4           
#>  [81] gtable_0.3.0           openssl_1.4.4          assertthat_0.2.1       clue_0.3-59           
#>  [85] cachem_1.0.5           ggplot2_3.3.5          xfun_0.24              eulerr_6.1.0          
#>  [89] xtable_1.8-4           skmeans_0.2-13         RSpectra_0.16-0        viridisLite_0.4.0     
#>  [93] survival_3.2-11        tibble_3.1.2           Polychrome_1.3.1       iterators_1.0.13      
#>  [97] AnnotationDbi_1.54.1   memoise_2.0.0          IRanges_2.26.0         cluster_2.1.2         
#> [101] ellipsis_0.3.2         brew_1.0-6
```




