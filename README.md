# LAMPLINK

The LAMPLINK can detect statistically significant epistatic interactions of two or more SNPs from GWAS data. 
This software can be used in the same way as the widely used GWAS analysis software [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/), but LAMPLINK has the additional options for the detection of epistatic interactions with [LAMP](http://a-terada.github.io/lamp/), which is a multiple testing procedure for combinatorial effects discovery.
You can apply LAMPLINK to an analysis pipeline with PLINK simply by replacing ``plink`` with ``lamplink`` and adding the ``--lamp`` option. 

Please see [the manual page](http://a-terada.github.io/lamplink/) for more details. 
