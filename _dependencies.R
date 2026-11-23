# Packages CNAscreen needs that renv cannot detect automatically.
# QDNAseq loads the genome annotation package internally via
# getBinAnnotations(), so it is listed here to keep it in renv.lock.
# This file is never run.
library(QDNAseq.hg38)