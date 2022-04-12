# Gene Ontology Fingerprint Construction

This repository mainly about codes and data of Gene Ontology Fingerprint  (GOF) construction of the study: **XGOF: an knowledge empowered gene ontology fingerprint approach to improve data interpretability and hypothesis generation**.

***

## 1. Introduction

For a case study, such as cancer, we construct gene ontology fingerprint with following steps:

### step 1. literature retrieval and corpus construction

For each case, we obtained its standard expression term and synonyms from MeSH database, then used these terms to retrieved PubMed database and get the relevant PubMed identifiers (PMIDs) and PubMed Central identifiers (PMCIDs). PubTator is a web-based system providing automatic annotations of biomedical concepts such as genes in PubMed abstracts and PMC full-text articles. We progarmatically retrieved annotations from PubTator with PMIDs or PMCIDs, and extratced gene annotations and corresponding sentences. The Gene Ontology (GO) knowledgebase curated the world’s largest source of information on the functions of genes. We downloaded the ontology file which is released in May 01 2021. Here we applied a lexicon look-up method with exact string matching strategy to recognize GO terms between the gene annotation sentences and ontology dictionary. We also considered word boundary detection and longest string matching rule.

#### sentence-gene mapping

```html
pmid_sentid,geneid
34682775_s4,6240
34682775_s5,1633;6240;978;2030
34682775_s1,1633;6240;978;2030
... ...
```

#### sentence-ontology mapping

```html
pmid_sentid,goid
34667161_s4,GO:0023052
34667161_s5,GO:0023052
34667161_s6,
34667161_s8,GO:0040007
... ...
```

### step 2. Generation of the GOF

**calculate P-value**

Given a gene and a GO term recognized in the corpus, we first calculated the statistical significance of the connectivity of each gene-GO pair using the unadjusted p-value generated from the hypergeometric test.

**adjust P-value**

Second, to convert the p-value into the adjusted p-value, we took into account several noised non-significant gene-GO relationships with the following three situations: i) abstracts that mentioned many genes; ii) abstracts that included general GO terms (e.g., “cell”, GO:0005623); and iii) hotspot genes or well-studied genes (e.g., TP53, and ERBB2) linked to massive biomedical studies. To eliminate these noised non-significant gene-GO relationships, we adopted the adjusted p-value defined by Tsoi LC.

### step 3. Calculate gene-gene similarity score (GGSS) 

The gene-gene similarity score (GGSS) was computed based on a modified inner product algorithm. The GGSS reflects the similarity between two genes, indicating their similar ontology fingerprints in a certain disease. Furthermore, we determined the GGSS percentile threshold (**top 1%**) .

To identify biologically justified subnetworks (highly dense interconnected clusters) in the GOF network,
we used the MCODE algorithm for network clustering and Cytoscape for visualization.

### step 4. Other application 

Finally, we applied GOF for multiple biological task, such as gene similarity network, tumor sample classification, gene-GO pattern discovery, and gene function prediction.

***

## 2. TCGA GOF

Descriptive statistics of GOF construction of 33 human related cancers from TCGA database. (PMID: the number of PubMed identifers, PMID-abstract: the number of PMIDs in which only annotations of abstracts can be accessed from PubTator, PMID-full text: the number of PMIDs in which we can download annotations of full text from PubTator. Gene sentence: the number of sentences annotated with gene concepts. Gene entities: the number of genes in GOF. GO terms: the number of GO terms in GOF.)

| **Cancers** | **PMID** | **PMID-abstract** | **PMID-full text** | **Gene setences** | **Gene entities** | **GO terms** |
| ----------- | -------- | ----------------- | ------------------ | ----------------- | ----------------- | ------------ |
| ACC         | 4113     | 3135              | 978                | 44110             | 706               | 294          |
| BLCA        | 69460    | 59954             | 9506               | 426480            | 3208              | 955          |
| BRCA        | 368476   | 284808            | 83668              | 4983147           | 8923              | 2565         |
| CESC        | 179930   | 143125            | 36805              | 1910774           | 6204              | 1736         |
| CHOL        | 47306    | 37873             | 9433               | 324880            | 2124              | 671          |
| COAD        | 245094   | 194048            | 51046              | 2591352           | 7244              | 1901         |
| DLBC        | 29098    | 22388             | 6710               | 282973            | 1884              | 551          |
| ESCA        | 73948    | 57429             | 16519              | 623866            | 3444              | 968          |
| GBM         | 45828    | 27874             | 17954              | 1244521           | 5262              | 1575         |
| HNSC        | 42017    | 31028             | 10989              | 650473            | 3648              | 1073         |
| KICH        | 3318     | 2838              | 480                | 20631             | 441               | 202          |
| KIRC        | 57426    | 43192             | 14234              | 701381            | 4011              | 1516         |
| KIRP        | 4810     | 3938              | 872                | 40152             | 558               | 248          |
| LAML        | 83768    | 66132             | 17636              | 1067931           | 4046              | 1179         |
| LGG         | 1931     | 1141              | 790                | 43285             | 809               | 325          |
| LIHC        | 320731   | 240511            | 80220              | 4380489           | 8461              | 2409         |
| LUAD        | 56632    | 40331             | 16301              | 986313            | 4655              | 1247         |
| LUSC        | 92005    | 69566             | 22439              | 1496063           | 5236              | 1380         |
| MESO        | 19388    | 15212             | 4176               | 169537            | 1617              | 520          |
| OV          | 131094   | 109642            | 21452              | 1318576           | 5311              | 1589         |
| PAAD        | 104126   | 80899             | 23227              | 1377115           | 5017              | 1334         |
| PCPG        | 31333    | 26585             | 4748               | 133009            | 748               | 331          |
| PRAD        | 143766   | 115072            | 28694              | 1827049           | 5415              | 1514         |
| READ        | 70386    | 61247             | 9139               | 188483            | 1333              | 497          |
| SARC        | 185006   | 151845            | 33161              | 1264890           | 5242              | 1555         |
| SKCM        | 55792    | 44857             | 10935              | 575319            | 2973              | 938          |
| STAD        | 120793   | 101385            | 19408              | 894986            | 4267              | 1147         |
| TGCT        | 15534    | 13696             | 1838               | 62495             | 963               | 330          |
| THCA        | 80917    | 66523             | 14394              | 635474            | 3407              | 1044         |
| THYM        | 11499    | 9317              | 2182               | 66106             | 756               | 330          |
| UCEC        | 35774    | 30359             | 5415               | 296419            | 2213              | 677          |
| UCS         | 7300     | 6299              | 1001               | 21607             | 301               | 151          |
| UVM         | 4478     | 2999              | 1479               | 72767             | 925               | 343          |

All gene and GO annotations and GOF results can be downloaded from zenodo.

## 3. Contacts

nadger_wang@139.com, Shanghai Eastern Hepatobiliary Surgery Hospital, Shanghai, 200438,
China

zonghui@tongji.edu.cn, Tongji University, Shanghai, 200092, China

