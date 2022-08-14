# PubMed+PubTator

Here we describe how to query PubMed.gov database, and how to download annotations from PubTator.

***

## 1. PubMed

For a simple example case named 'test', we curate the search term in PubMed.gov, and save as *./case/test/test.term.txt*. 
```txt
Adrenocortical carcinoma[Title/Abstract]
```
Run following code to obtain pmid and pmcid:

```python
cd lib
#python pubmed.py -c <case>
python pubmed.py -c test
``` 

The results are saved as *../case/test/test.pmid.txt* and *../case/test/test.pmcid.txt*. 
In this step, we used the pmcid-pmid convert file *PMC-ids.csv*, which is downloaded and saved in *./knol/PubMed/PMC-ids.csv*. By the way, the api is https://www.ncbi.nlm.nih.gov/pmc/tools/id-converter-api/ .
## 2. PubTator

Next, we download gene annotations from PubTator. The pmid is only annotated in abstract text, pmcid is annotated in full text.  

```python
cd lib
#python pubtator.py -c <case>
python pubtator.py -c test
```
The annotations are saved in folder *../case/test/PubTator/*. We further extract gene annotation and original sentences, which are saved as *../case/test/PubTator/test.pubtator.bioclist.anno.txt* and *../case/test/PubTator/test.pubtator.bioclist.sent.txt*. For each sentence, we set sentence id which save as *../case/test/test.sentid.txt*.

### 2.1. extract Gene annotations

```python
cd lib
#python pubtator2GENE.py -c <case>
python pubtator2GENE.py -c test
```

### 2.2 recognize Gene Ontology
```python
cd lib
#python pubtator2GO.py -c <case>
python pubtator2GO.py -c test
```


