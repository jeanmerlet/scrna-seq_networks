import pandas as pd
import mygene



genes_path = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/kidney/healthy/kpmp/meta/kidney_ensembl_genes.txt'
genes = pd.read_csv(genes_path)

mg = mygene.MyGeneInfo()
test = mg.querymany(genes.values[:10], scopes='ensembl.gene', returnall=True)

print(test)
print(test['symbol'])




