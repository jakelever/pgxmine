import re
import os
import argparse
import kindred

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Try to map star alleles to rs IDs with some basic text mining')
#	parser.add_argument('--inDir',required=True,type=str,help='Input directory of BioC xml files')
	parser.add_argument('--inFile',required=True,type=str,help='Input directory of BioC xml files')
	parser.add_argument('--genes',required=True,type=str,help='Gene file')
	parser.add_argument('--outFile',required=True,type=str,help='Output file')
	args = parser.parse_args()

	geneList = set()
	with open(args.genes) as f:
		for line in f:
			hugo_id,single,synonyms,entrez_id = line.strip('\n').split('\t')
			geneList.add(single.lower())

	#inputFiles = sorted( os.path.join(args.inDir,f) for f in os.listdir(args.inDir) if f.endswith('.xml') )

	with open(args.outFile,'w') as outF:
		for corpus in kindred.iterLoad('biocxml',args.inFile):

			corpus.documents = [ doc for doc in corpus.documents if 'rs' in doc.text and '*' in doc.text ]

			for doc in corpus.documents:
				for sentence in doc.text.split('.'):
					#matches = re.finditer('(?P<gene>\w+)\s*\*\s*(?P<star>\w+).{0,40}\((?P<rs>rs\d+)\)', doc.text)
					genes = []
					matches = re.finditer('(?P<gene>\w+)\s*\*\s*(?P<star>\w+)', sentence)
					for match in matches:
						startPos,endPos = match.span('gene')

						if match.group('gene').lower() in geneList:
							genes.append((startPos,match.group('gene')))

					genes = sorted(genes)

					matches = re.finditer('\*\s*(?P<star>\w+)(?P<inbetween>[^,]{0,40}?)\((?P<rs>rs\d+)\)', sentence)
					for match in matches:
						startPos = match.span()[0]
						gene = [ gene for p,gene in genes if p < startPos ]
						if len(gene) == 0:
							continue
						gene = gene[-1]
						#print("-"*30)
						#print(genes,gene,match.group(),match.groupdict())

						star = match.group('star')
						rsid = match.group('rs')
						inbetween = match.group('inbetween')

						out = [gene,star,rsid,inbetween]

						outF.write("\t".join(out) + '\n')
						


