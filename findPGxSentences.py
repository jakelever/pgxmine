import kindred
import argparse
import pickle
import re

def filterCorpus(corpus,filterTerms):
	filtered = kindred.Corpus()
	for doc in corpus.documents:
		termsFound = any( ft in doc.text.lower() for ft in filterTerms )
		if termsFound:
			filtered.addDocument(doc)
	return filtered

def getNextSourceEntityID(doc):
	sourceEntityIDs = set( e.sourceEntityID for e in doc.entities)
	for i in range(1,1000000):
		candidate = "T%d" % i
		if not candidate in sourceEntityIDs:
			return candidate
	raise RuntimeError("Run out of candidate source entity IDs")


def annotateStarAlleles(corpus):
	for doc in corpus.documents:

		genes = [ e for e in doc.entities if e.entityType == 'Gene' ]

		for gene in genes:
			geneEnd = gene.position[0][1]
			geneID = gene.metadata['conceptid']

			offset = geneEnd
			regex = '^(,|and|or|/|\s|\+)*(?P<main>\*\s*[0-9]([\w:]*\w+)?)'
			star = re.search(regex,doc.text[offset:])
			while star:
				_,length = star.span()

				startPos,endPos = star.span('main')
				text = star.group('main')

				sourceEntityID = getNextSourceEntityID(doc)
				alleleName = text.strip()[1:].strip()

				conceptid = '*%s' % alleleName

				starAllele = kindred.Entity('Mutation',text,[(offset+startPos,offset+endPos)],sourceEntityID=sourceEntityID,metadata={'conceptid':conceptid,'associated_gene':geneID})
				doc.addEntity(starAllele)

				offset += length
				star = re.search(regex,doc.text[offset:])


if __name__ == '__main__':
	parser = argparse.ArgumentParser('Identify sentences that are potentially relevant to PharmGKB')
	parser.add_argument('--inBioc',required=True,type=str,help='BioC file annotated with GNormPlus, tmChem and tmVar')
	parser.add_argument('--filterTermsFile',required=True,type=str,help='Terms used to filter sentences to enrich for pharmacogenomics relations')
	parser.add_argument('--outBioc',required=True,type=str,help='Output BioC file with identified sentences')
	args = parser.parse_args()

	sentenceCorpus = kindred.Corpus()

	with open(args.filterTermsFile,'r') as f:
		filterTerms = [ line.strip().lower() for line in f ]

	for corpus in kindred.iterLoad('biocxml',args.inBioc):
		corpus = filterCorpus(corpus,filterTerms)

		for doc in corpus.documents:
			doc.text = doc.text.replace('\t',' ').replace('\n', ' ').replace('\r',' ')
			for e in doc.entities:
				e.text = e.text.replace('\t',' ').replace('\n', ' ').replace('\r',' ')

		annotateStarAlleles(corpus)

		corpusEntityTypes = [ set( e.entityType for e in doc.entities ) for doc in corpus.documents ]

		corpus.documents = [ doc for doc,types in zip(corpus.documents,corpusEntityTypes) if "Mutation" in types and "Chemical" in types ]

		parser = kindred.Parser(model='en_core_sci_sm')
		parser.parse(corpus)

		for doc in corpus.documents:
			for sentence in doc.sentences:
				sentenceLower = sentence.text.lower()
				hasAnyFilterTerm = any ( ft in sentenceLower for ft in filterTerms )

				entityTypes = set([ entity.entityType for entity,tokenIndices in sentence.entityAnnotations ])
				entityInfo = [ (e.entityType,e.text) for e,tokenIndices in sentence.entityAnnotations ]

				hasMutation = "Mutation" in entityTypes
				hasChemical = "Chemical" in entityTypes

				if hasMutation and hasChemical:
					sentenceStart = sentence.tokens[0].startPos
					
					sentenceEntities = [ kindred.Entity(e.entityType,e.text,[(e.position[0][0]-sentenceStart,e.position[0][1]-sentenceStart)],e.sourceEntityID,e.externalID,metadata=e.metadata) for e,_ in sentence.entityAnnotations ]
					sentenceEntities = [ e for e in sentenceEntities if e.position[0][0] >= 0 and e.position[0][1] < len(sentence.text) ]
					newDoc = kindred.Document(sentence.text,sentenceEntities,metadata=doc.metadata)
					sentenceCorpus.addDocument(newDoc)

	kindred.save(sentenceCorpus,'biocxml',args.outBioc)

	print("Found %d candidate sentences" % len(sentenceCorpus.documents))

