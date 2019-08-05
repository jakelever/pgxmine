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
		#print(doc.text)

		genes = [ e for e in doc.entities if e.entityType == 'Gene' ]

		for gene in genes:
			geneEnd = gene.position[0][1]
			geneID = gene.metadata['conceptid']

			offset = geneEnd
			#print(doc.text[offset:])
			regex = '^(,|and|or|/|\s|\+)*(?P<main>\*\s*[0-9]([\w:]*\w+)?)'
			star = re.search(regex,doc.text[offset:])
			while star:
				#print(star)
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


def annotateStarAlleles2(corpus):
	for doc in corpus.documents:
		#print(doc.text)
		#print('-'*30)
		geneEnds = set([ e.position[0][1] for e in doc.entities if e.entityType == 'Gene' ])

		if geneEnds:
			for candidate in re.finditer('\s*\*\s*\w+', doc.text):
				#print(candidate)
				#print(dir(candidate))
				startPos,endPos = candidate.span()
				text = candidate.group()
				#print(startPos, geneEnds)

				if startPos in geneEnds:
					sourceEntityID = getNextSourceEntityID(doc)
					alleleName = text.strip()[1:].strip()
					conceptid = '* %s' % alleleName
					starAllele = kindred.Entity('Mutation',text,[(startPos,endPos)],sourceEntityID=sourceEntityID,metadata={'conceptid':conceptid})
					doc.addEntity(starAllele)
					#doc.entities.append(starAllele)

					#for sentence in doc.sentences:
					#	overlappingTokens = [ i for i,t in enumerate(sentence.tokens) if not (t.endPos < startPos or t.startPos > endPos) ]
					#	if overlappingTokens:
					#		sentence.addEntityAnnotation(starAllele,overlappingTokens)

			# Double check that there aren't any duplicate sourceEntityIDs
			assert len(doc.entities) == len(set( e.sourceEntityID for e in doc.entities) )

if __name__ == '__main__':
	parser = argparse.ArgumentParser('Identify sentences that are potentially relevant to PharmGKB')
	parser.add_argument('--inBioc',required=True,type=str,help='BioC file annotated with GNormPlus, tmChem and tmVar')
	parser.add_argument('--filterTermsFile',required=True,type=str,help='Terms used to filter sentences to enrich for pharmacogenomics relations')
	#parser.add_argument('--selectedChemicalsFile',required=True,type=str,help='List of selected chemicals to use')
	parser.add_argument('--outBioc',required=True,type=str,help='Output BioC file with identified sentences')
	args = parser.parse_args()

	sentenceCorpus = kindred.Corpus()

	with open(args.filterTermsFile,'r') as f:
		filterTerms = [ line.strip().lower() for line in f ]

	#with open(args.selectedChemicalsFile,'r') as f:
	#	selectedChemicals = set([ line.strip().lower() for line in f ])

	for corpus in kindred.iterLoad('biocxml',args.inBioc):
		corpus = filterCorpus(corpus,filterTerms)

		#for doc in corpus.documents:
		#	doc.entities = [ e for e in doc.entities if (e.entityType != 'Chemical' or e.text.lower() in selectedChemicals) ]
		#corpus.documents = [ doc for doc in corpus.documents if len(doc.entities) > 0 ]
		annotateStarAlleles(corpus)

		corpusEntityTypes = [ set( e.entityType for e in doc.entities ) for doc in corpus.documents ]

		corpus.documents = [ doc for doc,types in zip(corpus.documents,corpusEntityTypes) if "Mutation" in types and "Chemical" in types ]

		parser = kindred.Parser(model='en_core_sci_sm')
		parser.parse(corpus)

		for doc in corpus.documents:
			for sentence in doc.sentences:
				sentenceLower = sentence.text.lower()
				hasAnyFilterTerm = any ( ft in sentenceLower for ft in filterTerms )

				#if not hasAnyFilterTerm:
				#	continue

				#if sentence.text.count(';') >= 5:
				#	continue

				entityTypes = set([ entity.entityType for entity,tokenIndices in sentence.entityAnnotations ])
				entityInfo = [ (e.entityType,e.text) for e,tokenIndices in sentence.entityAnnotations ]
				#print(hasAnyFilterTerm, entityInfo)
				#print(sentence.text)
				#print('-'*30)

				hasMutation = "Mutation" in entityTypes
				hasChemical = "Chemical" in entityTypes
				#print(sentence.text, hasVariant, hasChemical)

				if hasMutation and hasChemical:
					sentenceStart = sentence.tokens[0].startPos
					
					sentenceEntities = [ kindred.Entity(e.entityType,e.text,[(e.position[0][0]-sentenceStart,e.position[0][1]-sentenceStart)],e.sourceEntityID,e.externalID,metadata=e.metadata) for e,_ in sentence.entityAnnotations ]
					newDoc = kindred.Document(sentence.text,sentenceEntities,metadata=doc.metadata)
					sentenceCorpus.addDocument(newDoc)

	kindred.save(sentenceCorpus,'biocxml',args.outBioc)

