import argparse
import json
import os
import kindred
import random
import pickle
from collections import defaultdict
from sklearn.linear_model import LogisticRegression
import numpy as np
import math
import time
import random
import itertools
import pgmine

def entitiesOverlap(candidateRelation):
	for e1,e2 in itertools.combinations(candidateRelation.entities,2):
		for start1,end1 in e1.position:
			for start2,end2 in e2.position:
				if not (end1 <= start2 or end2 <= start1 ):
					return True

	return False

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Merge sentence data and do some annotation')
	parser.add_argument('--mode',required=True,type=str,help='Whether to focus on star_allele, rs or other')
	parser.add_argument('--selectedChemicals',required=True,type=str,help='Which chemicals to filter for')
	parser.add_argument('--sentenceData',required=True,type=str,help='Directory with sentence data')
	parser.add_argument('--variantStopwords',required=True,type=str,help='Mutations to remove')
	parser.add_argument('--outPickle',required=True,type=str,help='Pickle to dump everything')
	parser.add_argument('--fileCount',required=False,type=int,help='Optionally limit the number of BioC file to load')
	parser.add_argument('--filterTerms',type=str,help='Optional file of filter terms')
	args = parser.parse_args()

	start = time.time()

	assert args.mode in ['star_allele','rs','other']

	with open(args.selectedChemicals) as f:
		drugData = json.load(f)
		drugMeshIDs = set([ 'MESH:'+d['MeSH'] for d in drugData ])

	filterTerms = set()
	if args.filterTerms:
		with open(args.filterTerms) as f:
			for line in f:
				filterTerms.add(line.strip().lower())

	with open(args.variantStopwords) as f:
		variantStopwords = set( line.strip().lower() for line in f )

	stopwords = {'dopamine','insulin','caffeine','nicotine','choline'}
	print("Loading sentences...")
	corpus = kindred.Corpus()
	filenames = sorted(os.listdir(args.sentenceData))

	if args.fileCount:
		filenames = random.sample(filenames, args.fileCount)

	for i,filename in enumerate(filenames):
		print(i, filename)

		tmpCorpus = kindred.load('biocxml',os.path.join(args.sentenceData,filename))
		if args.filterTerms:
			tmpCorpus.documents = [ doc for doc in tmpCorpus.documents if any ( filterTerm in doc.text.lower() for filterTerm in filterTerms ) ]
		tmpCorpus.documents = [ doc for doc in tmpCorpus.documents if not any ( stopword in doc.text.lower() for stopword in stopwords ) ]

		filtered = []
		for doc in tmpCorpus.documents:
			if args.mode == 'star_allele':
				doc.entities = [ e for e in doc.entities if not (e.entityType == 'Mutation' and not e.text.strip().startswith('*')) ]
			elif args.mode == 'rs':
				doc.entities = [ e for e in doc.entities if not (e.entityType == 'Mutation' and not e.text.strip().startswith('rs')) ]
			else:
				doc.entities = [ e for e in doc.entities if not (e.entityType == 'Mutation' and (e.text.strip().startswith('*') or e.text.strip().startswith('rs'))) ]

			doc.entities = [ e for e in doc.entities if e.position[0][0] >= 0 ]
			doc.entities = [ e for e in doc.entities if e.position[0][1] <= len(doc.text) ]
			doc.entities = [ e for e in doc.entities if not (e.entityType == 'Chemical' and not e.metadata['conceptid'] in drugMeshIDs) ]
			doc.entities = [ e for e in doc.entities if not (e.entityType == 'Chemical' and len(e.text) <= 4) ]
			doc.entities = [ e for e in doc.entities if not (e.entityType == 'Mutation' and pgmine.normalizeMutation(e.text) is None) ]
			doc.entities = [ e for e in doc.entities if not (e.entityType == 'Mutation' and e.text.lower() in variantStopwords) ]

			entityTypes = set( e.entityType for e in doc.entities )
			if 'Chemical' in entityTypes and 'Mutation' in entityTypes:
				filtered.append(doc)


		corpus.documents += filtered

	print("Parsing sentences...")
	parser = kindred.Parser(model='en_core_sci_sm')
	parser.parse(corpus)

	print("Generating candidates...")
	entityTypes = ('Chemical','Mutation')
	assert len(entityTypes) == 2, "Expected binary relations"
	candidateBuilder = kindred.CandidateBuilder(entityCount=2,acceptedEntityTypes=[entityTypes])
	candidateRelations = candidateBuilder.build(corpus)

	print("# candidateRelations = %d" % len(candidateRelations))
	candidateRelations = [ cr for cr in candidateRelations if not entitiesOverlap(cr) ]
	print("# candidateRelations = %d" % len(candidateRelations))

	print("Vectorizing...")
	vectorizer = kindred.Vectorizer()
	vectors = vectorizer.fit_transform(candidateRelations)
	vectors = vectors.tocsr()

	print("Writing pickle...")
	outData = corpus,candidateRelations,vectors#,allMetadata
	with open(args.outPickle,'wb') as outF:
		pickle.dump(outData,outF)

	end = time.time()
	duration = end-start

	print("TIMER\t%d\t%f" % (len(corpus.documents), duration))

