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

def now():
	return time.strftime("%Y-%m-%d %H:%M:%S")

class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'

class RESPONSE:
	POSITIVE = 1
	NEGATIVE = 0
	ENTITYERROR = -1
		
	TABLE = {'y':POSITIVE,'n':NEGATIVE,'x':ENTITYERROR}

def promptRelation(cr, num):
	assert isinstance(cr, kindred.CandidateRelation)
	entityIDToEntity = { entity.entityID:entity for entity,tokenIndices in cr.sentence.entityAnnotations }

	keyword,geneOrProtein = cr.entities

	keywordStart,keywordEnd = keyword.position[0]
	geneOrProteinStart,geneOrProteinEnd = geneOrProtein.position[0]

	charByChar = list(cr.sentence.text)
	charByChar[keywordStart] = bcolors.FAIL + charByChar[keywordStart]
	charByChar[keywordEnd-1] += bcolors.ENDC
	charByChar[geneOrProteinStart] = bcolors.HEADER + charByChar[geneOrProteinStart]
	charByChar[geneOrProteinEnd-1] += bcolors.ENDC
	
	sentence = "".join(charByChar)

	print()
	print(str(num) + ' ' + '#'*30)
	print(sentence)

	response = None
	while not response in RESPONSE.TABLE:
		response = input('Good? Positive=y,Negative=n,EntityError=x: ')

	return RESPONSE.TABLE[response]

def classifyForProbs(responses, vectors):
	trainIndices = sorted(list(responses.keys()))

	trainIndices = [ i for i in trainIndices if not responses[i] == RESPONSE.ENTITYERROR ]

	trainX = vectors[trainIndices,:]
	trainY = [ 1 if responses[i] == RESPONSE.POSITIVE else 0 for i in trainIndices ]

	clf = LogisticRegression(class_weight='balanced',random_state=1)
	clf.fit(trainX,trainY)
	assert list(clf.classes_) == [0,1]
	
	probs = clf.predict_proba(vectors)
	return probs[:,1].tolist()

def bootstrapProbs(allResponses, vectors):
	bootstrapCount = 10
	bootstrapFraction = 0.8
	countRequirement = 5

	probs = defaultdict(list)
	for _ in range(bootstrapCount):
		#print("bootstrapping...")
		subsetSize = int(round(bootstrapFraction*len(allResponses)))
		subsetKeys = random.sample(list(allResponses.keys()), subsetSize)
		subset = { k:allResponses[k] for k in subsetKeys }
		tmpProbs = classifyForProbs(subset, vectors)
		for index,prob in enumerate(tmpProbs):
			probs[index].append(prob)

	scores = [ (np.mean( [ abs(v-0.5) for v in vals ] ),index) for index,vals in probs.items() if len(vals) >= countRequirement ]

	ambig = [ 1 if (max(vals) > 0.5 and min(vals) < 0.5) else 0 for index,vals in probs.items() ]
	ambigCount = sum(ambig)
	ambigPerc = round(100*ambigCount/len(ambig),1)
	print()
	print("Ambig: %f%% (%d/%d)" % (ambigPerc,ambigCount,len(ambig)))

	scores = sorted(scores)
	return scores

def writeCorpus(allResponses, candidateRelations, corpus, annotationOrder, outCorpusFilename, threshold=0.5):

	entitiesToOrder = {}
	for i,sentence in enumerate(annotationOrder):
		for e,_ in sentence.entityAnnotations:
			entitiesToOrder[e] = i+1
			
	sentenceToIndex = { s:(i+1) for i,s in enumerate(annotationOrder) }

	known = [ i for i,response in allResponses.items() if response == RESPONSE.POSITIVE ]

	relations = []
	for i in known:
		cr = candidateRelations[i]
		relation = kindred.Relation('related',cr.entities)
		relations.append(relation)

	entitiesToDoc = {}
	for i,doc in enumerate(corpus.documents):
		for e in doc.entities:
			entitiesToDoc[e] = i
	
	annotatedDocs = set()
	entityErrorDocs = set()
	for i,response in allResponses.items():
		cr = candidateRelations[i]
		docID = entitiesToDoc[cr.entities[0]]
		doc = corpus.documents[docID]
		doc.metadata['annotationOrder'] = entitiesToOrder[cr.entities[0]]
		annotatedDocs.add(doc)

		if response==RESPONSE.ENTITYERROR:
			entityErrorDocs.add(corpus.documents[docID])

	annotatedDocs -= entityErrorDocs

	corpus.removeRelations()
	for relation in relations:
		docIDs = [ entitiesToDoc[e] for e in relation.entities ]
		docIDs = list(set(docIDs))
		assert len(docIDs) > 0, "Predicted relation contains entities that don't match any documents in corpus"
		assert len(docIDs) == 1, "Predicted relation contains entities that are spread across documents"

		docID = docIDs[0]
		if not relation in corpus.documents[docID].relations:
			corpus.documents[docID].addRelation(relation)

	outCorpus = kindred.Corpus()
	outCorpus.documents = list(annotatedDocs)
	kindred.save(outCorpus,'biocxml',outCorpusFilename)

	return outCorpus


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Merge sentence data and do some annotation')
	parser.add_argument('--inPickle',required=True,type=str,help='Pickle created by prepareToLearn')
	parser.add_argument('--mode',required=False,default='active',type=str,help='Mode for annotation')
	parser.add_argument('--inTest',required=True,type=str,help='Test set')
	parser.add_argument('--outAnnotated',required=True,type=str,help='BioC file with annotated sentences')
	args = parser.parse_args()

	assert args.mode in ['active','random']

	print("Loading data pickle...")
	with open(args.inPickle,'rb') as inF:
		corpus,candidateRelations,vectors = pickle.load(inF)
		#candidateRelations,vectors,allMetadata = pickle.load(inF)

	testCorpus = kindred.load('biocxml',args.inTest)


	print("Annotating and doing magic...")

	groupedBySentence = defaultdict(list)
	for i,candidateRelation in enumerate(candidateRelations):
		groupedBySentence[candidateRelation.sentence].append((i,candidateRelation))

	# Sort relations by their entity start locations
	groupedBySentence = { s:sorted(relations, key=lambda x : tuple( [ e.position[0] for e in x[1].entities ] ) ) for s,relations in groupedBySentence.items() }

	testTexts = set( doc.text for doc in testCorpus.documents )
	groupedBySentence = { s:relations for s,relations in groupedBySentence.items() if not s.text in testTexts }

	sentences = list(groupedBySentence.keys())
	annotationOrder = []
	
	allResponses = {}
	while True:
		sentence = random.choice(sentences)
		print("Chosen a sentence!")
		for i,cr in groupedBySentence[sentence]:
			#i,cr = random.choice(
			#i = random.randint(0,len(candidateRelations)-1)
			cr = candidateRelations[i]
			allResponses[i] = promptRelation(cr,len(annotationOrder)+1)
		sentences.remove(sentence)
		annotationOrder.append(sentence)

		positiveCount = sum ( 1 for _,response in allResponses.items() if response == RESPONSE.POSITIVE )
		negativeCount = sum ( 1 for _,response in allResponses.items() if response == RESPONSE.NEGATIVE )

		if positiveCount >=5 and negativeCount >= 5 and len(allResponses) > 10:
			if args.mode == 'active':
				uncertainty = bootstrapProbs(allResponses,vectors)
				#for std_dev,i in uncertainty[:1]:
				#	cr = candidateRelations[i]
				#	allResponses[i] = promptRelation(cr)
				std_dev,i = uncertainty[0]
				sentence = candidateRelations[i].sentence
			else:
				sentence = random.choice(sentences)

			for i,cr in groupedBySentence[sentence]:
				cr = candidateRelations[i]
				allResponses[i] = promptRelation(cr,len(annotationOrder)+1)
			sentences.remove(sentence)
			annotationOrder.append(sentence)


			#buildKnowledgebase(allResponses, vectors, candidateRelations, allMetadata, ID2Term, outFile, outPickle)
			writeCorpus(allResponses, candidateRelations, corpus, annotationOrder, args.outAnnotated)

			outCorpus = kindred.load('biocxml',args.outAnnotated)
			predictionCorpus = testCorpus.clone()
			predictionCorpus.removeRelations()

			parser = kindred.Parser(model='en_core_sci_sm')
			parser.parse(outCorpus)
			parser.parse(predictionCorpus)

			classifier = kindred.RelationClassifier()
			classifier.train(outCorpus)

			classifier.predict(predictionCorpus)

			kindred.evaluate(testCorpus,predictionCorpus,display=True)

