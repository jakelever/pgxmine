import kindred
import argparse
import numpy as np

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Create PR curves using a train and test set')
	parser.add_argument('--train',required=True,type=str,help='Directory containing stand-off training test')
	parser.add_argument('--test',required=True,type=str,help='Directory containing stand-off testing test')
	args = parser.parse_args()

	print("threshold\tprecision\trecall")
	for threshold in np.arange(0,1.01,0.01):

		trainCorpus = kindred.load('standoff',args.train)
		testCorpus = kindred.load('standoff',args.test)

		predCorpus = testCorpus.clone()
		predCorpus.removeRelations()

		parser = kindred.Parser(model='en_core_sci_sm')
		parser.parse(trainCorpus)
		parser.parse(testCorpus)
		parser.parse(predCorpus)

		classifier = kindred.RelationClassifier(classifierType='LogisticRegression',threshold=threshold,acceptedEntityTypes=[('Chemical','Mutation')])
		classifier.train(trainCorpus)
		classifier.predict(predCorpus)

		precision,recall,f1score = kindred.evaluate(testCorpus,predCorpus,metric='all',display=False)

		print("%f\t%f\t%f" % (threshold,precision,recall))

