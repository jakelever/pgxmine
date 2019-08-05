import argparse
import csv
import os
import hashlib
from collections import defaultdict

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Filter Cancermine for more conservative predictions')
	parser.add_argument('--inData',required=True,type=str,help='Input directory with TSV files to be filtered')
	parser.add_argument('--outUnfiltered',required=True,type=str,help='Output unfiltered data')
	parser.add_argument('--outCollated',required=True,type=str,help='Output collated and filtered data')
	parser.add_argument('--outSentences',required=True,type=str,help='Output filtered sentences that match collated data')
	args = parser.parse_args()

	assert os.path.isdir(args.inData)

	threshold = 0.75

	collated = defaultdict(set)
	collatedMatchingID = {}

	#sentenceKeyFields = 'pmid,title,journal,journal_short,year,month,day,section,subsection,chemical_mesh_id,chemical_pharmgkb_id,chemical_drugbank_id,chemical_text,chemical_normalized,variant_id,variant_type,variant_text,variant_normalized,gene_ids,gene_names,score,sentence'
	collatedKeyFields = 'chemical_mesh_id,chemical_pharmgkb_id,chemical_drugbank_id,chemical_normalized,variant_id,variant_normalized,variant_type,gene_ids,gene_names'

	inputFiles = sorted( [ os.path.join(args.inData,f) for f in os.listdir(args.inData) if f.endswith('.tsv') ] )

	inputFilesHeader = None

	alreadySeen = set()

	recordCount,filteredRecordCount = 0,0
	with open(args.outUnfiltered,'w') as outUnfiltered, open(args.outSentences,'w') as outSentences:
		for inputFile in inputFiles:
			with open(inputFile) as inF:
				#inTSV = csv.reader(inF,delimiter='\t')
				
				#headers = next(inTSV, None)
				headers = inF.readline().strip('\n').split('\t')
				if inputFilesHeader is None:
					inputFilesHeader = headers
					outUnfiltered.write("\t".join(headers) + '\n')
					outSentences.write("matching_id\t" + "\t".join(headers) + '\n')
				else:
					assert inputFilesHeader == headers, "Headers don't match expected in file %s" % inputFile

				#for i,row in enumerate(inTSV):
				for i,line in enumerate(inF):
					row = line.strip('\n').split('\t')

					assert len(row) == len(headers), "Got %d columns, expected %d in row %d, file %s" % (len(row),len(headers),i+1,inputFile)
					r = { h:v for h,v in zip(headers,row) }

					score = float(r['score'])
					recordCount += 1

					pmid = r['pmid']

					key = (r['pmid'],r['formatted_sentence'])
					if key in alreadySeen:
						continue
					alreadySeen.add(key)

					outUnfiltered.write("\t".join(r[h] for h in headers) + "\n")

					keepIt = score > threshold and pmid != 'None'
					if keepIt:
						collatedKey = tuple( [ r[k] for k in collatedKeyFields.split(',') ] )
						collated[collatedKey].add(pmid)

						# Make a field using the key data that can be used to match between tables
						matchingID = hashlib.md5("|".join(list(collatedKey)).encode('utf-8')).hexdigest()
						collatedMatchingID[collatedKey] = matchingID

						outSentences.write(matchingID + "\t" + "\t".join(r[h] for h in headers) + "\n")
						filteredRecordCount += 1


	with open(args.outCollated,'w') as outF:
		headers = 'matching_id,%s,citation_count' % collatedKeyFields
		headerCount = len(headers.split(','))
		outF.write(headers.replace(',','\t') + '\n')

		collatedCounts = [ (len(pmids),key) for key,pmids in collated.items() ]
		collatedCounts = sorted(collatedCounts,reverse=True)
		for citation_count,collatedKey in collatedCounts:

			matchingID = collatedMatchingID[collatedKey]

			outData = [matchingID] + list(collatedKey) + [str(citation_count)]
			assert len(outData) == headerCount

			outLine = "\t".join(outData)
			outF.write(outLine + "\n")

	print("%d records filtered to %d sentences and collated to %d chemical/variant associations" % (recordCount, filteredRecordCount, len(collated)))
	print("Written to %s and %s" % (args.outSentences, args.outCollated))

