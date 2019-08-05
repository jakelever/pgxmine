import argparse
import xml.etree.cElementTree as etree
import re
import SPARQLWrapper
import codecs
import six
from collections import defaultdict
import json
import csv
import sys

if __name__ == '__main__':
	parser = argparse.ArgumentParser('Process the DrugBank XML file for the data we need')
	parser.add_argument('--meshC',required=True,type=str,help='Mesh c*.bin ASCII file')
	parser.add_argument('--meshD',required=True,type=str,help='Mesh d*.bin ASCII file')
	parser.add_argument('--drugbank',required=True,type=str,help='DrugBank XML file')
	parser.add_argument('--pharmgkb',required=True,type=str,help='PharmGKB Drugs File')
	parser.add_argument('--outFile',required=True,type=str,help='Output JSON')
	args = parser.parse_args()

	mesh = {}
	with open(args.meshC) as f:
		record = {}
		synonyms = []
		for line in f:
			line = line.strip()
			if line == '*NEWRECORD' and record:
				mesh[record['name'].lower()] = record['id']
				for synonym in synonyms:
					mesh[synonym.lower()] = record['id']
				record = {}
				synonyms = []
			elif line.startswith('NM = '):
				record['name'] = line[len('NM = '):]
			elif line.startswith('UI = '):
				record['id'] = line[len('UI = '):]
			elif line.startswith('SY = '):
				synonymLine = line[len('SY = '):]
				synonyms.append( synonymLine.split('|')[0].strip() )

		if record:
			mesh[record['name'].lower()] = record['id']
			for synonym in synonyms:
				mesh[synonym.lower()] = record['id']

	with open(args.meshD) as f:
		record = {}
		for line in f:
			line = line.strip()
			if line == '*NEWRECORD' and record:
				mesh[record['name'].lower()] = record['id']
				for synonym in synonyms:
					mesh[synonym.lower()] = record['id']
				record = {}
				synonyms = []
			elif line.startswith('MH = '):
				record['name'] = line[len('MH = '):]
			elif line.startswith('UI = '):
				record['id'] = line[len('UI = '):]
			elif line.startswith('ENTRY = '):
				synonymLine = line[len('ENTRY = '):]
				synonyms.append( synonymLine.split('|')[0].strip() )
			elif line.startswith('PRINT ENTRY = '):
				synonymLine = line[len('PRINT ENTRY = '):]
				synonyms.append( synonymLine.split('|')[0].strip() )

		if record:
			mesh[record['name'].lower()] = record['id']
			for synonym in synonyms:
				mesh[synonym.lower()] = record['id']
			

	print("Loaded MeSH")

	pharmGKB_id2Name = {}
	pharmGKB_id2Mesh = {}
	pharmGKB_name2id = {}

	pharmGKB_hasMesh = {}
	#MeSH:C400278(erlotinib)
	csv.field_size_limit(sys.maxsize)
	with open(args.pharmgkb) as f:
		reader = csv.DictReader(f, delimiter='\t')
		for row in reader:
			pharmgkbID = row['PharmGKB Accession Id']
			name = row['Name']
			externalIDs = row['External Vocabulary']

			externalIDs = externalIDs.replace('"','').split(',')
			meshID = None
			for externalID in externalIDs:
				if externalID.startswith('MeSH:'):
					meshID = externalID[len('MeSH:'):].split('(')[0]
					break

			if not meshID is None:
				pharmGKB_id2Mesh[pharmgkbID] = meshID
			pharmGKB_name2id[name.lower()] = pharmgkbID

			pharmGKB_id2Name[pharmgkbID] = name
			pharmGKB_hasMesh[pharmgkbID] = False

		
	output = []

	pharmgkbInDrugbank = set()

	stopCategories = {'Elements','Adenine Nucleotides','Sweetening Agents','Salt Solutions','Supplements','Solvents','Electrolyte Solutions','Food Additives','Food','Lactates','Diluents','Gases','Mineral Supplements','Basic Lotions and Liniments','Phosphate salts','Potassium Salt'}
	allowed = {'fludarabine''nilotinib', 'diamorphine', 'cocaine', 'lopinavir', 'flucloxacillin', 'isoniazid', 'hydrochlorothiazide', 'tolbutamide', 'streptomycin', 'ampicillin', 'phenobarbital', 'azithromycin', 'tenofovir', 'peramivir', 'artemisinin', 'lapatinib', 'rociletinib'}
	not_allowed = {'Cholesterol','Carbon dioxide','Adenine','Guanine','Thymine','Uracil','Cytosine','Adenosine','Guanosine','5-Methyluridine','Uridine','Cytidine','Deoxyadenosine','Deoxyguanosine','Thymidine','Deoxyuridine','Deoxycytidine','Heparin','Hydrocortisone','Estradiol','Tretinoin','Testosterone','Progesterone','Melatonin'}

	nameRemap = {'Metamfetamine':'Methamphetamine'}

	namespace = '{http://www.drugbank.ca}'
	for event, elem in etree.iterparse(args.drugbank, events=('start', 'end', 'start-ns', 'end-ns')):
		if (event=='end' and elem.tag==namespace+'drug'):

			if len(elem) < 5: # This isn't a main drug record
				elem.clear()
				continue

			name = elem.find('./%sname' % namespace).text

			if name == nameRemap:
				name = nameRemap[name]

			nameLower = name.lower()
			drugbankID = elem.find("./%sdrugbank-id[@primary='true']" % namespace).text
			groups = [ group.text for group in elem.findall('./%sgroups/%sgroup' % (namespace,namespace)) ]

			externalIDs = {}
			for externalIDElem in elem.findall('./%sexternal-identifiers/%sexternal-identifier' % (namespace,namespace)):
				resource = externalIDElem.find('./%sresource' % namespace).text
				identifier = externalIDElem.find('./%sidentifier' % namespace).text
				externalIDs[resource] = identifier

			categories = [ c.text for c in elem.findall('./%scategories/%scategory/%scategory' % (namespace,namespace,namespace)) ]
			productNames = [ c.text for c in elem.findall('./%sproducts/%sproduct/%sname' % (namespace,namespace,namespace)) ]
			fdaLabelElems = elem.findall('./%sfda-label' % namespace)
			
			exclude = False

			if name in not_allowed:
				exclude = True
			elif any (c in stopCategories for c in categories):
				exclude = True
			elif len(productNames) == 0:
				exclude = True
			elif len(fdaLabelElems) == 0:
				exclude = True

			if exclude and not nameLower in allowed:
				print("EXCLUDED\t%s\t%s" % (drugbankID,name))
				elem.clear()
				continue

			if not 'PharmGKB' in externalIDs and nameLower in pharmGKB_name2id:
				externalIDs['PharmGKB'] = pharmGKB_name2id[nameLower]

			isCancerDrug = any('Antineoplastic' in c for c in categories)

			meshIDs = []
			nameLowerHydrochloride = "%s hydrochloride" % nameLower
			if nameLower in mesh:
				meshIDs.append( mesh[nameLower] )
			elif nameLowerHydrochloride in mesh:
				meshIDs.append( mesh[nameLowerHydrochloride] )

			if 'PharmGKB' in externalIDs:
				pharmGKBID = externalIDs['PharmGKB']
				if pharmGKBID in pharmGKB_id2Mesh:
					meshIDs.append( pharmGKB_id2Mesh[pharmGKBID] )

			meshIDs = sorted(list(set(meshIDs)))

			if 'PharmGKB' in externalIDs and len(meshIDs) > 0:
				pharmGKB_hasMesh[externalIDs['PharmGKB']] = True

			tmp = dict(externalIDs)
			tmp['name'] = name
			tmp['DrugBank'] = drugbankID
			tmp['MeSH'] = meshIDs
			tmp['isCancerDrug'] = isCancerDrug

			output.append(tmp)

			elem.clear()

	print("Loaded DrugBank")


	for pharmGKBID,hasMesh in pharmGKB_hasMesh.items():
		if not hasMesh:
			name = pharmGKB_id2Name[pharmGKBID]
			print("No MESH: %s (%s)" % (name, pharmGKBID))


	with open(args.outFile,'w') as outF:
		json.dump(output,outF,indent=2)

