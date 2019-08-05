import kindred
import re
import html

aminoAcidInfo = [('ALA','A'),('ARG','R'),('ASN','N'),('ASP','D'),('CYS','C'),('GLU','E'),('GLN','Q'),('GLY','G'),('HIS','H'),('ILE','I'),('LEU','L'),('LYS','K'),('MET','M'),('PHE','F'),('PRO','P'),('SER','S'),('THR','T'),('TRP','W'),('TYR','Y'),('VAL','V'),('ALANINE','A'), ('CYSTEINE','C'), ('ASPARTICACID','D'), ('GLUTAMICACID','E'), ('PHENYLALANINE','F'), ('GLYCINE','G'), ('HISTIDINE','H'), ('ISOLEUCINE','I'), ('LYSINE','K'), ('LEUCINE','L'), ('METHIONINE','M'), ('ASPARAGINE','N'), ('PROLINE','P'), ('GLUTAMINE','Q'), ('ARGININE','R'), ('SERINE','S'), ('THREONINE','T'), ('VALINE','V'), ('TRYPTOPHAN','W'), ('TYROSINE','Y'),('STOP','X'),('TER','X')]
aminoAcidMap = { big:small for big,small in aminoAcidInfo }
#for letter in string.ascii_lowercase:
#	aminoAcidMap[letter] = letter.upper()
for letter in 'ABCDEFGHIKLMNPQRSTVWYZX':
	aminoAcidMap[letter] = letter
aminoAcidMap['*'] = '*'
#aminoAcidMap['STOP'] = 'stop'

def normalizeMutation(mention):
	if mention.strip().startswith('*'):
		return mention.replace(' ','')
	elif mention.startswith('rs'):
		return mention.replace(' ','')
		

	examples = [
		('p.T790M',	'T790M'),
		('c.93G>A',	'G93A'),
		('c.93G>A',	'c.G93A'),
		('c.93G>A',	'c.93G>A'),
		('c.93G>A',	'c.93G/A'),
		('c.93G>A',	'93G>A'),
		('c.93G>A',	'G/A-93'),
		('c.93G>A',	'93G->A'),
		('c.93G>A',	'93G-->A'),
		('c.93G>A',	'G93->A'),
		('c.93G>A',	'G93-->A'),
		('c.93G>A',	'93G-A'),
		('c.93G>A',	'G modified A 93'),
		('c.93G>A',	'93G/A'),
		('c.93G>A',	'93,G/A'),
		('c.93G>A',	'(93) G/A'),
		('c.93G>A',	'93 (G/A)'),
		('c.93G>A',	'G to A substitution at nucleotide 93'),
		('c.93G>A',	'G to A substitution at position 93'),
		('c.93G>A',	'G to A at nucleotide 93'),
		('c.93G>A',	'G to A at position 93'),
		('c.93G>A',	'g+93G>A'),
		('c.93delG',	'c.93delG'),
		('c.93delG',	'c.93Gdel'),
		('c.93delG',	'93delG'),
		('c.93delG',	'93Gdel'),
		('c.93GGC>GAC',	'GGC93GAC'),
		('c.93_94del',	'c.93-94del'),
		('c.93_94del',	'c.93_94del'),
		('c.93_94del',	'93-94del'),
		('c.93_94del',	'93_94del'),
		('c.93dup',	'c.93dup'),
		('c.93_94dup',	'c.93-94dup'),
		('c.93_94dup',	'c.93_94dup'),
		('c.93_94dup',	'93-94dup'),
		('c.93_94dup',	'93_94dup'),
		('g.93G>A',	'g.93G>A'),
		('m.93G>A',	'm.93G>A'),
		('p.T790M',	'THR790MET'),
		('p.T790M',	'THR790/MET'),
		('p.T790M',	'THR790 to MET'),
		('p.T790M',	'THR-790 to MET'),
		('p.T790M',	'THR790-to-MET'),
		('p.T790M',	'THR790->MET'),
		('p.T790M',	'THR790-->MET'),
		('p.T790M',	'THR790-MET'),
		('p.T790M',	'THR790----MET'),
		('p.T790M',	'790THR----MET'),
		('p.T790M',	'THR-790-MET'),
		('p.T790M',	'THR-790MET'),
		('p.T790M',	'THR-790 -> MET'),
		('p.T790M',	'THR-790 --> MET'),
		('p.T790M',	'THR(790)MET'),
		('p.T790M',	'p.THR790MET'),
		('p.T790M',	'THR-to-MET substitution at position 790'),
		('p.T790M',	'THR 790 is replaced by MET'),
		('p.T790M',	'THR 790 mutated to MET'),
		('p.T790M',	'THR 790 was mutated to MET'),
		('p.T790M',	'THREONINE-to-METHIONINE mutation at residue 790'),
		('p.T790M',	'THREONINE-to-METHIONINE mutation at amino acid 790'),
		('p.T790M',	'THREONINE-to-METHIONINE mutation at amino acid position 790'),
		('p.T790M',	'THREONINE-to-METHIONINE mutation at position 790'),
		('p.T790M',	'THREONINE-to-METHIONINE mutation in position 790'),
		('p.T790M',	'THREONINE-to-METHIONINE substitution at residue 790'),
		('p.T790M',	'THREONINE-to-METHIONINE substitution at amino acid 790'),
		('p.T790M',	'THREONINE-to-METHIONINE substitution at amino acid position 790'),
		('p.T790M',	'THREONINE-to-METHIONINE substitution at position 790'),
		('p.T790M',	'THREONINE-to-METHIONINE substitution in position 790'),
		('p.T790M',	'THREONINE-to-METHIONINE alteration at residue 790'),
		('p.T790M',	'THREONINE-to-METHIONINE alteration at amino acid 790'),
		('p.T790M',	'THREONINE-to-METHIONINE alteration at amino acid position 790'),
		('p.T790M',	'THREONINE-to-METHIONINE alteration at position 790'),
		('p.T790M',	'THREONINE-to-METHIONINE alteration in position 790'),
		('p.T790M',	'THREONINE-to-METHIONINE change at residue 790'),
		('p.T790M',	'THREONINE-to-METHIONINE change at amino acid 790'),
		('p.T790M',	'THREONINE-to-METHIONINE change at amino acid position 790'),
		('p.T790M',	'THREONINE-to-METHIONINE change at position 790'),
		('p.T790M',	'THREONINE-to-METHIONINE change in position 790'),
		('p.T790M',	'THREONINE-to-METHIONINE at residue 790'),
		('p.T790M',	'THREONINE-to-METHIONINE at amino acid 790'),
		('p.T790M',	'THREONINE to METHIONINE mutation at residue 790'),
		('p.T790M',	'THREONINE to METHIONINE mutation at amino acid 790'),
		('p.T790M',	'THREONINE to METHIONINE mutation at amino acid position 790'),
		('p.T790M',	'THREONINE to METHIONINE mutation at position 790'),
		('p.T790M',	'THREONINE to METHIONINE mutation in position 790'),
		('p.T790M',	'THREONINE to METHIONINE substitution at residue 790'),
		('p.T790M',	'THREONINE to METHIONINE substitution at amino acid 790'),
		('p.T790M',	'THREONINE to METHIONINE substitution at amino acid position 790'),
		('p.T790M',	'THREONINE to METHIONINE substitution at position 790'),
		('p.T790M',	'THREONINE to METHIONINE substitution in position 790'),
		('p.T790M',	'THREONINE to METHIONINE alteration at residue 790'),
		('p.T790M',	'THREONINE to METHIONINE alteration at amino acid 790'),
		('p.T790M',	'THREONINE to METHIONINE alteration at amino acid position 790'),
		('p.T790M',	'THREONINE to METHIONINE alteration at position 790'),
		('p.T790M',	'THREONINE to METHIONINE alteration in position 790'),
		('p.T790M',	'THREONINE to METHIONINE change at residue 790'),
		('p.T790M',	'THREONINE to METHIONINE change at amino acid 790'),
		('p.T790M',	'THREONINE to METHIONINE change at amino acid position 790'),
		('p.T790M',	'THREONINE to METHIONINE change at position 790'),
		('p.T790M',	'THREONINE to METHIONINE change in position 790'),
		('p.T790M',	'THREONINE to METHIONINE at residue 790'),
		('p.T790M',	'THREONINE to METHIONINE at amino acid 790'),
		('p.T790M',	'THREONINE by METHIONINE at position 790'),
		('p.T790M',	'THREONINE-790-METHIONINE'),
		('p.T790M',	'THREONINE-790 -> METHIONINE'),
		('p.T790M',	'THREONINE-790 --> METHIONINE'),
		('p.T790M',	'THREONINE 790 METHIONINE'),
		('p.T790M',	'THREONINE 790 changed to METHIONINE'),
		('p.T790M',	'THREONINE-790 METHIONINE'),
		('p.T790M',	'THREONINE 790-METHIONINE'),
		('p.T790M',	'THREONINE 790 to METHIONINE'),
		('p.T790M',	'THREONINE 790 by METHIONINE'),
		('p.T790M',	'790 THREONINE to METHIONINE'),
		('p.T790M',	'METHIONINE for THREONINE at amino acid 790'),
		('p.T790M',	'METHIONINE for THREONINE at position 790'),
		('p.T790M',	'METHIONINE for THREONINE 790'),
		('p.T790M',	'METHIONINE-for-THREONINE at position 790'),
		('p.T790M',	'METHIONINE for THREONINE substitution at position 790'),
		('p.T790M',	'METHIONINE-for-THREONINE substitution at position 790'),
		('p.T790M',	'METHIONINE for a THREONINE at position 790'),
		('p.T790M',	'METHIONINE for an THREONINE at position 790'),

		('p.T790M',	'p.T790M'),
		('p.T790M',	'p.(T790M)'),
		('p.T790M',	'790T>M'),
		('p.T790M',	'790T->M'),
		('p.T790M',	'790T-->M'),
		('p.T790M',	'T790->M'),
		('p.T790M',	'T790-->M'),
		('p.T790fsX',	'T790fs'),
		('p.T790fsX791','p.T790fsX791'),
		('p.T790fsX791','p.THR790fsx791'),
		('p.T790fsX791','THR790fsx791'),
		('p.790delT',	'THR790del'),
		('p.790delT',	'p.T790del'),
		('p.790delT',	'p.790delT'),
		('p.790delT',	'T790del'),
		('p.790delT',	'790delT'),
	]

	mention = mention.replace(' ' ,'')

	for patternOut,patternIn in examples:
		regex = "^%s$" % re.escape(patternIn.replace(' ',''))

		mapping = [
			('THREONINE','(?P<from>Alanine|Cysteine|AsparticAcid|GlutamicAcid|Phenylalanine|Glycine|Histidine|Isoleucine|Lysine|Leucine|Methionine|Asparagine|Proline|Glutamine|Arginine|Serine|Threonine|Valine|Tryptophan|Tyrosine)'),
			('METHIONINE','(?P<to1>Alanine|Cysteine|AsparticAcid|GlutamicAcid|Phenylalanine|Glycine|Histidine|Isoleucine|Lysine|Leucine|Methionine|Asparagine|Proline|Glutamine|Arginine|Serine|Threonine|Valine|Tryptophan|Tyrosine)'),
			('THR','(?P<from>Ala|Arg|Asn|Asp|Cys|Glu|Gln|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)'),
			('MET','(?P<to1>Ala|Arg|Asn|Asp|Cys|Glu|Gln|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|X|\*|Ter|Stop)'),
			('790','(?P<num>[1-9][0-9]*)'),
			('791','(?P<num2>[1-9][0-9]*)'),
			('T','(?P<from>[ABCDEFGHIKLMNPQRSTVWYZ])'),
			('M','(?P<to1>([ABCDEFGHIKLMNPQRSTVWYZX\*]|stop))'),
			('E','(?P<to2>[ABCDEFGHIKLMNPQRSTVWYZX\*])'),
			('V','(?P<to3>[ABCDEFGHIKLMNPQRSTVWYZX\*])'),
			('GGC','(?P<from>[acgt]+)'),
			('GAC','(?P<to1>[acgt]+)'),
			('G','(?P<from>[acgt])'),
			('A','(?P<to1>[acgt])'),
			('C','(?P<to2>[acgt])'),
			('93','(?P<num>[\+\-]?[1-9][0-9\-\+]*)'),
			('94','(?P<num2>[\+\-]?[1-9][0-9\-]*)')
			]

		unique = {}
		for mapFrom,mapTo in mapping:
			unique[mapFrom] = "!!!%04d" % len(unique)
			regex = regex.replace(mapFrom, unique[mapFrom])
		for mapFrom,mapTo in mapping:
			regex = regex.replace(unique[mapFrom], mapTo)

		#if patternIn == '93G>A':
		#	print(regex)
		#	assert False
		#print(regex)
		match = re.match(regex, mention, re.IGNORECASE)
		if match:
			d = { key:value.upper() for key,value in match.groupdict().items() }
			if 'num' in d:
				d['num'] = d['num'].rstrip('-+')
			#	d['num'] = d['num'][0] + d['num'][1:].replace('-','_')

			if patternOut == 'c.G>A':
				return "c.%s>%s" % (d['from'],d['to1'])
			elif patternOut == 'c.93G>A':
				return "c.%s%s>%s" % (d['num'],d['from'],d['to1'])
			elif patternOut == 'c.93delG':
				return "c.%sdel%s" % (d['num'],d['from'])
			elif patternOut == 'c.GGC>GAC':
				return "c.%s>%s" % (d['from'],d['to1'])
			elif patternOut == 'c.93GGC>GAC':
				return "c.%s%s>%s" % (d['num'],d['from'],d['to1'])
			elif patternOut == 'c.93G>A,C':
				return "c.%s%s>%s,%s" % (d['num'],d['from'],d['to1'],d['to2'])
			elif patternOut == 'c.93_94del':
				return "c.%s_%sdel" % (d['num'],d['num2'])
			elif patternOut == 'c.93_94dup':
				return "c.%s_%sdup" % (d['num'],d['num2'])
			elif patternOut == 'c.93dup':
				return "c.%sdup" % d['num']
			elif patternOut == 'g.93G>A':
				return "g.%s%s>%s" % (d['num'],d['from'],d['to1'])
			elif patternOut == 'm.93G>A':
				return "m.%s%s>%s" % (d['num'],d['from'],d['to1'])
			elif patternOut == 'p.TM':
				return "p.%s%s" % (aminoAcidMap[d['from']],aminoAcidMap[d['to1']])
			elif patternOut == 'p.T790M':
				return "p.%s%s%s" % (aminoAcidMap[d['from']],d['num'],aminoAcidMap[d['to1']])
			elif patternOut == 'p.T790M/E':
				return "p.%s%s%s,%s" % (aminoAcidMap[d['from']],d['num'],aminoAcidMap[d['to1']],aminoAcidMap[d['to2']])
			elif patternOut == 'p.T790M/E/V':
				return "p.%s%s%s,%s,%s" % (aminoAcidMap[d['from']],d['num'],aminoAcidMap[d['to1']],aminoAcidMap[d['to2']],aminoAcidMap[d['to3']])
			elif patternOut == 'p.T790fsX':
				return "p.%s%sfsX" % (aminoAcidMap[d['from']],d['num'])
			elif patternOut == 'p.T790fsX791':
				return "p.%s%sfsX%s" % (aminoAcidMap[d['from']],d['num'],d['num2'])
			elif patternOut == 'p.790delT':
				return "p.%sdel%s" % (d['num'],aminoAcidMap[d['from']])


	return None

def getFormattedDoc(doc,entitiesToHighlight):
	charArray = [ html.escape(c) for c in doc.text ]

	for e in entitiesToHighlight:
		for startPos,endPos in e.position:
			try:
				charArray[startPos] = '<b>' + charArray[startPos]
				charArray[endPos-1] = charArray[endPos-1] + '</b>'
			except:
				print("ERROR in getFormattedSentence")
				print(doc.text)
				print(e.text)
				print(e.position)
				sys.exit(1)

	return "".join(charArray)
