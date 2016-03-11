#!/usr/bin/python
# -*- coding: utf-8 -*-

import re, itertools

def test_regex_mirna( regex_mirna ):
	unmatched = []
	rxmirna   = re.compile( regex_mirna, re.IGNORECASE )
	with open( 'data/mirnames.tsv', 'r' ) as f:
		for line in f:
			mirname = line.strip()
			tt      = rxmirna.match( mirname )
			if tt.start() != 0 or tt.end() != len( mirname ):
				unmatched.append( mirname )
	return unmatched;

def mirna_detector( strin ):
	mirna_pieces = '(hsa[- ])?((?:let|(?:mi(?:r(?:na)?|cro[ -]?rna))))?([- ]\d+[a-z]*)?([- ](?:(?:[0-9.]+[pl]*)|(?:x(?:.[12])?)|[a]?s))?((?:[- ][35]p)|\*)?'
	first_mirna  = '(hsa[- ])?((?:let|(?:mi(?:r(?:na)?|cro[ -]?rna))))([- ]\d+\w*)([- ](?:(?:[0-9.]+[pl]*)|(?:x(?:.[12])?)|[a]?s))?((?:[- ][35]p)|\*)?'
	second_mirna = '(hsa[- ])?((?:let|(?:mi(?:r(?:na)?|cro[ -]?rna))))?([- ]\d+\w*)?([- ](?:(?:[0-9.]+[pl]*)|(?:x(?:.[12])?)|[a]?s))?((?:[- ][35]p)|\*)?'
	andcomma     = '(?:[ ]*[,]?[ ]+(?:and[ ]+)?)'
	enumeration  = '((?:'+first_mirna+')(?:(?:'+andcomma+'(?:'+second_mirna+'))*))'

	#assert( len( test_regex_mirna( mirna_pieces ) ) == 0 )
	#assert( len( test_regex_mirna( first_mirna  ) ) == 0 )
	#assert( len( test_regex_mirna( second_mirna ) ) == 0 )

	rx1     = re.compile(enumeration, re.IGNORECASE)
	rx2     = re.compile(',|and', re.IGNORECASE)
	rx3     = re.compile(mirna_pieces, re.IGNORECASE)	

	def regularize_mir( pieces ):
		pieces[0] = 'hsa-'
		if pieces[1].lower() != 'let':
			pieces[1] = 'mir'
		else:
			pieces[1] = 'let'
		for i in [ 2,3,4 ]:
			if pieces[i]:
				pieces[i] = re.sub(' ', '-', pieces[i].lower() )
		return ''.join( [ s for s in pieces if s!=None ] )

	def complete_enum( enum ):
		comp_enum = []
		for mirna in enum:
			pieces = list( rx3.match( mirna ).groups() )
			if not any( pieces ):
				pieces = rx3.match( '-'+mirna ).groups()
			if pieces[1] == None:
				last_filled  = [ i for i,p in enumerate( last_pieces ) if p != None ]
				old_values   = [ p for p in last_pieces if p != None ]
				new_values   = [ p for p in pieces if p != None ]
				values       = old_values[:len(old_values)-len(new_values)] + new_values
				pieces       = [ None ] * len( last_pieces )
				for i in range(len(last_filled)):
					pieces[ last_filled[i] ] = values[i]
			last_pieces = pieces	
			comp_enum.append( regularize_mir( pieces ) )
		return comp_enum


	enums = [ [ s.strip() for s in rx2.split( enum[0] ) if len( s.strip() )>0 ] for enum in rx1.findall(strin) ]
	crnums= [ complete_enum( e ) for e in enums ]
	mirs  = list(itertools.chain.from_iterable(crnums))
	return sorted(list(set(mirs)))

def main():
    for sentence in read_input(sys.stdin):
        print ','.join( mirna_detector(sentence) )

if __name__ == '__main__':
    main()