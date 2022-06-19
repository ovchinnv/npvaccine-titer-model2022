% nucleotide and amino acid code conversion functions

aas={ 'ALA' 'ARG' 'ASP' 'GLN' 'LEU' 'THR' 'GLU' 'ILE' 'PHE' 'LYS' 'SER' 'VAL' 'MET' 'ASN' 'PRO' 'TYR' 'HIS' 'HSD' 'HSE' 'GLY' 'TRP' 'CYS' ;
      'A'   'R'   'D'   'Q'   'L'   'T'   'E'   'I'   'F'   'K'   'S'   'V'   'M'   'N'   'P'   'Y'   'H'   'H'   'H'   'G'   'W'   'C' } ;
aa1= @(x) char(aas(find(ismember(aas(:), x),1)+1)); % get 1 letter code from 3 letter code
aa3= @(x) char(aas(find(ismember(aas(:), x),1)-1)); % get 3 letter code from 1 letter code

nts={ 'ADE' 'THY' 'CYT' 'GUA' 'URA' ;
      'A'   'T'   'C'   'G'   'U' } ;
nt1= @(x) char(nts(find(ismember(nts(:), x))+1)); % get 1 letter code from 3 letter code
nt3= @(x) char(nts(find(ismember(nts(:), x))-1)); % get 3 letter code from 1 letter code
