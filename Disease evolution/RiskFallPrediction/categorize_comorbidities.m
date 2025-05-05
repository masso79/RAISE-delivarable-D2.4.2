function Tcom=categorize_comorbidities(diseasedescr)


Wmetabolic=readtable('./metabolico.xlsx');
Wcardiovascular=readtable('./cardio.xlsx');
Wpsichic=readtable('./psico.xlsx');
Wonconeuro=readtable('./onco.xlsx');
Wother=readtable('./other.xlsx');

Words={Wmetabolic Wcardiovascular Wpsichic Wonconeuro Wother};
C=cell(numel(Words),1);

for w=1:numel(Words)
	C{w}=contains(diseasedescr,Words{w}.malattie,'ignorecase',true);
end

Tcom=table(C{:},'variablenames',{'metabolic' 'cardiovascular' 'psychic' 'onco-neuro' 'other'});
