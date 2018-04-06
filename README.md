##Intro
This is Flora's Code repo for ROP299 project.

[Script/Parser](./Script/Parser) contains parsers inplimented in slightly different algorithms. [Script/Analyser](./Script/Analyser) contains R scripts used for enumerating the pased sequences.
    
[Makefile](./Makefile) and [parser.py](./Script/Parser/parser.py) are currently under development. New version would automate fastq file splitting and self assembly of submission file. It would also allow easier extension and modification of the parser.

##Dependency
`simple_alignment.py` and `alignment_junction.py` is dependent on Biopython/1.70. If those two scipts are needed, please follow the [instruction](https://wiki.scinet.utoronto.ca/wiki/index.php/Python) to install biopython in home directory.



