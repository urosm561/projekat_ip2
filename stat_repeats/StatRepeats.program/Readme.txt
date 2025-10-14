Please cite:
Jelovic Ana M., Mitic Nenad S., Eshafah Samira, and Beljanski Milos V.. Journal of Computational Biology. December 2017, ahead of print. https://doi.org/10.1089/cmb.2017.0046

===============================================================================================================================================================

StatRepeats is a tool for finding maximal repeats in an input (nucleotide or protein) sequence. It can extract four types of repeats: direct non-complementary, direct complementary, inverse non-complementary and inverse complementary repeats. It can use statistical estimation to extract only statistically significant repeats, or all repeats. It also can filter the set of found repeats based on the specified p-value. Results can be redirected to the standard output, to files in a few different formats, or directly into an ODBC-compatible database.

Usage:

StatRepeats < input file name > < minimal length > [optional arguments] 
where
Input file name and minimal length are mandatory arguments.
Input file name is the name of the FASTA file to be processed. The file can contain multiple FASTA sequences, which are by default processed one by one.
Minimal length specifies the minimal length of the extracted repeats. If the specified length is to low for statistical filtering to work correctly the program will suggest a minimal length. If zero is specified as the minimal length StatRepeats will automatically choose the lowest possible value for which the probabability estimation used for statistical filtering is valid.
StatRepeats by default can handle collection of multiple fasta sequences in one file. Every sequence is processed separately; program options are applied on all fasta sequences. Calculation was done as series of calculation for each sequence. Output for next calculation is append at the end of previous one.
Example 1: calculate (dn) repeats with minimal length 10  Results 
Note: In this and all further examples, it is assumed that the StatRepeats program is in the same directory as the input data

Optional arguments
File output options
By default StatRepeats will print all results to the standard output.

General output format: [<prefix>],<LP start position>, <LP end position>, <RP start position>, <RP end position>, <sequence length>, <LP>, <RP> 
where LP and RP are the left and right parts of repeat sequence. Prefix in result is taken from fasta header in input file. If fasta header is not present prefix is null (i.e. not exists). Both start and end position of left and right parts of repeat sequence are included in order to make furher processing more efficient.

Options related to file output
-o[utput] filename	Write results to a specified file. Default for output is writing results on console. When -output option is specified, other output options are ignored
-loa[d] name	Produce three output files. Sequence name is written in file name.id, statistics are written in file name.stat, and results are written to a file name.load. If -load name option is combined with multiple fasta sequences in one file, then output for all sequences is placed in three files name.id, name.stat and name.load
-sp[lit]	Produce three output files. Sequence name is written in file named sequencename.id, statistics are written in file named sequencename.stat, and results are written in file named sequencename.load. Input file must have fasta header. In opposite case, generated files with results have no prefix and start with "." (.id, .load, .stat). If -split option is combined with multiple fasta sequences in one file, then output produce, for each sequence in input file, three files with sufixes .id, .stat and .load and prefix equal with id from fasta header
-pri[ntinstanes] <true | false>	Specifies whether to output found pairs of sequences. -pri[ntinstanes] <false> can be used as a quick way to count the number of pairs found. -pri[ntinstanes] <true> is the default option
Example: calculate (dn) repeats and write output to file  Results

Repeat types options
-dn	Search for direct non-complementary repeats
-dc	Search for direct complementary repeats
-in	Search for inverse non-complementary repeats
-ic	Search for inverse complementary repeats
If none of these options is specified the default is -dn.
Example: calculate different types of repeats  Results

Probability estimation
StatRepeats uses probability estimation as explained in paper to estimate if a result is statistically significant. The following options can be specified:
-pv[alue] value	Specify value of p-value parameter. Value is a number between 0 and 1. The default is 0.05
-n[oprobability]	extract all repeat sequences without probability estimation
Example: calculate different types of repeats  Results

Option related to input sequence characteristics
Notes:
StatRepeats can handle an arbitrary alphabet. In the case of complementary repeat search letter complements must be specified. Input alphabet must include even number of characters (even - character and its complement)
If input material includes characters not defined in the alphabet, then number of such character is reported, and each of these characters is assumed to be its own complement.
In options aliphatic, sulphur, tiny, aromatic, hydrophobic, charged, positive, polar, acidic, small, and hydroxylic amino acids are groupped according to Esquivel et al:Decoding the Building Blocks of Life from the Perspective of Quantum Information in "Advances in Quantum Mechanics", Paul Bracken (ed.), INTECH, Rijeka, pp: 641--669, (2013). 
If some of this options is used options -protein must be specified for correct results. Multiple options from this group can be specified only if intersection of their AA groups is empty.
Options:
-dna	Specifies that the input sequence is DNA with the letters ACGT having the complements TGCA. This option is the default
-rna	Specifies that the input sequence is RNA with the letters ACGU having the complements UGCA
-compl[ement]<filename>	Specifies the file with the letters and their complements. The file should be a text file with an even number of letters. The first two letters form a pair with the first letter having the second letter as a complement, then the second two letters form another pair, etc. Statistical filtering is less precise for alphabets that do not consist of four or twenty letters in cases of complementary repeats and is not supported for non-complementary repeats.
-ex[clude] len	Exclude every consecutive len or more letters N from repeat searching. The letter N is genomic sequences is used to specify an unknown nucleid acid. Large islands of N's typically lead to an explosion of found pairs. This option instructs StatRepeats to skip over any island of len or more consecutive N's. If -exclude 0 is specified than all N will be excluded.
-msr	Treat a multi-sequence FASTA file as a single large sequence. Both sequences in the output pair are characterized by the sequence name and intra-sequence position
-max[gap] number	Specifies the maximal gap between two repeats in results. If max=0 then program is search for tandem repeats. If not, specified maximal gap is not set.
-protein	Specifies that the input sequence(s) are proteins. Proteins don't support extraction of complementary repeats as they don't have complements so this option just informs StatRepeats of the alphabet cardinality so that probabilities can be computed. Alphabet cardinality is 20, 21 or 22 depending on the input sequence. Different mappings are also provided when using proteins.
-acidic	Specifies that in the input sequence amino acids N,Q will be mapped into simbol 8. Statistical estimation with new alphabet cardinality is enabled.
-aliphatic	Specifies that in the input sequence amino acids I,V,L will be mapped into simbol 0. Statistical estimation with new alphabet cardinality is enabled.
-aromatic	Specifies that in the input sequence amino acids F,Y,W,H will be mapped into simbol 3. Statistical estimation with new alphabet cardinality is enabled.
-charged	Specifies that in the input sequence amino acids D,E,H,K,R will be mapped into simbol 5. Statistical estimation with new alphabet cardinality is enabled.
-hydrophobic	Specifies that in the input sequence amino acids A,C,T,K,H,W,Y,F,M,I,L,V will be mapped into simbol 4. Statistical estimation with new alphabet cardinality is enabled.
-hydroxylic	Specifies that in the input sequence amino acids S,T will be mapped into simbol @. Statistical estimation with new alphabet cardinality is enabled.
-polar	Specifies that in the input sequence amino acids N,Q,S,T,C,D,E,H,K,R,Y,W will be mapped into simbol 7. Statistical estimation with new alphabet cardinality is enabled.
-positive	Specifies that in the input sequence amino acids H,K,R will be mapped into simbol 6. Statistical estimation with new alphabet cardinality is enabled.
-small	Specifies that in the input sequence amino acids V,P,A,G,C,T,S,D,N will be mapped into simbol 9. Statistical estimation with new alphabet cardinality is enabled.
-sulphur	Specifies that in the input sequence amino acids M,C will be mapped into simbol 1. Statistical estimation with new alphabet cardinality is enabled.
-tiny	Specifies that in the input sequence amino acids A,G,C,S will be mapped into simbol 2. Statistical estimation with new alphabet cardinality is enabled.
Example: Options -dna, -rna, -compl  Results         Options -exclude, -proteins, -msr, -maxgap  Results         Options -proteins -aromatic -aliphatic -charged -tiny  Results


Database Output
Results can be written directly either to DB2 or to other relational database with an ODBC driver installed. A sample database schema can be found in dbschema. Order of SQL files execution is shown in create dbschema. Options are:

-dat[abase] <ODBC connection string>	Specify ODBC connection string used to connect to database
-db2[output] <connection string>	Specify connection string (not necessary ODBC) used to connect to DB2 database
-com[mitcount] number	Specifies how many inserts are done in one transaction. Default is 0 (only one commit is done at the end of the program
-log[ging] true/false	Used only with DB2 row oriented database. Default is true. If false, than alter table ... activate not logged initially is executed in the beginning to the transaction on fragment and match tables
Example: load results directly to relational database  Database schema


Download
Source code of StatRepeats can be download from site Bioinformatics Research Group at Faculty of Mathematics, University of Belgrade. 
Distribution directory includes precompiled version for 64 bit Linux (Ubuntu), 64-bit Windows and license (StatRepeats and Divsufsort) files. Two verson o program are available: with (StatRepeats) or without (StatRepeatsNoDB) interface to relational database. Suffix in program name denotes curent version and release (for example, StatRepeats.v1.r1).

Setup Instructions
Linux
To compile StatRepeats execute the following command in the directory where you unpacked the archive:

     ./StatRepeats.compile.sh for program with interface to relational database

     ./StatRepeatsNoDB.compile.sh for program without interface to relational database

These commands produce the StatRepeats.version.release or StatRepeatsNoDB.version.release binaries. In order for the above command to succeed you will need gcc 4.8.4 or later and libodbc (for version with relational database interface).

In case you are using Debian or a derived distribution the typical way to install ODBC support which is needed to compile StatRepeats is by running the following commands:

    sudo apt-get install odbcinst

    sudo apt-get install unixodbc-dev

For other distributions please check your package manager documentation or build from source.

Windows
StatRepeats executable version (StatRepeats.version.release.exe) is provided in this archive. You will need to install the Visual C++ 2015 runtime in order to run it. It can be downloaded from Microsoft site.

To compile from source you will need the free Visual Studio 2015 Community Edition. Open StatRepeats.sln which is included in the archive and build the solution.