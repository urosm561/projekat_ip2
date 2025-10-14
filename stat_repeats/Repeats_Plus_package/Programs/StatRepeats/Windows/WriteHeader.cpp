#include "WriteHeader.h"

using namespace std;

void WriteHeader(const CommandLine& c, FILE* writer) {
	fprintf(writer, "Input file is %s.\n", c.GetInputFileName().c_str());
	fprintf(writer, "Minimal fragment length is %d.\n", (int)c.GetFragmentLength());
	if (c.GetRna()) {
		fprintf(writer, "Working with rna sequence.\n");
	}
	else if (c.GetDna()) {
		fprintf(writer, "Working with dna sequence.\n");
	}
	else if (c.GetProteins()) {
		fprintf(writer, "Working with protein sequence.\n");
	}

	if (c.GetProteins() && !c.GetNewLetters().empty()) {
		for (const char x: c.GetNewLetters()) {
			switch (x) {
			case '0':
				fprintf(writer, "Simbol 0 is used for aliphatic amino acids.\n");
				break;
			case '1':
				fprintf(writer, "Simbol 1 is used for sulphur amino acids.\n");
				break;
			case '2':
				fprintf(writer, "Simbol 2 is used for tiny amino acids.\n");
				break;
			case '3':
				fprintf(writer, "Simbol 3 is used for aromatic amino acids.\n");
				break;
			case '4':
				fprintf(writer, "Simbol 4 is used for hydrophobic amino acids.\n");
				break;
			case '5':
				fprintf(writer, "Simbol 5 is used for charged amino acids.\n");
				break;
			case '6':
				fprintf(writer, "Simbol 6 is used for positive amino acids.\n");
				break;
			case '7':
				fprintf(writer, "Simbol 7 is used for polar amino acids.\n");
				break;
			case '8':
				fprintf(writer, "Simbol 8 is used for acidic amino acids.\n");
				break;
			case '9':
				fprintf(writer, "Simbol 9 is used for small amino acids.\n");
				break;
			case '@':
				fprintf(writer, "Simbol @ is used for hydroxylic amino acids.\n");
				break;
			}
		}
	}

	if (c.GetProbability()) {
		fprintf(writer, "Working with probability estimation.\n");
		fprintf(writer, "P value is : %f.\n", c.GetPVal());
	}
	else fprintf(writer, "Working without probability estimation.\n");

	if (c.GetPrintInstances()) fprintf(writer, "Printing all instances.\n");
	else fprintf(writer, "Not printing all instances.\n");

	if (c.GetDatabase() != "") fprintf(writer, "Writing to database.\n");

	if (c.GetDb2Output() != "") fprintf(writer, "Writing to database.\n");

	fprintf(writer, "Looking for ");
	c.GetIsRepeat() ? fprintf(writer, "direct ") : fprintf(writer, "inverse ");
	c.GetIsMathematical() ? fprintf(writer, "non-complementary repeats.\n") : fprintf(writer, "complementary repeats.\n");
}

void WriteHeader(const CommandLine& c, ostream& writer)
{
	writer << "Input file is " << c.GetInputFileName() << ". " << '\n';
	writer << "Minimal fragment length is " << c.GetFragmentLength() << ". " << '\n';
	if (c.GetRna()) {
		writer << "Working with rna sequence.\n";
	}
	else if (c.GetDna()) {
		writer << "Working with dna sequence.\n";
	}
	else if (c.GetProteins()) {
		writer << "Working with protein sequence.\n";
	}

	if (c.GetProteins() && !c.GetNewLetters().empty()) {
		for (const char x : c.GetNewLetters()) {
			switch (x) {
			case '0':
				writer << "Simbol 0 is used for aliphatic amino acids.\n";
					break;
			case '1':
				writer << "Simbol 1 is used for sulphur amino acids.\n";
				break;
			case '2':
				writer << "Simbol 2 is used for tiny amino acids.\n";
				break;
			case '3':
				writer << "Simbol 3 is used for aromatic amino acids.\n";
				break;
			case '4':
				writer <<"Simbol 4 is used for hydrophobic amino acids.\n";
				break;
			case '5':
				writer << "Simbol 5 is used for charged amino acids.\n";
				break;
			case '6':
				writer << "Simbol 6 is used for positive amino acids.\n";
				break;
			case '7':
				writer << "Simbol 7 is used for polar amino acids.\n";
				break;
			case '8':
				writer << "Simbol 8 is used for acidic amino acids.\n";
				break;
			case '9':
				writer << "Simbol 9 is used for small amino acids.\n";
				break;
			case '@':
				writer << "Simbol @ is used for hydroxylic amino acids.\n";
				break;
			}
		}
	}

	if (c.GetProbability()) {
		writer << "Working with probabilty estimation." << '\n';
		writer << "P value is : " << c.GetPVal() << "." << '\n';
	}
	else writer << "Working without probability estimation." << '\n';

	if (c.GetPrintInstances()) writer << "Printing all instances." << '\n';
	else writer << "Not printing all instances." << '\n';

	if (c.GetDatabase() != "") writer << "Writing to database." << '\n';

	if (c.GetDb2Output() != "") writer << "Writing to database." << '\n';



	writer << "Looking for "
		<< (c.GetIsRepeat() ? "direct " : "inverse ")
		<< (c.GetIsMathematical() ? "non-complementary repeats." : "complementary repeats.") << '\n';
}
