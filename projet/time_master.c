
int main(int argc, char const *argv[])
{
	char * man = "Give me a valid path, nega #WestCoast"
	if (argc<2){
		printf("%s\n", man);
		return 0;
	}
	unsigned long resultat = 0;

	FILE * fp = fopen(argv[1], "r");
	if (fp!=NULL)
	{
		unsigned long tmp;
		fseek(fp, 0, SEEK_END);
		long lSize = ftell(fp);
		rewind(fp);
		while(fread(&tmp, sizeof(char))){

		}


	}

	fclose(fp);
	return 0;
}