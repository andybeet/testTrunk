/**************************************************************************//**
 \file
 \brief The code to interface with the CLAM model.
 \ingroup atCLAMLink

 File:       atCLAMLink.c
 Created:    4/9/2008
 Author:     Bec Gorton,   CSIRO Marine Laboratories,  2009
 Purpose:    The code to link to the CLAM model.
 Note:

 Revisions:

 *************/
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <sjwlib.h>
#include <atlantisboxmodel.h>
#include <atUtilLib.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include "atCLAMLinkPrivate.h"

static void Create_CLAM_Input_File(MSEBoxModel *bm, char *inputFileName);

/**
 * \brief Routine to start the CLAM Process.
 *
 * This will create the CLAM input file containing information about which scenarios to run and which output
 * values to generate.
 * If our operating system is windows (if _WIN32 is defined in the preprocessor) then we will do the actual
 * call to the CLAM by starting a new thread.
 * If we are not in windows then it will just generate the CLAM input file for test purposes.
 *
 */
void Start_CLAM(MSEBoxModel *bm, char *inputFileName) {
#ifdef _WIN32
	FILE *fp;

	STARTUPINFO si;
	PROCESS_INFORMATION pi;
#endif

	/* Might need to delete path contents for each call */

	/* Create the CLAM input file */
	Create_CLAM_Input_File(bm, inputFileName);

#ifdef _WIN32

	ZeroMemory( &si, sizeof(si) );
	si.cb = sizeof(si);
	ZeroMemory( &pi, sizeof(pi) );

	printf("atCLAMLink: Creating CLAM process.\n");

	/* Start the child process. */
	if( !CreateProcess( NULL, // No module name (use command line)
					commandLinesStr, // Command line
					NULL, // Process handle not inheritable
					NULL, // Thread handle not inheritable
					FALSE, // Set handle inheritance to FALSE
					0, // No creation flags
					NULL, // Use parent's environment block
					NULL, // Use parent's starting directory
					&si, // Pointer to STARTUPINFO structure
					&pi ) // Pointer to PROCESS_INFORMATION structure
	) {
		printf( "CreateProcess failed (%d)\n", GetLastError() );
		return;
	}

	/* Add a check to see that the CLAM model has started */
	/* Want to wait a small amount of time - check that the 'started' file has been created - if not then quit */

	/* Wait for a little period of time to see if the 'started' file has been created */
	WaitForSingleObject( pi.hProcess, startedWaitTime);

	/* Check to see that the file has been created */
	if( (fp = Open_Input_File(bm->inputFolder,startedFileName, "r")) == NULL) {
		/* Terminate the process */
		TerminateProcess(pi.hProcess, 0);
		quit("CLAM Model did not start\n");
	}
	fclose(fp);

	printf("atCLAMLink: CLAM process successfully started.\n");

	/* Wait until child process exits. */
	WaitForSingleObject( pi.hProcess, processWaitTime );

	printf("atCLAMLink: CLAM process finished - checking to see if that was sucessful.\n");

	/* Now see if there is a file called error.txt*/
	if( (fp = Open_Input_File(bm->inputFolder,errorFileName, "r")) != NULL) {
		/* Terminate the process */
		TerminateProcess(pi.hProcess, 0);
		fclose(fp);
		quit("There was an error with the CLAM model. Please check the error file %s%s for more information.\n", bm->inputFolder, errorFileName);
	}

	printf("atCLAMLink: No CLAM error file found. Thats a good thing!\n");

	/* Now see if there is a file called finished.txt*/
	if( (fp = Open_Input_File(bm->inputFolder,finishedFileName, "r")) == NULL) {
		/* Terminate the process */
		TerminateProcess(pi.hProcess, 0);
		quit("CLAM Model did not finish successfully. Please check the error file %s%s.\n", bm->inputFolder, errorFileName);
	}
	fclose(fp);

	printf("atCLAMLink: CLAM finished file found.\n");

	printf("atCLAMLink: Finished reading in the CLAM output data.\n");

	/* Close process and thread handles. */
	CloseHandle( pi.hProcess );
	CloseHandle( pi.hThread );

	printf("atCLAMLink: All Done. Freeing memory and continuing with the Atlantis model run.\n");
#endif
}

/**
 *  /brief Create the CLAM input file that specifies the management scenario and output options.
 *
 *
 */
static void Create_CLAM_Input_File(MSEBoxModel *bm, char *inputFileName) {
	int i;
	FILE *fp;

	/* Create file */
	if ((fp = fopen(inputFileName, "w")) == NULL)
		quit("Create_CLAM_Input_File: Can't open input file '%s'\n", inputFileName);

	/* Write out information about the scenarios to run */
	fprintf(fp, "INPUT RULES\n");
	for (i = 0; i < numScenarioOptions; i++) {
		if (CLAMScenarioOptions[i].optionSelected > CLAMScenarioOptions[i].maxOptionValue)
			quit("Create_CLAM_Input_File The CLAM Scenario option for %s is greater than the allowed value. You selected %d but the max index is %d\n",
					CLAMScenarioOptions[i].scenarioOptionName, CLAMScenarioOptions[i].optionSelected, CLAMScenarioOptions[i].maxOptionValue);
		fprintf(fp, "%s,%d\n", CLAMScenarioOptions[i].scenarioOptionName, CLAMScenarioOptions[i].optionSelected + 1);
		//fprintf(fp, "%s,%d\n", CLAMScenarioOptions[i].scenarioOptionName, 3);
	}

	/* Tell the CLAM which indicators we are interest in */
	fprintf(fp, "OUTPUT RULES\n");
	for (i = 0; i < numClamIndicators; i++) {
		if (CLAMOutputOptions[i].optionSelected != 0 && CLAMOutputOptions[i].optionSelected != 1)
			quit("Create_CLAM_Input_File The CLAM Output option for %s is greater than the allowed value. You selected %d but the value should be 0 or 1.\n",
					CLAMOutputOptions[i].outputOptionName, CLAMOutputOptions[i].optionSelected);
		fprintf(fp, "%s,%d\n", CLAMOutputOptions[i].outputOptionName, CLAMOutputOptions[i].optionSelected);
	}
	fclose(fp);

	printf("atCLAMLink: CLAM input file sucessfully created.\n");

}

/**
 *  \brief Create the ini file that is used to init the CLAM model.
 *
 *  This file contains things like where the input csv file is,
 *  the path of the folder to store the output file.
 *
 *
 */
void Create_CLAM_Ini_File(MSEBoxModel *bm) {
	FILE *fp;

	/* Create file */
	if ((fp = fopen(iniFileName, "w")) == NULL)
		quit("Create_CLAM_Ini_File: Can't open ini file - '%s'\n", iniFileName);

	fprintf(fp, "[autoPluginCLAM]\n");

	fprintf(fp, ";Folder to write the output to (do not add final )\n");
	fprintf(fp, "where_to_write=%s\n", outputFolderPath);

	fprintf(fp, ";File that drives the automation\n");
	fprintf(fp, "driver_file=%s\n", inputFileName);

	fclose(fp);

	printf("atCLAMLink: CLAM ini file sucessfully created.\n");
}

/**
 *  \brief Parse the CLAM output file.
 *
 *
 *
 *
 */
void Parse_CLAM_OutputFile(MSEBoxModel *bm, int scenarioIndex, int strategyIndex) {
	char ch, buf[200], seps[] = ",";
	char *valueStr;
	int buflen = 200;
	int i;
	FILE *infile;
	int indicatorIndex;

	/* Create file */
	if ((infile = Open_Input_File(bm->inputFolder,outputFileName, "r")) == NULL)
		quit("createCLAMiniFile: Can't open %s%s\n", bm->inputFolder, outputFileName);

	/* Loop for all lines */
	fseek(infile, 0L, 0);
	while (fgets(buf, buflen, infile) != NULL) {
		ch = buf[0];

		if (ch == '\n')
			continue;

		printf("buf = %s\n", buf);
		valueStr = strtok(buf, seps);
		/* Ignore teh first token */
		valueStr = strtok(NULL, seps);

		/* Strip off the text after the '.' */
		*(strstr(valueStr, ".")) = '\0';
		printf("Output Parameter = %s\n", valueStr);

		/* Find the index of this output parameter  */
		indicatorIndex = Get_CLAM_IndicatorID(bm, valueStr);

		// indicatorIndex = 1;

		/* ignore the next token */
		//valueStr = strtok(NULL, seps);

		/* Ignore indicators we are not interested in */
		if (indicatorIndex != -1) {
			printf("indicatorIndex = %d\n", indicatorIndex);
			fprintf(bm->logFile, "Output Parameter = %s, indicatorIndex = %d, scenarioIndex = %d, strategyIndex = %d\n", valueStr, indicatorIndex,
					scenarioIndex, strategyIndex);

			i = 0;
			while ((valueStr = strtok(NULL, seps)) != NULL) {
				outcomeValueArray[indicatorIndex][scenarioIndex][strategyIndex][i] = atof(valueStr);
				printf("outcomeValueArray[%d][%d][%d][%d] = %e, \n", indicatorIndex, scenarioIndex, strategyIndex, i,
						outcomeValueArray[indicatorIndex][scenarioIndex][strategyIndex][i]);
				fprintf(bm->logFile, "outcomeValueArray[%d][%d][%d][%d] = %e, \n", indicatorIndex, scenarioIndex, strategyIndex, i,
						outcomeValueArray[indicatorIndex][scenarioIndex][strategyIndex][i]);

				i++;
			}

			printf("\n\n");
		}
		/*else{
		 quit("Could not find indicators %s\n", valueStr);
		 }*/
	}
	fclose(infile);
}

