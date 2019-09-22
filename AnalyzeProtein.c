/**
 * @file AnalyzeProtein.c
 * @author  Ran Hadar <ran.hadar1@mail.huji.ac.il>
 * @version 1.0
 * @date 31/10/2018
 *
 * @brief Program that reads atoms from a file and calculates equations.
 *
 * @section DESCRIPTION
 * The system keeps track of the cooking times.
 * Input  : pdb files.
 * Process: reads the files responsible for parsing the lines gets the
 * relevants parameters and calculates equations with them.
 * Output : prints the results to the screen.
 */

//************************************  includes ***********************************************
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

//************************************* define **********************************************
#define MAX_ATOMS 20000
#define MAX_LINE 80
#define COORDINATES 3
#define END_OF_STRING '\0'
#define LINE_STARTER "ATOM  "
#define LINE_STARTER_LEN 6
#define COORDINATE_LEN 8 //the max len of the coordinate
#define MIN_LINE_LEN 60
#define MIN_ARGS 2
#define ERROR_ARGS


//*********************************** functions declarations ********************************
int startsWith(char* line);
float getFloat(char* input);
int parseLine(char *fileLine, int atomCounter, float fileAtoms[MAX_ATOMS][COORDINATES]);
void calCenterOfGravity(int atomCounter, float fileAtoms[MAX_ATOMS][COORDINATES], float gravityCenter[3]);
void calRotationRadious(int atomCounter, float fileAtoms[MAX_ATOMS][COORDINATES], float gravityCenter[3]);
void calMaxDistance(int atomCounter, float fileAtoms[MAX_ATOMS][COORDINATES]);

//*******************************************************************************************


/**
 * Gets relevant arguments for running the program.
 * Passing over the files calculates the relevant equations by calling to the relevant functions.
 * Return 0 if the program ran successfuly.
 * @param argc
 * @param argv
 * @return int 0 for success
 */
int main(int argc, char *argv[])
{
    FILE *fp;
    char fileLine[MAX_LINE] = {0};
    float fileAtoms[MAX_ATOMS][COORDINATES] = {{0}};
    int atomCounter;
    float gravityCenter[3];

    if(argc < MIN_ARGS) //checks if no arguments entered
    {
        printf("Usage: AnalyzeProtein <pdb1> <pdb2>");
        return 1;

    }

    for (int i = 1; i < argc; i++)
    {
        atomCounter = 0;
        fp = fopen(argv[i], "r");
        if(fp == NULL) // if the file couldn't open
        {
            printf("Error opening file: %s", argv[i]);
            return 1;
        }

        while(fgets(fileLine, MAX_LINE, fp))
        {

            if (startsWith(fileLine))
            {
                if(strlen(fileLine) > MIN_LINE_LEN)
                {
                    parseLine(fileLine, atomCounter, fileAtoms);
                    atomCounter ++;
                }
                else
                {
                    printf("ATOM line is too short %lu characters", strlen(fileLine));
                    exit(EXIT_FAILURE);
                }

            }

        }
        //checks if problem accoured while reading a line
        if(feof(fp) == 0 || atomCounter == 0)
        {
            printf("Error - 0 atoms were found in the file %s", argv[i]);
            return 1;
        }

        fclose(fp);

        printf("PDB file %s, %d atoms were read\n", argv[i], atomCounter);
        calCenterOfGravity(atomCounter, fileAtoms, gravityCenter);
        calRotationRadious(atomCounter, fileAtoms, gravityCenter);
        calMaxDistance(atomCounter, fileAtoms);

    }

    return 0;
}


/**
 * Gets a string line and checks if the line starts with the substring "ATOM"
 * return 0 for true 1 otherwise
 * @param line
 * @return the result
 */
int startsWith(char* line)
{
    return (strncmp(line, LINE_STARTER, LINE_STARTER_LEN) == 0);

}

/**
 * Gets a string line, number of atoms, and a two dimensional array
 * Responsible for parsing the line, gets the relevant substrings
 * and enters them to the atoms array.
 * @param fileLine
 * @param atomCounter
 * @param fileAtoms
 * @return 0 for success
 */
int parseLine(char *fileLine, int atomCounter, float fileAtoms[MAX_ATOMS][COORDINATES])
{
    char cgString[COORDINATE_LEN + 1], radString[COORDINATE_LEN + 1], disString[COORDINATE_LEN + 1];


    memcpy(cgString, fileLine + 30, COORDINATE_LEN);
    cgString[COORDINATE_LEN] = END_OF_STRING;

    memcpy(radString, fileLine + 38, COORDINATE_LEN);
    radString[COORDINATE_LEN] = END_OF_STRING;

    memcpy(disString, fileLine + 46, COORDINATE_LEN);
    disString[COORDINATE_LEN] = END_OF_STRING;

    fileAtoms[atomCounter][0] = getFloat(cgString);
    fileAtoms[atomCounter][1] = getFloat(radString);
    fileAtoms[atomCounter][2] = getFloat(disString);


    return 0;
}


/**
 * Gets a pointer to the input
 * Responsible for converting the input value to a float.
 * @param input
 * @return the result of the float.
 */
float getFloat(char* input)
{
    char *end;
    float result = 0;
    errno = 0;
    result = strtof(input, &end);
    if(result == 0 && (errno != 0 || end == input))
    {
        fprintf(stderr, "Error in coordinate conversion %s", input);
        exit(EXIT_FAILURE);
    }
    return result;
}


/**
 * Gets the number of atoms, two dimensional atoms array ,empty array of gravity
 * center with 3 coordinates.
 * Responsible for calculation of the gravity center, and fills the results in
 * the initialize gravity center array.
 * @param atomCounter
 * @param fileAtoms
 * @param gravityCenter
 */
void calCenterOfGravity(int atomCounter, float fileAtoms[MAX_ATOMS][COORDINATES], float gravityCenter[3])
{
    float xCorSum, yCorSum, zCorSum;
    xCorSum = yCorSum = zCorSum = 0.0;
    for (int i = 0; i < atomCounter; i++)
    {
        xCorSum += fileAtoms[i][0];
        yCorSum += fileAtoms[i][1];
        zCorSum += fileAtoms[i][2];

    }
    printf("Cg = %.3f %.3f %.3f\n", (xCorSum / atomCounter), (yCorSum / atomCounter), (zCorSum / atomCounter));
    gravityCenter[0] = (xCorSum / atomCounter);
    gravityCenter[1] = (yCorSum / atomCounter);
    gravityCenter[2] = (zCorSum / atomCounter);

}


/**
 * Gets 6 coordinates.
 * Responsible for calculation of the distance between them.
 * Return float distance result.
 * @param xCor1
 * @param yCor1
 * @param zCor1
 * @param xCor2
 * @param yCor2
 * @param zCor2
 * @return result
 */
float calDistance(float xCor1, float yCor1, float zCor1, float xCor2, float yCor2, float zCor2)
{
    return ((xCor1-xCor2)*(xCor1-xCor2) + ((yCor1-yCor2) * (yCor1-yCor2))
                  + ((zCor1-zCor2) * (zCor1-zCor2)));
}

/**
 * Gets number of atoms, two dimensional atoms array, an array
 * represents the gravity center.
 * Responsible for calculation of the distance og each atom in the array
 * comparing to the gravity center coordinadtes.
 * Prints the rotation radious that was found by formula.
 * @param atomCounter
 * @param fileAtoms
 * @param gravityCenter
 */
void calRotationRadious(int atomCounter, float fileAtoms[MAX_ATOMS][COORDINATES], float gravityCenter[3])
{
    float rotationSum = 0.0;
    float xCor1 = gravityCenter[0];
    float yCor1 = gravityCenter[1];
    float zCor1 = gravityCenter[2];
    float tempdis;

    for(int i = 0; i < atomCounter ; i++)
    {
        tempdis = calDistance(xCor1, yCor1, zCor1, fileAtoms[i][0], fileAtoms[i][1], fileAtoms[i][2]);
        rotationSum += tempdis;

    }
    rotationSum = sqrt(rotationSum / atomCounter);
    printf("Rg = %.3f\n", rotationSum);

}

/**
 *  Gets number of atoms, two dimensional atoms array.
 * Responsible for calculation of the max distance between two atoms
 * by 3 coordinates.
 * Prints the max distance.
 * @param atomCounter
 * @param fileAtoms
 */
void calMaxDistance(int atomCounter, float fileAtoms[MAX_ATOMS][COORDINATES])
{
    float maxDistance = 0.0;
    float xCor1, yCor1, zCor1;
    float tempDistance;

    //Passing all over the atoms in the array
    for(int i = 0; i < atomCounter -1; i++)
    {
        xCor1 = fileAtoms[i][0];
        yCor1 = fileAtoms[i][1];
        zCor1 = fileAtoms[i][2];
        for(int j = 1; j < atomCounter; j++)
        {
            tempDistance = sqrt(calDistance(xCor1, yCor1, zCor1, fileAtoms[j][0], fileAtoms[j][1], fileAtoms[j][2]));
            if(maxDistance < tempDistance)
            {//checks if the distance is bigger
                maxDistance = tempDistance;
            }
        }

    }
    printf("Dmax = %.3f\n", maxDistance);
}

