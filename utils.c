#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

void params_from_file(char *p, char *eA, char *eB, char *PA0, char *PA1,
                      char *QA0, char *QA1, char *RA0, char *RA1, char *PB0,
                      char *PB1, char *QB0, char *QB1, char *RB0, char *RB1,
                      char *file) {
    FILE *fr;
    fr = fopen(file, "rt");
    int lineMAX = 1000;
    char line[lineMAX];

    fgets(line, lineMAX, fr);
    sscanf(line, "%s", p);
    fgets(line, lineMAX, fr);
    sscanf(line, "%s", eA);
    fgets(line, lineMAX, fr);
    sscanf(line, "%s", eB);
    fgets(line, lineMAX, fr);
    sscanf(line, "%s", PA0);
    fgets(line, lineMAX, fr);
    sscanf(line, "%s", PA1);
    fgets(line, lineMAX, fr);
    sscanf(line, "%s", QA0);
    fgets(line, lineMAX, fr);
    sscanf(line, "%s", QA1);
    fgets(line, lineMAX, fr);
    sscanf(line, "%s", RA0);
    fgets(line, lineMAX, fr);
    sscanf(line, "%s", RA1);
    fgets(line, lineMAX, fr);
    sscanf(line, "%s", PB0);
    fgets(line, lineMAX, fr);
    sscanf(line, "%s", PB1);
    fgets(line, lineMAX, fr);
    sscanf(line, "%s", QB0);
    fgets(line, lineMAX, fr);
    sscanf(line, "%s", QB1);
    fgets(line, lineMAX, fr);
    sscanf(line, "%s", RB0);
    fgets(line, lineMAX, fr);
    sscanf(line, "%s", RB1);

    return;
}
