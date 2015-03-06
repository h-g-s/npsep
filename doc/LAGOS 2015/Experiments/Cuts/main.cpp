#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <iostream>
#include <cmath>
#include <cassert>

using namespace std;

void getGaps(const string filename, double *gaps, int maxTime)
{
    FILE *file = fopen(filename.c_str(), "r");
    if(!file)
    {
        perror("Cant open this file!\n");
        exit(EXIT_FAILURE);
    }

    char line[128];
    int rnd, clq;
    double sepTime, obj, opt, gap, compTime = 0.0;
    while(fscanf(file, "%lf %d %d %lf %lf %lf", &sepTime, &rnd, &clq, &obj, &opt, &gap) == 6)
    {
        int intTime;
        double closedGap = (1.0 - gap) * 100.0;
        compTime += sepTime;
        intTime = (int)floor(compTime);
        if(intTime > maxTime)
        	break;
        assert(gaps[intTime] <= closedGap);
        for(int i = intTime; i <= maxTime; i++)
            gaps[i] = closedGap;
    }

    fclose(file);
}

int main(int argc, char const *argv[])
{
    DIR *dp;
    struct dirent *ep;
    int maxTime;

    if(argc != 2)
    {
        perror("Invalid number of argumets. Required: 2.\n");
        exit(EXIT_FAILURE);
    }

    maxTime = atoi(argv[1]); //argument: time limit used in experiments
    dp = opendir("./");

    double avg[maxTime+1];
    int nFiles = 0;
    fill(avg, avg + maxTime + 1, 0.0);

    if(dp != NULL)
    {
        while(ep = readdir (dp))
        {
            string filename = ep->d_name;
            if(filename.find(".log") != string::npos)
            {
                double gaps[maxTime+1];
                fill(gaps, gaps + maxTime + 1, 0.0);
                getGaps(filename, gaps, maxTime);

                for(int i = 0; i <= maxTime; i++)
                    avg[i] += gaps[i];
                nFiles++;
            }
        }
        closedir (dp);
        if(nFiles == 0) return 0;

        avg[0] /= (double)nFiles;
        int prev = (int)floor(avg[0]);
        printf("%d %lf\n", 0, avg[0]);
        for(int i = 1; i < maxTime; i++)
        {
            avg[i] /= (double)nFiles;
            if(fabs(floor(avg[i]) - floor(prev)) <= 0.0001)
            	continue;
            printf("%d %lf\n", i, avg[i]);
            prev = floor(avg[i]);
        }
        avg[maxTime] /= (double)nFiles;
        printf("%d %lf\n", maxTime, avg[maxTime]);
    }
    else
        perror ("Couldn't open the directory");

	return 0;
}