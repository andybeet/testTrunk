/**
\file
\brief The code to interface with the CLAM model.
\ingroup atCLAMLink

    File:		atCLAMLink.c
	Created:	4/9/2008
	Author:		Bec Gorton,   CSIRO Marine Laboratories,  2009
	Purpose:	The code to link to the CLAM model.
    Note:

    Revisions:



    This page is not used but kept for info purposes.

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
#include "atCLAMLinkPrivate.h"

CLAMScenarioOptionStructure CLAMScenarioOptions[] = {
    {"Sea Level Rise", 4, 1},
    {"Climate Change", 5, 1},
    {"New Urban Development", 5, 1},
    {"Agricultural Land Management", 7, 1},
    {"Foreshore Management", 6, 1},
    {"Public Education", 3, 1},
    {"Stormwater Management", 3, 1},
    {"Commercial Fishing Activities", 4, 1},
    {"Dredging", 3, 1},
    {"Aquaculture", 6, 1},
    {"STP Upgrade", 6, 1},
    {"National Park Management", 6, 1},
    {"Flood Mitigation Works", 7, 1}
};


CLAMOutputOptionStructure CLAMOutputOptions[numOutputOptions] = {
    {"Lake Total Nitrogen",1, "-1/2sd, -1/5sd, 0, 1/3sd", 4},
    {"Lake Water Levels",1, "<0.1, 0.1-0.5, >0.5", 3},
    {"Lake Total Phosphorus",1, "-1.2sd, -1/5sd, 0, 1/3sd", 4},
    {"Lake Total Suspended Solids",1, "-1/2sd, -1/5sd, 0, 1/3sd", 4},
    {"Lake Salinity",1, "1/3sd, 0, 1/3sd", 3},
    {"Lake Pathogens",1, "1/2sd, 1/5sd, 0, 1/3sd", 4},
    {"Lake Temperature",1, "-1/3sd, 0, 1/5sd, 1/2sd", 4},
    {"Change in Seagrass Area",1, "decrease, no change, increase", 3},
    {"Change in Saltmarsh Area",1, "decrease, no change, small increase, large increase", 4},
    {"Change in Mangrove Area",1, "decrease, no change, small increase, large increase", 4},
    {"Change in Terrestrial Habitat Area",1, "", 0},
    {"Erosion",1, "", 0},
    {"Change in Agricultural Land Area",1, "", 0},
    {"Social Amenity",1, "", 0},
    {"ASS runoff",1, "", 0},
    {"Fish & Prawn Catch",1, "", 0},
    {"Sedimentation",1, "", 0},
    {"Lake pH",1, "", 0},
    {"Total Seagrass Area",1, "", 0},
    {"Total Saltmarsh Area",1, "", 0},
    {"Total Mangrove Area",1, "", 0},
    {"Total Terrestrial Habitat Area",1, "", 0},
    {"Agricultural Productivity",1, "decrease, no change, increase", 3},
    {"Social Value",1, "", 0},
    {"Lake Algae",1, "", 0},
    {"Total Agricultural Land Area",1, "", 0},
    {"Change in Foreshore Vegetation Area",1, "large decrease, small decrease, no change, small increase, large increase", 5},
    {"Total Foreshore Vegetation Area",1, "", 0},
    {"Lake Dissolved Oxygen",1, "", 0},
    {"Long Term Fish & Prawn Stocks",1, "", 0},
    {"Bird Habitat Value",1, "", 0},
    {"Lake Pesticides",1, "", 0},
    {"Drain Flushing",1, "", 0},
    {"Introduction of Pest Aquatic Species",1, "", 0},
    {"Aquaculture Economy",1, "", 0},
    {"Local Economy",1, "", 0},
    {"Commercial Fishing Economy",1, "", 0},
    {"Aboriginal Values",1, "large decrease, small decrease, no change, small increase, large increase", 5}
};
