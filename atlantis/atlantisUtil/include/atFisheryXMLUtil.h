/**
 * \ingroup atUtil
 * \file atFisheryXMLUtil.h
 *
 *  Created on: 01/02/2010
 *      Author: bec
 */

#ifndef ATFISHERYXMLUTIL_H_
#define ATFISHERYXMLUTIL_H_

#define numFisheryParams 79
#define numGroupFisheryParams 44
extern SpeciesParamStruct FisheryParamsArray[numFisheryParams];
extern SpeciesParamStruct FisheryGroupParamsArray[numGroupFisheryParams];

int Util_XML_Get_Fishery_Param_Index(SpeciesParamStruct array[], int size, int paramID);
void Util_XML_Read_Fishery_Param(MSEBoxModel *bm, char *fileName, xmlNodePtr parent, int paramID);
void Util_XML_Read_Fishery_Param_Double(MSEBoxModel *bm, char *fileName, xmlNodePtr parent, char *paramName, int checkType, double **targetArray, int size);
void Util_XML_Read_Fishery_Param_Int(MSEBoxModel *bm, char *fileName, xmlNodePtr parent, char *paramName, int checkType, int **targetArray, int size);
void Util_XML_Read_Fished_Group_Param(MSEBoxModel *bm, char *fileName, xmlNodePtr parent, char *paramName, int check_type, double **targetArray, int size);

void Util_XML_Read_Impacted_Group_Param(MSEBoxModel *bm, char *fileName, xmlNodePtr parent, int paramID);
void Util_XML_Read_Fishery_Group_Param(MSEBoxModel *bm, char *fileName, xmlNodePtr parent, int paramID);
void Util_XML_Read_Change_Values(MSEBoxModel *bm, char *fileName, xmlNodePtr attributeGroup, char *paramName, int paramID, double ***array);
void Util_XML_Read_Change_Fished_Groups(MSEBoxModel *bm, char *fileName, xmlNodePtr attributeGroup, char *paramName, int paramID, double ****array);
void Util_XML_Read_Vert_Fishery_Param(MSEBoxModel *bm, char *fileName, xmlNodePtr parent, int paramID);

#endif /* ATFISHERYXMLUTIL_H_ */
