/**************** gf_imp_xyzt.cc ****************************/
/*
 *
 * Calculate symzatic improved gauge force on each site
 * including rectangulars and bentchairs.
 *
 * Update momentum fields.
 *
 */
/********************************************************/

#include "gf_imp_xyzt.h"

/*
static string inHHBase[8] = { "inHxback", "inHxforw", 
			     "inHyback", "inHyforw", 
			     "inHzback", "inHzforw", 
			     "inHtback", "inHtforw" };
*/
static string inGBase[8] = { "inGxback", "inGxforw", 
                             "inGyback", "inGyforw", 
                             "inGzback", "inGzforw", 
                             "inGtback", "inGtforw" };
static string outHBase[8] = { "outHxback", "outHxforw", 
                               "outHyback", "outHyforw", 
                               "outHzback", "outHzforw", 
                               "outHtback", "outHtforw" };
static string inHBase("inHBase");
static string iprefHBD("hiprefdist");
static string prefHBD[4][2] = { { "hxyBase", "hprefdist1" },
				{ "hpfBase2", "hprefdist2" },
				{ "hpfBase3", "hprefdist3" },
				{ "hpfBase4", "hprefdist4" } };

static string prefGBD[4][2] = { { "gxyBase", "gprefdist1" },
    				{ "gpfBase2", "gprefdist2" },
    				{ "gpfBase3", "gprefdist3" },
    				{ "gpfBase4", "gprefdist4" } };
static string epsilonH("epsilonH");
static string half_epsilonH("(epsilonH/2.0)");
static string epsilonU("epsilonU");
static string kappaR("kappaR");
static string kappaS("kappaS");
static string kappaB("kappaB");
string mask;
static FVec tmpMtr[3][3][2] = { { {FVec("tmpMtr_00RE"), FVec("tmpMtr_00IM")}, 
				  {FVec("tmpMtr_01RE"), FVec("tmpMtr_01IM")}, 
				  {FVec("tmpMtr_02RE"), FVec("tmpMtr_02IM")} },
				{ {FVec("tmpMtr_10RE"), FVec("tmpMtr_10IM")},  
                                  {FVec("tmpMtr_11RE"), FVec("tmpMtr_11IM")},
                                  {FVec("tmpMtr_12RE"), FVec("tmpMtr_12IM")} },
				{ {FVec("tmpMtr_20RE"), FVec("tmpMtr_20IM")},  
                                  {FVec("tmpMtr_21RE"), FVec("tmpMtr_21IM")},
                                  {FVec("tmpMtr_22RE"), FVec("tmpMtr_22IM")} } }; 

static FVec tmpMtr2[3][3][2] = { { {FVec("tmpMtr2_00RE"), FVec("tmpMtr2_00IM")},  
                                  {FVec("tmpMtr2_01RE"), FVec("tmpMtr2_01IM")},
                                  {FVec("tmpMtr2_02RE"), FVec("tmpMtr2_02IM")} },
                                { {FVec("tmpMtr2_10RE"), FVec("tmpMtr2_10IM")},
                                  {FVec("tmpMtr2_11RE"), FVec("tmpMtr2_11IM")},
                                  {FVec("tmpMtr2_12RE"), FVec("tmpMtr2_12IM")} },
                                { {FVec("tmpMtr2_20RE"), FVec("tmpMtr2_20IM")},
                                  {FVec("tmpMtr2_21RE"), FVec("tmpMtr2_21IM")},
                                  {FVec("tmpMtr2_22RE"), FVec("tmpMtr2_22IM")} } };

FVec tmpMtr3[3][3][2] = { { {FVec("tmpMtr3_00RE"), FVec("tmpMtr3_00IM")},
                                  {FVec("tmpMtr3_01RE"), FVec("tmpMtr3_01IM")},
                                  {FVec("tmpMtr3_02RE"), FVec("tmpMtr3_02IM")} },
                                { {FVec("tmpMtr3_10RE"), FVec("tmpMtr3_10IM")},
                                  {FVec("tmpMtr3_11RE"), FVec("tmpMtr3_11IM")},
                                  {FVec("tmpMtr3_12RE"), FVec("tmpMtr3_12IM")} },
                                { {FVec("tmpMtr3_20RE"), FVec("tmpMtr3_20IM")},
                                  {FVec("tmpMtr3_21RE"), FVec("tmpMtr3_21IM")},
                                  {FVec("tmpMtr3_22RE"), FVec("tmpMtr3_22IM")} } };

static FVec Out1[3][3][2] = { { {FVec("Out1_00RE"), FVec("Out1_00IM")},
                                  {FVec("Out1_01RE"), FVec("Out1_01IM")},
                                  {FVec("Out1_02RE"), FVec("Out1_02IM")} },
                                { {FVec("Out1_10RE"), FVec("Out1_10IM")},
                                  {FVec("Out1_11RE"), FVec("Out1_11IM")},
                                  {FVec("Out1_12RE"), FVec("Out1_12IM")} },
                                { {FVec("Out1_20RE"), FVec("Out1_20IM")},
                                  {FVec("Out1_21RE"), FVec("Out1_21IM")},
                                  {FVec("Out1_22RE"), FVec("Out1_22IM")} } };

static FVec Out2[3][3][2] = { { {FVec("Out2_00RE"), FVec("Out2_00IM")},
                                  {FVec("Out2_01RE"), FVec("Out2_01IM")},
                                  {FVec("Out2_02RE"), FVec("Out2_02IM")} },
                                { {FVec("Out2_10RE"), FVec("Out2_10IM")},
                                  {FVec("Out2_11RE"), FVec("Out2_11IM")},
                                  {FVec("Out2_12RE"), FVec("Out2_12IM")} },
                                { {FVec("Out2_20RE"), FVec("Out2_20IM")},
                                  {FVec("Out2_21RE"), FVec("Out2_21IM")},
                                  {FVec("Out2_22RE"), FVec("Out2_22IM")} } };

static FVec staple1[6][16][3][3][2] = { { { { {FVec("staple1xy00_00RE"), FVec("staple1xy00_00IM")},
                                  {FVec("staple1xy00_01RE"), FVec("staple1xy00_01IM")},
                                  {FVec("staple1xy00_02RE"), FVec("staple1xy00_02IM")} },
                                { {FVec("staple1xy00_10RE"), FVec("staple1xy00_10IM")},
                                  {FVec("staple1xy00_11RE"), FVec("staple1xy00_11IM")},
                                  {FVec("staple1xy00_12RE"), FVec("staple1xy00_12IM")} },
                                { {FVec("staple1xy00_20RE"), FVec("staple1xy00_20IM")},
                                  {FVec("staple1xy00_21RE"), FVec("staple1xy00_21IM")},
                                  {FVec("staple1xy00_22RE"), FVec("staple1xy00_22IM")} } }, 
					{ { {FVec("staple1xy01_00RE"), FVec("staple1xy01_00IM")},
                                  {FVec("staple1xy01_01RE"), FVec("staple1xy01_01IM")},
                                  {FVec("staple1xy01_02RE"), FVec("staple1xy01_02IM")} },
                                { {FVec("staple1xy01_10RE"), FVec("staple1xy01_10IM")},
                                  {FVec("staple1xy01_11RE"), FVec("staple1xy01_11IM")},
                                  {FVec("staple1xy01_12RE"), FVec("staple1xy01_12IM")} },
                                { {FVec("staple1xy01_20RE"), FVec("staple1xy01_20IM")},
                                  {FVec("staple1xy01_21RE"), FVec("staple1xy01_21IM")},
                                  {FVec("staple1xy01_22RE"), FVec("staple1xy01_22IM")} } }, 
					{ { {FVec("staple1xy02_00RE"), FVec("staple1xy02_00IM")},
                                  {FVec("staple1xy02_01RE"), FVec("staple1xy02_01IM")},
                                  {FVec("staple1xy02_02RE"), FVec("staple1xy02_02IM")} },
                                { {FVec("staple1xy02_10RE"), FVec("staple1xy02_10IM")},
                                  {FVec("staple1xy02_11RE"), FVec("staple1xy02_11IM")},
                                  {FVec("staple1xy02_12RE"), FVec("staple1xy02_12IM")} },
                                { {FVec("staple1xy02_20RE"), FVec("staple1xy02_20IM")},
                                  {FVec("staple1xy02_21RE"), FVec("staple1xy02_21IM")},
                                  {FVec("staple1xy02_22RE"), FVec("staple1xy02_22IM")} } }, 
					{ { {FVec("staple1xy03_00RE"), FVec("staple1xy03_00IM")},
                                  {FVec("staple1xy03_01RE"), FVec("staple1xy03_01IM")},
                                  {FVec("staple1xy03_02RE"), FVec("staple1xy03_02IM")} },
                                { {FVec("staple1xy03_10RE"), FVec("staple1xy03_10IM")},
                                  {FVec("staple1xy03_11RE"), FVec("staple1xy03_11IM")},
                                  {FVec("staple1xy03_12RE"), FVec("staple1xy03_12IM")} },
                                { {FVec("staple1xy03_20RE"), FVec("staple1xy03_20IM")},
                                  {FVec("staple1xy03_21RE"), FVec("staple1xy03_21IM")},
                                  {FVec("staple1xy03_22RE"), FVec("staple1xy03_22IM")} } }, 
					{ { {FVec("staple1xy04_00RE"), FVec("staple1xy04_00IM")},
                                  {FVec("staple1xy04_01RE"), FVec("staple1xy04_01IM")},
                                  {FVec("staple1xy04_02RE"), FVec("staple1xy04_02IM")} },
                                { {FVec("staple1xy04_10RE"), FVec("staple1xy04_10IM")},
                                  {FVec("staple1xy04_11RE"), FVec("staple1xy04_11IM")},
                                  {FVec("staple1xy04_12RE"), FVec("staple1xy04_12IM")} },
                                { {FVec("staple1xy04_20RE"), FVec("staple1xy04_20IM")},
                                  {FVec("staple1xy04_21RE"), FVec("staple1xy04_21IM")},
                                  {FVec("staple1xy04_22RE"), FVec("staple1xy04_22IM")} } }, 
					{ { {FVec("staple1xy05_00RE"), FVec("staple1xy05_00IM")},
                                  {FVec("staple1xy05_01RE"), FVec("staple1xy05_01IM")},
                                  {FVec("staple1xy05_02RE"), FVec("staple1xy05_02IM")} },
                                { {FVec("staple1xy05_10RE"), FVec("staple1xy05_10IM")},
                                  {FVec("staple1xy05_11RE"), FVec("staple1xy05_11IM")},
                                  {FVec("staple1xy05_12RE"), FVec("staple1xy05_12IM")} },
                                { {FVec("staple1xy05_20RE"), FVec("staple1xy05_20IM")},
                                  {FVec("staple1xy05_21RE"), FVec("staple1xy05_21IM")},
                                  {FVec("staple1xy05_22RE"), FVec("staple1xy05_22IM")} } }, 
					{ { {FVec("staple1xy06_00RE"), FVec("staple1xy06_00IM")},
                                  {FVec("staple1xy06_01RE"), FVec("staple1xy06_01IM")},
                                  {FVec("staple1xy06_02RE"), FVec("staple1xy06_02IM")} },
                                { {FVec("staple1xy06_10RE"), FVec("staple1xy06_10IM")},
                                  {FVec("staple1xy06_11RE"), FVec("staple1xy06_11IM")},
                                  {FVec("staple1xy06_12RE"), FVec("staple1xy06_12IM")} },
                                { {FVec("staple1xy06_20RE"), FVec("staple1xy06_20IM")},
                                  {FVec("staple1xy06_21RE"), FVec("staple1xy06_21IM")},
                                  {FVec("staple1xy06_22RE"), FVec("staple1xy06_22IM")} } }, 
					{ { {FVec("staple1xy07_00RE"), FVec("staple1xy07_00IM")},
                                  {FVec("staple1xy07_01RE"), FVec("staple1xy07_01IM")},
                                  {FVec("staple1xy07_02RE"), FVec("staple1xy07_02IM")} },
                                { {FVec("staple1xy07_10RE"), FVec("staple1xy07_10IM")},
                                  {FVec("staple1xy07_11RE"), FVec("staple1xy07_11IM")},
                                  {FVec("staple1xy07_12RE"), FVec("staple1xy07_12IM")} },
                                { {FVec("staple1xy07_20RE"), FVec("staple1xy07_20IM")},
                                  {FVec("staple1xy07_21RE"), FVec("staple1xy07_21IM")},
                                  {FVec("staple1xy07_22RE"), FVec("staple1xy07_22IM")} } }, 
					{ { {FVec("staple1xy10_00RE"), FVec("staple1xy10_00IM")},
                                  {FVec("staple1xy10_01RE"), FVec("staple1xy10_01IM")},
                                  {FVec("staple1xy10_02RE"), FVec("staple1xy10_02IM")} },
                                { {FVec("staple1xy10_10RE"), FVec("staple1xy10_10IM")},
                                  {FVec("staple1xy10_11RE"), FVec("staple1xy10_11IM")},
                                  {FVec("staple1xy10_12RE"), FVec("staple1xy10_12IM")} },
                                { {FVec("staple1xy10_20RE"), FVec("staple1xy10_20IM")},
                                  {FVec("staple1xy10_21RE"), FVec("staple1xy10_21IM")},
                                  {FVec("staple1xy10_22RE"), FVec("staple1xy10_22IM")} } }, 
					{ { {FVec("staple1xy11_00RE"), FVec("staple1xy11_00IM")},
                                  {FVec("staple1xy11_01RE"), FVec("staple1xy11_01IM")},
                                  {FVec("staple1xy11_02RE"), FVec("staple1xy11_02IM")} },
                                { {FVec("staple1xy11_10RE"), FVec("staple1xy11_10IM")},
                                  {FVec("staple1xy11_11RE"), FVec("staple1xy11_11IM")},
                                  {FVec("staple1xy11_12RE"), FVec("staple1xy11_12IM")} },
                                { {FVec("staple1xy11_20RE"), FVec("staple1xy11_20IM")},
                                  {FVec("staple1xy11_21RE"), FVec("staple1xy11_21IM")},
                                  {FVec("staple1xy11_22RE"), FVec("staple1xy11_22IM")} } }, 
					{ { {FVec("staple1xy12_00RE"), FVec("staple1xy12_00IM")},
                                  {FVec("staple1xy12_01RE"), FVec("staple1xy12_01IM")},
                                  {FVec("staple1xy12_02RE"), FVec("staple1xy12_02IM")} },
                                { {FVec("staple1xy12_10RE"), FVec("staple1xy12_10IM")},
                                  {FVec("staple1xy12_11RE"), FVec("staple1xy12_11IM")},
                                  {FVec("staple1xy12_12RE"), FVec("staple1xy12_12IM")} },
                                { {FVec("staple1xy12_20RE"), FVec("staple1xy12_20IM")},
                                  {FVec("staple1xy12_21RE"), FVec("staple1xy12_21IM")},
                                  {FVec("staple1xy12_22RE"), FVec("staple1xy12_22IM")} } }, 
					{ { {FVec("staple1xy13_00RE"), FVec("staple1xy13_00IM")},
                                  {FVec("staple1xy13_01RE"), FVec("staple1xy13_01IM")},
                                  {FVec("staple1xy13_02RE"), FVec("staple1xy13_02IM")} },
                                { {FVec("staple1xy13_10RE"), FVec("staple1xy13_10IM")},
                                  {FVec("staple1xy13_11RE"), FVec("staple1xy13_11IM")},
                                  {FVec("staple1xy13_12RE"), FVec("staple1xy13_12IM")} },
                                { {FVec("staple1xy13_20RE"), FVec("staple1xy13_20IM")},
                                  {FVec("staple1xy13_21RE"), FVec("staple1xy13_21IM")},
                                  {FVec("staple1xy13_22RE"), FVec("staple1xy13_22IM")} } }, 
					{ { {FVec("staple1xy14_00RE"), FVec("staple1xy14_00IM")},
                                  {FVec("staple1xy14_01RE"), FVec("staple1xy14_01IM")},
                                  {FVec("staple1xy14_02RE"), FVec("staple1xy14_02IM")} },
                                { {FVec("staple1xy14_10RE"), FVec("staple1xy14_10IM")},
                                  {FVec("staple1xy14_11RE"), FVec("staple1xy14_11IM")},
                                  {FVec("staple1xy14_12RE"), FVec("staple1xy14_12IM")} },
                                { {FVec("staple1xy14_20RE"), FVec("staple1xy14_20IM")},
                                  {FVec("staple1xy14_21RE"), FVec("staple1xy14_21IM")},
                                  {FVec("staple1xy14_22RE"), FVec("staple1xy14_22IM")} } }, 
					{ { {FVec("staple1xy15_00RE"), FVec("staple1xy15_00IM")},
                                  {FVec("staple1xy15_01RE"), FVec("staple1xy15_01IM")},
                                  {FVec("staple1xy15_02RE"), FVec("staple1xy15_02IM")} },
                                { {FVec("staple1xy15_10RE"), FVec("staple1xy15_10IM")},
                                  {FVec("staple1xy15_11RE"), FVec("staple1xy15_11IM")},
                                  {FVec("staple1xy15_12RE"), FVec("staple1xy15_12IM")} },
                                { {FVec("staple1xy15_20RE"), FVec("staple1xy15_20IM")},
                                  {FVec("staple1xy15_21RE"), FVec("staple1xy15_21IM")},
                                  {FVec("staple1xy15_22RE"), FVec("staple1xy15_22IM")} } }, 
					{ { {FVec("staple1xy16_00RE"), FVec("staple1xy16_00IM")},
                                  {FVec("staple1xy16_01RE"), FVec("staple1xy16_01IM")},
                                  {FVec("staple1xy16_02RE"), FVec("staple1xy16_02IM")} },
                                { {FVec("staple1xy16_10RE"), FVec("staple1xy16_10IM")},
                                  {FVec("staple1xy16_11RE"), FVec("staple1xy16_11IM")},
                                  {FVec("staple1xy16_12RE"), FVec("staple1xy16_12IM")} },
                                { {FVec("staple1xy16_20RE"), FVec("staple1xy16_20IM")},
                                  {FVec("staple1xy16_21RE"), FVec("staple1xy16_21IM")},
                                  {FVec("staple1xy16_22RE"), FVec("staple1xy16_22IM")} } }, 
					{ { {FVec("staple1xy17_00RE"), FVec("staple1xy17_00IM")},
                                  {FVec("staple1xy17_01RE"), FVec("staple1xy17_01IM")},
                                  {FVec("staple1xy17_02RE"), FVec("staple1xy17_02IM")} },
                                { {FVec("staple1xy17_10RE"), FVec("staple1xy17_10IM")},
                                  {FVec("staple1xy17_11RE"), FVec("staple1xy17_11IM")},
                                  {FVec("staple1xy17_12RE"), FVec("staple1xy17_12IM")} },
                                { {FVec("staple1xy17_20RE"), FVec("staple1xy17_20IM")},
                                  {FVec("staple1xy17_21RE"), FVec("staple1xy17_21IM")},
                                  {FVec("staple1xy17_22RE"), FVec("staple1xy17_22IM")} } } },

					{ { { {FVec("staple1xz00_00RE"), FVec("staple1xz00_00IM")},
                                  {FVec("staple1xz00_01RE"), FVec("staple1xz00_01IM")},
                                  {FVec("staple1xz00_02RE"), FVec("staple1xz00_02IM")} },
                                { {FVec("staple1xz00_10RE"), FVec("staple1xz00_10IM")},
                                  {FVec("staple1xz00_11RE"), FVec("staple1xz00_11IM")},
                                  {FVec("staple1xz00_12RE"), FVec("staple1xz00_12IM")} },
                                { {FVec("staple1xz00_20RE"), FVec("staple1xz00_20IM")},
                                  {FVec("staple1xz00_21RE"), FVec("staple1xz00_21IM")},
                                  {FVec("staple1xz00_22RE"), FVec("staple1xz00_22IM")} } },
                                        { { {FVec("staple1xz01_00RE"), FVec("staple1xz01_00IM")},
                                  {FVec("staple1xz01_01RE"), FVec("staple1xz01_01IM")},
                                  {FVec("staple1xz01_02RE"), FVec("staple1xz01_02IM")} },
                                { {FVec("staple1xz01_10RE"), FVec("staple1xz01_10IM")},
                                  {FVec("staple1xz01_11RE"), FVec("staple1xz01_11IM")},
                                  {FVec("staple1xz01_12RE"), FVec("staple1xz01_12IM")} },
                                { {FVec("staple1xz01_20RE"), FVec("staple1xz01_20IM")},
                                  {FVec("staple1xz01_21RE"), FVec("staple1xz01_21IM")},
                                  {FVec("staple1xz01_22RE"), FVec("staple1xz01_22IM")} } },
                                        { { {FVec("staple1xz02_00RE"), FVec("staple1xz02_00IM")},
                                  {FVec("staple1xz02_01RE"), FVec("staple1xz02_01IM")},
                                  {FVec("staple1xz02_02RE"), FVec("staple1xz02_02IM")} },
                                { {FVec("staple1xz02_10RE"), FVec("staple1xz02_10IM")},
                                  {FVec("staple1xz02_11RE"), FVec("staple1xz02_11IM")},
                                  {FVec("staple1xz02_12RE"), FVec("staple1xz02_12IM")} },
                                { {FVec("staple1xz02_20RE"), FVec("staple1xz02_20IM")},
                                  {FVec("staple1xz02_21RE"), FVec("staple1xz02_21IM")},
                                  {FVec("staple1xz02_22RE"), FVec("staple1xz02_22IM")} } },
                                        { { {FVec("staple1xz03_00RE"), FVec("staple1xz03_00IM")},
                                  {FVec("staple1xz03_01RE"), FVec("staple1xz03_01IM")},
                                  {FVec("staple1xz03_02RE"), FVec("staple1xz03_02IM")} },
                                { {FVec("staple1xz03_10RE"), FVec("staple1xz03_10IM")},
                                  {FVec("staple1xz03_11RE"), FVec("staple1xz03_11IM")},
                                  {FVec("staple1xz03_12RE"), FVec("staple1xz03_12IM")} },
                                { {FVec("staple1xz03_20RE"), FVec("staple1xz03_20IM")},
                                  {FVec("staple1xz03_21RE"), FVec("staple1xz03_21IM")},
                                  {FVec("staple1xz03_22RE"), FVec("staple1xz03_22IM")} } },
                                        { { {FVec("staple1xz04_00RE"), FVec("staple1xz04_00IM")},
                                  {FVec("staple1xz04_01RE"), FVec("staple1xz04_01IM")},
                                  {FVec("staple1xz04_02RE"), FVec("staple1xz04_02IM")} },
                                { {FVec("staple1xz04_10RE"), FVec("staple1xz04_10IM")},
                                  {FVec("staple1xz04_11RE"), FVec("staple1xz04_11IM")},
                                  {FVec("staple1xz04_12RE"), FVec("staple1xz04_12IM")} },
                                { {FVec("staple1xz04_20RE"), FVec("staple1xz04_20IM")},
                                  {FVec("staple1xz04_21RE"), FVec("staple1xz04_21IM")},
                                  {FVec("staple1xz04_22RE"), FVec("staple1xz04_22IM")} } },
                                        { { {FVec("staple1xz05_00RE"), FVec("staple1xz05_00IM")},
                                  {FVec("staple1xz05_01RE"), FVec("staple1xz05_01IM")},
                                  {FVec("staple1xz05_02RE"), FVec("staple1xz05_02IM")} },
                                { {FVec("staple1xz05_10RE"), FVec("staple1xz05_10IM")},
                                  {FVec("staple1xz05_11RE"), FVec("staple1xz05_11IM")},
                                  {FVec("staple1xz05_12RE"), FVec("staple1xz05_12IM")} },
                                { {FVec("staple1xz05_20RE"), FVec("staple1xz05_20IM")},
                                  {FVec("staple1xz05_21RE"), FVec("staple1xz05_21IM")},
                                  {FVec("staple1xz05_22RE"), FVec("staple1xz05_22IM")} } },
                                        { { {FVec("staple1xz06_00RE"), FVec("staple1xz06_00IM")},
                                  {FVec("staple1xz06_01RE"), FVec("staple1xz06_01IM")},
                                  {FVec("staple1xz06_02RE"), FVec("staple1xz06_02IM")} },
                                { {FVec("staple1xz06_10RE"), FVec("staple1xz06_10IM")},
                                  {FVec("staple1xz06_11RE"), FVec("staple1xz06_11IM")},
                                  {FVec("staple1xz06_12RE"), FVec("staple1xz06_12IM")} },
                                { {FVec("staple1xz06_20RE"), FVec("staple1xz06_20IM")},
                                  {FVec("staple1xz06_21RE"), FVec("staple1xz06_21IM")},
                                  {FVec("staple1xz06_22RE"), FVec("staple1xz06_22IM")} } },
                                        { { {FVec("staple1xz07_00RE"), FVec("staple1xz07_00IM")},
                                  {FVec("staple1xz07_01RE"), FVec("staple1xz07_01IM")},
                                  {FVec("staple1xz07_02RE"), FVec("staple1xz07_02IM")} },
                                { {FVec("staple1xz07_10RE"), FVec("staple1xz07_10IM")},
                                  {FVec("staple1xz07_11RE"), FVec("staple1xz07_11IM")},
                                  {FVec("staple1xz07_12RE"), FVec("staple1xz07_12IM")} },
                                { {FVec("staple1xz07_20RE"), FVec("staple1xz07_20IM")},
                                  {FVec("staple1xz07_21RE"), FVec("staple1xz07_21IM")},
                                  {FVec("staple1xz07_22RE"), FVec("staple1xz07_22IM")} } },
                                        { { {FVec("staple1xz10_00RE"), FVec("staple1xz10_00IM")},
                                  {FVec("staple1xz10_01RE"), FVec("staple1xz10_01IM")},
                                  {FVec("staple1xz10_02RE"), FVec("staple1xz10_02IM")} },
                                { {FVec("staple1xz10_10RE"), FVec("staple1xz10_10IM")},
                                  {FVec("staple1xz10_11RE"), FVec("staple1xz10_11IM")},
                                  {FVec("staple1xz10_12RE"), FVec("staple1xz10_12IM")} },
                                { {FVec("staple1xz10_20RE"), FVec("staple1xz10_20IM")},
                                  {FVec("staple1xz10_21RE"), FVec("staple1xz10_21IM")},
                                  {FVec("staple1xz10_22RE"), FVec("staple1xz10_22IM")} } },
                                        { { {FVec("staple1xz11_00RE"), FVec("staple1xz11_00IM")},
                                  {FVec("staple1xz11_01RE"), FVec("staple1xz11_01IM")},
                                  {FVec("staple1xz11_02RE"), FVec("staple1xz11_02IM")} },
                                { {FVec("staple1xz11_10RE"), FVec("staple1xz11_10IM")},
                                  {FVec("staple1xz11_11RE"), FVec("staple1xz11_11IM")},
                                  {FVec("staple1xz11_12RE"), FVec("staple1xz11_12IM")} },
                                { {FVec("staple1xz11_20RE"), FVec("staple1xz11_20IM")},
                                  {FVec("staple1xz11_21RE"), FVec("staple1xz11_21IM")},
                                  {FVec("staple1xz11_22RE"), FVec("staple1xz11_22IM")} } },
                                        { { {FVec("staple1xz12_00RE"), FVec("staple1xz12_00IM")},
                                  {FVec("staple1xz12_01RE"), FVec("staple1xz12_01IM")},
                                  {FVec("staple1xz12_02RE"), FVec("staple1xz12_02IM")} },
                                { {FVec("staple1xz12_10RE"), FVec("staple1xz12_10IM")},
                                  {FVec("staple1xz12_11RE"), FVec("staple1xz12_11IM")},
                                  {FVec("staple1xz12_12RE"), FVec("staple1xz12_12IM")} },
                                { {FVec("staple1xz12_20RE"), FVec("staple1xz12_20IM")},
                                  {FVec("staple1xz12_21RE"), FVec("staple1xz12_21IM")},
                                  {FVec("staple1xz12_22RE"), FVec("staple1xz12_22IM")} } },
                                        { { {FVec("staple1xz13_00RE"), FVec("staple1xz13_00IM")},
                                  {FVec("staple1xz13_01RE"), FVec("staple1xz13_01IM")},
                                  {FVec("staple1xz13_02RE"), FVec("staple1xz13_02IM")} },
                                { {FVec("staple1xz13_10RE"), FVec("staple1xz13_10IM")},
                                  {FVec("staple1xz13_11RE"), FVec("staple1xz13_11IM")},
                                  {FVec("staple1xz13_12RE"), FVec("staple1xz13_12IM")} },
                                { {FVec("staple1xz13_20RE"), FVec("staple1xz13_20IM")},
                                  {FVec("staple1xz13_21RE"), FVec("staple1xz13_21IM")},
                                  {FVec("staple1xz13_22RE"), FVec("staple1xz13_22IM")} } },
                                        { { {FVec("staple1xz14_00RE"), FVec("staple1xz14_00IM")},
                                  {FVec("staple1xz14_01RE"), FVec("staple1xz14_01IM")},
                                  {FVec("staple1xz14_02RE"), FVec("staple1xz14_02IM")} },
                                { {FVec("staple1xz14_10RE"), FVec("staple1xz14_10IM")},
                                  {FVec("staple1xz14_11RE"), FVec("staple1xz14_11IM")},
                                  {FVec("staple1xz14_12RE"), FVec("staple1xz14_12IM")} },
                                { {FVec("staple1xz14_20RE"), FVec("staple1xz14_20IM")},
                                  {FVec("staple1xz14_21RE"), FVec("staple1xz14_21IM")},
                                  {FVec("staple1xz14_22RE"), FVec("staple1xz14_22IM")} } },
                                        { { {FVec("staple1xz15_00RE"), FVec("staple1xz15_00IM")},
                                  {FVec("staple1xz15_01RE"), FVec("staple1xz15_01IM")},
                                  {FVec("staple1xz15_02RE"), FVec("staple1xz15_02IM")} },
                                { {FVec("staple1xz15_10RE"), FVec("staple1xz15_10IM")},
                                  {FVec("staple1xz15_11RE"), FVec("staple1xz15_11IM")},
                                  {FVec("staple1xz15_12RE"), FVec("staple1xz15_12IM")} },
                                { {FVec("staple1xz15_20RE"), FVec("staple1xz15_20IM")},
                                  {FVec("staple1xz15_21RE"), FVec("staple1xz15_21IM")},
                                  {FVec("staple1xz15_22RE"), FVec("staple1xz15_22IM")} } },
                                        { { {FVec("staple1xz16_00RE"), FVec("staple1xz16_00IM")},
                                  {FVec("staple1xz16_01RE"), FVec("staple1xz16_01IM")},
                                  {FVec("staple1xz16_02RE"), FVec("staple1xz16_02IM")} },
                                { {FVec("staple1xz16_10RE"), FVec("staple1xz16_10IM")},
                                  {FVec("staple1xz16_11RE"), FVec("staple1xz16_11IM")},
                                  {FVec("staple1xz16_12RE"), FVec("staple1xz16_12IM")} },
                                { {FVec("staple1xz16_20RE"), FVec("staple1xz16_20IM")},
                                  {FVec("staple1xz16_21RE"), FVec("staple1xz16_21IM")},
                                  {FVec("staple1xz16_22RE"), FVec("staple1xz16_22IM")} } },
                                        { { {FVec("staple1xz17_00RE"), FVec("staple1xz17_00IM")},
                                  {FVec("staple1xz17_01RE"), FVec("staple1xz17_01IM")},
                                  {FVec("staple1xz17_02RE"), FVec("staple1xz17_02IM")} },
                                { {FVec("staple1xz17_10RE"), FVec("staple1xz17_10IM")},
                                  {FVec("staple1xz17_11RE"), FVec("staple1xz17_11IM")},
                                  {FVec("staple1xz17_12RE"), FVec("staple1xz17_12IM")} },
                                { {FVec("staple1xz17_20RE"), FVec("staple1xz17_20IM")},
                                  {FVec("staple1xz17_21RE"), FVec("staple1xz17_21IM")},
                                  {FVec("staple1xz17_22RE"), FVec("staple1xz17_22IM")} } } },

					{ { { {FVec("staple1xt00_00RE"), FVec("staple1xt00_00IM")},
                                  {FVec("staple1xt00_01RE"), FVec("staple1xt00_01IM")},
                                  {FVec("staple1xt00_02RE"), FVec("staple1xt00_02IM")} },
                                { {FVec("staple1xt00_10RE"), FVec("staple1xt00_10IM")},
                                  {FVec("staple1xt00_11RE"), FVec("staple1xt00_11IM")},
                                  {FVec("staple1xt00_12RE"), FVec("staple1xt00_12IM")} },
                                { {FVec("staple1xt00_20RE"), FVec("staple1xt00_20IM")},
                                  {FVec("staple1xt00_21RE"), FVec("staple1xt00_21IM")},
                                  {FVec("staple1xt00_22RE"), FVec("staple1xt00_22IM")} } },
                                        { { {FVec("staple1xt01_00RE"), FVec("staple1xt01_00IM")},
                                  {FVec("staple1xt01_01RE"), FVec("staple1xt01_01IM")},
                                  {FVec("staple1xt01_02RE"), FVec("staple1xt01_02IM")} },
                                { {FVec("staple1xt01_10RE"), FVec("staple1xt01_10IM")},
                                  {FVec("staple1xt01_11RE"), FVec("staple1xt01_11IM")},
                                  {FVec("staple1xt01_12RE"), FVec("staple1xt01_12IM")} },
                                { {FVec("staple1xt01_20RE"), FVec("staple1xt01_20IM")},
                                  {FVec("staple1xt01_21RE"), FVec("staple1xt01_21IM")},
                                  {FVec("staple1xt01_22RE"), FVec("staple1xt01_22IM")} } },
                                        { { {FVec("staple1xt02_00RE"), FVec("staple1xt02_00IM")},
                                  {FVec("staple1xt02_01RE"), FVec("staple1xt02_01IM")},
                                  {FVec("staple1xt02_02RE"), FVec("staple1xt02_02IM")} },
                                { {FVec("staple1xt02_10RE"), FVec("staple1xt02_10IM")},
                                  {FVec("staple1xt02_11RE"), FVec("staple1xt02_11IM")},
                                  {FVec("staple1xt02_12RE"), FVec("staple1xt02_12IM")} },
                                { {FVec("staple1xt02_20RE"), FVec("staple1xt02_20IM")},
                                  {FVec("staple1xt02_21RE"), FVec("staple1xt02_21IM")},
                                  {FVec("staple1xt02_22RE"), FVec("staple1xt02_22IM")} } },
                                        { { {FVec("staple1xt03_00RE"), FVec("staple1xt03_00IM")},
                                  {FVec("staple1xt03_01RE"), FVec("staple1xt03_01IM")},
                                  {FVec("staple1xt03_02RE"), FVec("staple1xt03_02IM")} },
                                { {FVec("staple1xt03_10RE"), FVec("staple1xt03_10IM")},
                                  {FVec("staple1xt03_11RE"), FVec("staple1xt03_11IM")},
                                  {FVec("staple1xt03_12RE"), FVec("staple1xt03_12IM")} },
                                { {FVec("staple1xt03_20RE"), FVec("staple1xt03_20IM")},
                                  {FVec("staple1xt03_21RE"), FVec("staple1xt03_21IM")},
                                  {FVec("staple1xt03_22RE"), FVec("staple1xt03_22IM")} } },
                                        { { {FVec("staple1xt04_00RE"), FVec("staple1xt04_00IM")},
                                  {FVec("staple1xt04_01RE"), FVec("staple1xt04_01IM")},
                                  {FVec("staple1xt04_02RE"), FVec("staple1xt04_02IM")} },
                                { {FVec("staple1xt04_10RE"), FVec("staple1xt04_10IM")},
                                  {FVec("staple1xt04_11RE"), FVec("staple1xt04_11IM")},
                                  {FVec("staple1xt04_12RE"), FVec("staple1xt04_12IM")} },
                                { {FVec("staple1xt04_20RE"), FVec("staple1xt04_20IM")},
                                  {FVec("staple1xt04_21RE"), FVec("staple1xt04_21IM")},
                                  {FVec("staple1xt04_22RE"), FVec("staple1xt04_22IM")} } },
                                        { { {FVec("staple1xt05_00RE"), FVec("staple1xt05_00IM")},
                                  {FVec("staple1xt05_01RE"), FVec("staple1xt05_01IM")},
                                  {FVec("staple1xt05_02RE"), FVec("staple1xt05_02IM")} },
                                { {FVec("staple1xt05_10RE"), FVec("staple1xt05_10IM")},
                                  {FVec("staple1xt05_11RE"), FVec("staple1xt05_11IM")},
                                  {FVec("staple1xt05_12RE"), FVec("staple1xt05_12IM")} },
                                { {FVec("staple1xt05_20RE"), FVec("staple1xt05_20IM")},
                                  {FVec("staple1xt05_21RE"), FVec("staple1xt05_21IM")},
                                  {FVec("staple1xt05_22RE"), FVec("staple1xt05_22IM")} } },
                                        { { {FVec("staple1xt06_00RE"), FVec("staple1xt06_00IM")},
                                  {FVec("staple1xt06_01RE"), FVec("staple1xt06_01IM")},
                                  {FVec("staple1xt06_02RE"), FVec("staple1xt06_02IM")} },
                                { {FVec("staple1xt06_10RE"), FVec("staple1xt06_10IM")},
                                  {FVec("staple1xt06_11RE"), FVec("staple1xt06_11IM")},
                                  {FVec("staple1xt06_12RE"), FVec("staple1xt06_12IM")} },
                                { {FVec("staple1xt06_20RE"), FVec("staple1xt06_20IM")},
                                  {FVec("staple1xt06_21RE"), FVec("staple1xt06_21IM")},
                                  {FVec("staple1xt06_22RE"), FVec("staple1xt06_22IM")} } },
                                        { { {FVec("staple1xt07_00RE"), FVec("staple1xt07_00IM")},
                                  {FVec("staple1xt07_01RE"), FVec("staple1xt07_01IM")},
                                  {FVec("staple1xt07_02RE"), FVec("staple1xt07_02IM")} },
                                { {FVec("staple1xt07_10RE"), FVec("staple1xt07_10IM")},
                                  {FVec("staple1xt07_11RE"), FVec("staple1xt07_11IM")},
                                  {FVec("staple1xt07_12RE"), FVec("staple1xt07_12IM")} },
                                { {FVec("staple1xt07_20RE"), FVec("staple1xt07_20IM")},
                                  {FVec("staple1xt07_21RE"), FVec("staple1xt07_21IM")},
                                  {FVec("staple1xt07_22RE"), FVec("staple1xt07_22IM")} } },
                                        { { {FVec("staple1xt10_00RE"), FVec("staple1xt10_00IM")},
                                  {FVec("staple1xt10_01RE"), FVec("staple1xt10_01IM")},
                                  {FVec("staple1xt10_02RE"), FVec("staple1xt10_02IM")} },
                                { {FVec("staple1xt10_10RE"), FVec("staple1xt10_10IM")},
                                  {FVec("staple1xt10_11RE"), FVec("staple1xt10_11IM")},
                                  {FVec("staple1xt10_12RE"), FVec("staple1xt10_12IM")} },
                                { {FVec("staple1xt10_20RE"), FVec("staple1xt10_20IM")},
                                  {FVec("staple1xt10_21RE"), FVec("staple1xt10_21IM")},
                                  {FVec("staple1xt10_22RE"), FVec("staple1xt10_22IM")} } },
                                        { { {FVec("staple1xt11_00RE"), FVec("staple1xt11_00IM")},
                                  {FVec("staple1xt11_01RE"), FVec("staple1xt11_01IM")},
                                  {FVec("staple1xt11_02RE"), FVec("staple1xt11_02IM")} },
                                { {FVec("staple1xt11_10RE"), FVec("staple1xt11_10IM")},
                                  {FVec("staple1xt11_11RE"), FVec("staple1xt11_11IM")},
                                  {FVec("staple1xt11_12RE"), FVec("staple1xt11_12IM")} },
                                { {FVec("staple1xt11_20RE"), FVec("staple1xt11_20IM")},
                                  {FVec("staple1xt11_21RE"), FVec("staple1xt11_21IM")},
                                  {FVec("staple1xt11_22RE"), FVec("staple1xt11_22IM")} } },
                                        { { {FVec("staple1xt12_00RE"), FVec("staple1xt12_00IM")},
                                  {FVec("staple1xt12_01RE"), FVec("staple1xt12_01IM")},
                                  {FVec("staple1xt12_02RE"), FVec("staple1xt12_02IM")} },
                                { {FVec("staple1xt12_10RE"), FVec("staple1xt12_10IM")},
                                  {FVec("staple1xt12_11RE"), FVec("staple1xt12_11IM")},
                                  {FVec("staple1xt12_12RE"), FVec("staple1xt12_12IM")} },
                                { {FVec("staple1xt12_20RE"), FVec("staple1xt12_20IM")},
                                  {FVec("staple1xt12_21RE"), FVec("staple1xt12_21IM")},
                                  {FVec("staple1xt12_22RE"), FVec("staple1xt12_22IM")} } },
                                        { { {FVec("staple1xt13_00RE"), FVec("staple1xt13_00IM")},
                                  {FVec("staple1xt13_01RE"), FVec("staple1xt13_01IM")},
                                  {FVec("staple1xt13_02RE"), FVec("staple1xt13_02IM")} },
                                { {FVec("staple1xt13_10RE"), FVec("staple1xt13_10IM")},
                                  {FVec("staple1xt13_11RE"), FVec("staple1xt13_11IM")},
                                  {FVec("staple1xt13_12RE"), FVec("staple1xt13_12IM")} },
                                { {FVec("staple1xt13_20RE"), FVec("staple1xt13_20IM")},
                                  {FVec("staple1xt13_21RE"), FVec("staple1xt13_21IM")},
                                  {FVec("staple1xt13_22RE"), FVec("staple1xt13_22IM")} } },
                                        { { {FVec("staple1xt14_00RE"), FVec("staple1xt14_00IM")},
                                  {FVec("staple1xt14_01RE"), FVec("staple1xt14_01IM")},
                                  {FVec("staple1xt14_02RE"), FVec("staple1xt14_02IM")} },
                                { {FVec("staple1xt14_10RE"), FVec("staple1xt14_10IM")},
                                  {FVec("staple1xt14_11RE"), FVec("staple1xt14_11IM")},
                                  {FVec("staple1xt14_12RE"), FVec("staple1xt14_12IM")} },
                                { {FVec("staple1xt14_20RE"), FVec("staple1xt14_20IM")},
                                  {FVec("staple1xt14_21RE"), FVec("staple1xt14_21IM")},
                                  {FVec("staple1xt14_22RE"), FVec("staple1xt14_22IM")} } },
                                        { { {FVec("staple1xt15_00RE"), FVec("staple1xt15_00IM")},
                                  {FVec("staple1xt15_01RE"), FVec("staple1xt15_01IM")},
                                  {FVec("staple1xt15_02RE"), FVec("staple1xt15_02IM")} },
                                { {FVec("staple1xt15_10RE"), FVec("staple1xt15_10IM")},
                                  {FVec("staple1xt15_11RE"), FVec("staple1xt15_11IM")},
                                  {FVec("staple1xt15_12RE"), FVec("staple1xt15_12IM")} },
                                { {FVec("staple1xt15_20RE"), FVec("staple1xt15_20IM")},
                                  {FVec("staple1xt15_21RE"), FVec("staple1xt15_21IM")},
                                  {FVec("staple1xt15_22RE"), FVec("staple1xt15_22IM")} } },
                                        { { {FVec("staple1xt16_00RE"), FVec("staple1xt16_00IM")},
                                  {FVec("staple1xt16_01RE"), FVec("staple1xt16_01IM")},
                                  {FVec("staple1xt16_02RE"), FVec("staple1xt16_02IM")} },
                                { {FVec("staple1xt16_10RE"), FVec("staple1xt16_10IM")},
                                  {FVec("staple1xt16_11RE"), FVec("staple1xt16_11IM")},
                                  {FVec("staple1xt16_12RE"), FVec("staple1xt16_12IM")} },
                                { {FVec("staple1xt16_20RE"), FVec("staple1xt16_20IM")},
                                  {FVec("staple1xt16_21RE"), FVec("staple1xt16_21IM")},
                                  {FVec("staple1xt16_22RE"), FVec("staple1xt16_22IM")} } },
                                        { { {FVec("staple1xt17_00RE"), FVec("staple1xt17_00IM")},
                                  {FVec("staple1xt17_01RE"), FVec("staple1xt17_01IM")},
                                  {FVec("staple1xt17_02RE"), FVec("staple1xt17_02IM")} },
                                { {FVec("staple1xt17_10RE"), FVec("staple1xt17_10IM")},
                                  {FVec("staple1xt17_11RE"), FVec("staple1xt17_11IM")},
                                  {FVec("staple1xt17_12RE"), FVec("staple1xt17_12IM")} },
                                { {FVec("staple1xt17_20RE"), FVec("staple1xt17_20IM")},
                                  {FVec("staple1xt17_21RE"), FVec("staple1xt17_21IM")},
                                  {FVec("staple1xt17_22RE"), FVec("staple1xt17_22IM")} } } }, 

					{ { { {FVec("staple1yz00_00RE"), FVec("staple1yz00_00IM")},
                                  {FVec("staple1yz00_01RE"), FVec("staple1yz00_01IM")},
                                  {FVec("staple1yz00_02RE"), FVec("staple1yz00_02IM")} },
                                { {FVec("staple1yz00_10RE"), FVec("staple1yz00_10IM")},
                                  {FVec("staple1yz00_11RE"), FVec("staple1yz00_11IM")},
                                  {FVec("staple1yz00_12RE"), FVec("staple1yz00_12IM")} },
                                { {FVec("staple1yz00_20RE"), FVec("staple1yz00_20IM")},
                                  {FVec("staple1yz00_21RE"), FVec("staple1yz00_21IM")},
                                  {FVec("staple1yz00_22RE"), FVec("staple1yz00_22IM")} } },
                                        { { {FVec("staple1yz01_00RE"), FVec("staple1yz01_00IM")},
                                  {FVec("staple1yz01_01RE"), FVec("staple1yz01_01IM")},
                                  {FVec("staple1yz01_02RE"), FVec("staple1yz01_02IM")} },
                                { {FVec("staple1yz01_10RE"), FVec("staple1yz01_10IM")},
                                  {FVec("staple1yz01_11RE"), FVec("staple1yz01_11IM")},
                                  {FVec("staple1yz01_12RE"), FVec("staple1yz01_12IM")} },
                                { {FVec("staple1yz01_20RE"), FVec("staple1yz01_20IM")},
                                  {FVec("staple1yz01_21RE"), FVec("staple1yz01_21IM")},
                                  {FVec("staple1yz01_22RE"), FVec("staple1yz01_22IM")} } },
                                        { { {FVec("staple1yz02_00RE"), FVec("staple1yz02_00IM")},
                                  {FVec("staple1yz02_01RE"), FVec("staple1yz02_01IM")},
                                  {FVec("staple1yz02_02RE"), FVec("staple1yz02_02IM")} },
                                { {FVec("staple1yz02_10RE"), FVec("staple1yz02_10IM")},
                                  {FVec("staple1yz02_11RE"), FVec("staple1yz02_11IM")},
                                  {FVec("staple1yz02_12RE"), FVec("staple1yz02_12IM")} },
                                { {FVec("staple1yz02_20RE"), FVec("staple1yz02_20IM")},
                                  {FVec("staple1yz02_21RE"), FVec("staple1yz02_21IM")},
                                  {FVec("staple1yz02_22RE"), FVec("staple1yz02_22IM")} } },
                                        { { {FVec("staple1yz03_00RE"), FVec("staple1yz03_00IM")},
                                  {FVec("staple1yz03_01RE"), FVec("staple1yz03_01IM")},
                                  {FVec("staple1yz03_02RE"), FVec("staple1yz03_02IM")} },
                                { {FVec("staple1yz03_10RE"), FVec("staple1yz03_10IM")},
                                  {FVec("staple1yz03_11RE"), FVec("staple1yz03_11IM")},
                                  {FVec("staple1yz03_12RE"), FVec("staple1yz03_12IM")} },
                                { {FVec("staple1yz03_20RE"), FVec("staple1yz03_20IM")},
                                  {FVec("staple1yz03_21RE"), FVec("staple1yz03_21IM")},
                                  {FVec("staple1yz03_22RE"), FVec("staple1yz03_22IM")} } },
                                        { { {FVec("staple1yz04_00RE"), FVec("staple1yz04_00IM")},
                                  {FVec("staple1yz04_01RE"), FVec("staple1yz04_01IM")},
                                  {FVec("staple1yz04_02RE"), FVec("staple1yz04_02IM")} },
                                { {FVec("staple1yz04_10RE"), FVec("staple1yz04_10IM")},
                                  {FVec("staple1yz04_11RE"), FVec("staple1yz04_11IM")},
                                  {FVec("staple1yz04_12RE"), FVec("staple1yz04_12IM")} },
                                { {FVec("staple1yz04_20RE"), FVec("staple1yz04_20IM")},
                                  {FVec("staple1yz04_21RE"), FVec("staple1yz04_21IM")},
                                  {FVec("staple1yz04_22RE"), FVec("staple1yz04_22IM")} } },
                                        { { {FVec("staple1yz05_00RE"), FVec("staple1yz05_00IM")},
                                  {FVec("staple1yz05_01RE"), FVec("staple1yz05_01IM")},
                                  {FVec("staple1yz05_02RE"), FVec("staple1yz05_02IM")} },
                                { {FVec("staple1yz05_10RE"), FVec("staple1yz05_10IM")},
                                  {FVec("staple1yz05_11RE"), FVec("staple1yz05_11IM")},
                                  {FVec("staple1yz05_12RE"), FVec("staple1yz05_12IM")} },
                                { {FVec("staple1yz05_20RE"), FVec("staple1yz05_20IM")},
                                  {FVec("staple1yz05_21RE"), FVec("staple1yz05_21IM")},
                                  {FVec("staple1yz05_22RE"), FVec("staple1yz05_22IM")} } },
                                        { { {FVec("staple1yz06_00RE"), FVec("staple1yz06_00IM")},
                                  {FVec("staple1yz06_01RE"), FVec("staple1yz06_01IM")},
                                  {FVec("staple1yz06_02RE"), FVec("staple1yz06_02IM")} },
                                { {FVec("staple1yz06_10RE"), FVec("staple1yz06_10IM")},
                                  {FVec("staple1yz06_11RE"), FVec("staple1yz06_11IM")},
                                  {FVec("staple1yz06_12RE"), FVec("staple1yz06_12IM")} },
                                { {FVec("staple1yz06_20RE"), FVec("staple1yz06_20IM")},
                                  {FVec("staple1yz06_21RE"), FVec("staple1yz06_21IM")},
                                  {FVec("staple1yz06_22RE"), FVec("staple1yz06_22IM")} } },
                                        { { {FVec("staple1yz07_00RE"), FVec("staple1yz07_00IM")},
                                  {FVec("staple1yz07_01RE"), FVec("staple1yz07_01IM")},
                                  {FVec("staple1yz07_02RE"), FVec("staple1yz07_02IM")} },
                                { {FVec("staple1yz07_10RE"), FVec("staple1yz07_10IM")},
                                  {FVec("staple1yz07_11RE"), FVec("staple1yz07_11IM")},
                                  {FVec("staple1yz07_12RE"), FVec("staple1yz07_12IM")} },
                                { {FVec("staple1yz07_20RE"), FVec("staple1yz07_20IM")},
                                  {FVec("staple1yz07_21RE"), FVec("staple1yz07_21IM")},
                                  {FVec("staple1yz07_22RE"), FVec("staple1yz07_22IM")} } },
                                        { { {FVec("staple1yz10_00RE"), FVec("staple1yz10_00IM")},
                                  {FVec("staple1yz10_01RE"), FVec("staple1yz10_01IM")},
                                  {FVec("staple1yz10_02RE"), FVec("staple1yz10_02IM")} },
                                { {FVec("staple1yz10_10RE"), FVec("staple1yz10_10IM")},
                                  {FVec("staple1yz10_11RE"), FVec("staple1yz10_11IM")},
                                  {FVec("staple1yz10_12RE"), FVec("staple1yz10_12IM")} },
                                { {FVec("staple1yz10_20RE"), FVec("staple1yz10_20IM")},
                                  {FVec("staple1yz10_21RE"), FVec("staple1yz10_21IM")},
                                  {FVec("staple1yz10_22RE"), FVec("staple1yz10_22IM")} } },
                                        { { {FVec("staple1yz11_00RE"), FVec("staple1yz11_00IM")},
                                  {FVec("staple1yz11_01RE"), FVec("staple1yz11_01IM")},
                                  {FVec("staple1yz11_02RE"), FVec("staple1yz11_02IM")} },
                                { {FVec("staple1yz11_10RE"), FVec("staple1yz11_10IM")},
                                  {FVec("staple1yz11_11RE"), FVec("staple1yz11_11IM")},
                                  {FVec("staple1yz11_12RE"), FVec("staple1yz11_12IM")} },
                                { {FVec("staple1yz11_20RE"), FVec("staple1yz11_20IM")},
                                  {FVec("staple1yz11_21RE"), FVec("staple1yz11_21IM")},
                                  {FVec("staple1yz11_22RE"), FVec("staple1yz11_22IM")} } },
                                        { { {FVec("staple1yz12_00RE"), FVec("staple1yz12_00IM")},
                                  {FVec("staple1yz12_01RE"), FVec("staple1yz12_01IM")},
                                  {FVec("staple1yz12_02RE"), FVec("staple1yz12_02IM")} },
                                { {FVec("staple1yz12_10RE"), FVec("staple1yz12_10IM")},
                                  {FVec("staple1yz12_11RE"), FVec("staple1yz12_11IM")},
                                  {FVec("staple1yz12_12RE"), FVec("staple1yz12_12IM")} },
                                { {FVec("staple1yz12_20RE"), FVec("staple1yz12_20IM")},
                                  {FVec("staple1yz12_21RE"), FVec("staple1yz12_21IM")},
                                  {FVec("staple1yz12_22RE"), FVec("staple1yz12_22IM")} } },
                                        { { {FVec("staple1yz13_00RE"), FVec("staple1yz13_00IM")},
                                  {FVec("staple1yz13_01RE"), FVec("staple1yz13_01IM")},
                                  {FVec("staple1yz13_02RE"), FVec("staple1yz13_02IM")} },
                                { {FVec("staple1yz13_10RE"), FVec("staple1yz13_10IM")},
                                  {FVec("staple1yz13_11RE"), FVec("staple1yz13_11IM")},
                                  {FVec("staple1yz13_12RE"), FVec("staple1yz13_12IM")} },
                                { {FVec("staple1yz13_20RE"), FVec("staple1yz13_20IM")},
                                  {FVec("staple1yz13_21RE"), FVec("staple1yz13_21IM")},
                                  {FVec("staple1yz13_22RE"), FVec("staple1yz13_22IM")} } },
                                        { { {FVec("staple1yz14_00RE"), FVec("staple1yz14_00IM")},
                                  {FVec("staple1yz14_01RE"), FVec("staple1yz14_01IM")},
                                  {FVec("staple1yz14_02RE"), FVec("staple1yz14_02IM")} },
                                { {FVec("staple1yz14_10RE"), FVec("staple1yz14_10IM")},
                                  {FVec("staple1yz14_11RE"), FVec("staple1yz14_11IM")},
                                  {FVec("staple1yz14_12RE"), FVec("staple1yz14_12IM")} },
                                { {FVec("staple1yz14_20RE"), FVec("staple1yz14_20IM")},
                                  {FVec("staple1yz14_21RE"), FVec("staple1yz14_21IM")},
                                  {FVec("staple1yz14_22RE"), FVec("staple1yz14_22IM")} } },
                                        { { {FVec("staple1yz15_00RE"), FVec("staple1yz15_00IM")},
                                  {FVec("staple1yz15_01RE"), FVec("staple1yz15_01IM")},
                                  {FVec("staple1yz15_02RE"), FVec("staple1yz15_02IM")} },
                                { {FVec("staple1yz15_10RE"), FVec("staple1yz15_10IM")},
                                  {FVec("staple1yz15_11RE"), FVec("staple1yz15_11IM")},
                                  {FVec("staple1yz15_12RE"), FVec("staple1yz15_12IM")} },
                                { {FVec("staple1yz15_20RE"), FVec("staple1yz15_20IM")},
                                  {FVec("staple1yz15_21RE"), FVec("staple1yz15_21IM")},
                                  {FVec("staple1yz15_22RE"), FVec("staple1yz15_22IM")} } },
                                        { { {FVec("staple1yz16_00RE"), FVec("staple1yz16_00IM")},
                                  {FVec("staple1yz16_01RE"), FVec("staple1yz16_01IM")},
                                  {FVec("staple1yz16_02RE"), FVec("staple1yz16_02IM")} },
                                { {FVec("staple1yz16_10RE"), FVec("staple1yz16_10IM")},
                                  {FVec("staple1yz16_11RE"), FVec("staple1yz16_11IM")},
                                  {FVec("staple1yz16_12RE"), FVec("staple1yz16_12IM")} },
                                { {FVec("staple1yz16_20RE"), FVec("staple1yz16_20IM")},
                                  {FVec("staple1yz16_21RE"), FVec("staple1yz16_21IM")},
                                  {FVec("staple1yz16_22RE"), FVec("staple1yz16_22IM")} } },
                                        { { {FVec("staple1yz17_00RE"), FVec("staple1yz17_00IM")},
                                  {FVec("staple1yz17_01RE"), FVec("staple1yz17_01IM")},
                                  {FVec("staple1yz17_02RE"), FVec("staple1yz17_02IM")} },
                                { {FVec("staple1yz17_10RE"), FVec("staple1yz17_10IM")},
                                  {FVec("staple1yz17_11RE"), FVec("staple1yz17_11IM")},
                                  {FVec("staple1yz17_12RE"), FVec("staple1yz17_12IM")} },
                                { {FVec("staple1yz17_20RE"), FVec("staple1yz17_20IM")},
                                  {FVec("staple1yz17_21RE"), FVec("staple1yz17_21IM")},
                                  {FVec("staple1yz17_22RE"), FVec("staple1yz17_22IM")} } } }, 

					{ { { {FVec("staple1yt00_00RE"), FVec("staple1yt00_00IM")},
                                  {FVec("staple1yt00_01RE"), FVec("staple1yt00_01IM")},
                                  {FVec("staple1yt00_02RE"), FVec("staple1yt00_02IM")} },
                                { {FVec("staple1yt00_10RE"), FVec("staple1yt00_10IM")},
                                  {FVec("staple1yt00_11RE"), FVec("staple1yt00_11IM")},
                                  {FVec("staple1yt00_12RE"), FVec("staple1yt00_12IM")} },
                                { {FVec("staple1yt00_20RE"), FVec("staple1yt00_20IM")},
                                  {FVec("staple1yt00_21RE"), FVec("staple1yt00_21IM")},
                                  {FVec("staple1yt00_22RE"), FVec("staple1yt00_22IM")} } },
                                        { { {FVec("staple1yt01_00RE"), FVec("staple1yt01_00IM")},
                                  {FVec("staple1yt01_01RE"), FVec("staple1yt01_01IM")},
                                  {FVec("staple1yt01_02RE"), FVec("staple1yt01_02IM")} },
                                { {FVec("staple1yt01_10RE"), FVec("staple1yt01_10IM")},
                                  {FVec("staple1yt01_11RE"), FVec("staple1yt01_11IM")},
                                  {FVec("staple1yt01_12RE"), FVec("staple1yt01_12IM")} },
                                { {FVec("staple1yt01_20RE"), FVec("staple1yt01_20IM")},
                                  {FVec("staple1yt01_21RE"), FVec("staple1yt01_21IM")},
                                  {FVec("staple1yt01_22RE"), FVec("staple1yt01_22IM")} } },
                                        { { {FVec("staple1yt02_00RE"), FVec("staple1yt02_00IM")},
                                  {FVec("staple1yt02_01RE"), FVec("staple1yt02_01IM")},
                                  {FVec("staple1yt02_02RE"), FVec("staple1yt02_02IM")} },
                                { {FVec("staple1yt02_10RE"), FVec("staple1yt02_10IM")},
                                  {FVec("staple1yt02_11RE"), FVec("staple1yt02_11IM")},
                                  {FVec("staple1yt02_12RE"), FVec("staple1yt02_12IM")} },
                                { {FVec("staple1yt02_20RE"), FVec("staple1yt02_20IM")},
                                  {FVec("staple1yt02_21RE"), FVec("staple1yt02_21IM")},
                                  {FVec("staple1yt02_22RE"), FVec("staple1yt02_22IM")} } },
                                        { { {FVec("staple1yt03_00RE"), FVec("staple1yt03_00IM")},
                                  {FVec("staple1yt03_01RE"), FVec("staple1yt03_01IM")},
                                  {FVec("staple1yt03_02RE"), FVec("staple1yt03_02IM")} },
                                { {FVec("staple1yt03_10RE"), FVec("staple1yt03_10IM")},
                                  {FVec("staple1yt03_11RE"), FVec("staple1yt03_11IM")},
                                  {FVec("staple1yt03_12RE"), FVec("staple1yt03_12IM")} },
                                { {FVec("staple1yt03_20RE"), FVec("staple1yt03_20IM")},
                                  {FVec("staple1yt03_21RE"), FVec("staple1yt03_21IM")},
                                  {FVec("staple1yt03_22RE"), FVec("staple1yt03_22IM")} } },
                                        { { {FVec("staple1yt04_00RE"), FVec("staple1yt04_00IM")},
                                  {FVec("staple1yt04_01RE"), FVec("staple1yt04_01IM")},
                                  {FVec("staple1yt04_02RE"), FVec("staple1yt04_02IM")} },
                                { {FVec("staple1yt04_10RE"), FVec("staple1yt04_10IM")},
                                  {FVec("staple1yt04_11RE"), FVec("staple1yt04_11IM")},
                                  {FVec("staple1yt04_12RE"), FVec("staple1yt04_12IM")} },
                                { {FVec("staple1yt04_20RE"), FVec("staple1yt04_20IM")},
                                  {FVec("staple1yt04_21RE"), FVec("staple1yt04_21IM")},
                                  {FVec("staple1yt04_22RE"), FVec("staple1yt04_22IM")} } },
                                        { { {FVec("staple1yt05_00RE"), FVec("staple1yt05_00IM")},
                                  {FVec("staple1yt05_01RE"), FVec("staple1yt05_01IM")},
                                  {FVec("staple1yt05_02RE"), FVec("staple1yt05_02IM")} },
                                { {FVec("staple1yt05_10RE"), FVec("staple1yt05_10IM")},
                                  {FVec("staple1yt05_11RE"), FVec("staple1yt05_11IM")},
                                  {FVec("staple1yt05_12RE"), FVec("staple1yt05_12IM")} },
                                { {FVec("staple1yt05_20RE"), FVec("staple1yt05_20IM")},
                                  {FVec("staple1yt05_21RE"), FVec("staple1yt05_21IM")},
                                  {FVec("staple1yt05_22RE"), FVec("staple1yt05_22IM")} } },
                                        { { {FVec("staple1yt06_00RE"), FVec("staple1yt06_00IM")},
                                  {FVec("staple1yt06_01RE"), FVec("staple1yt06_01IM")},
                                  {FVec("staple1yt06_02RE"), FVec("staple1yt06_02IM")} },
                                { {FVec("staple1yt06_10RE"), FVec("staple1yt06_10IM")},
                                  {FVec("staple1yt06_11RE"), FVec("staple1yt06_11IM")},
                                  {FVec("staple1yt06_12RE"), FVec("staple1yt06_12IM")} },
                                { {FVec("staple1yt06_20RE"), FVec("staple1yt06_20IM")},
                                  {FVec("staple1yt06_21RE"), FVec("staple1yt06_21IM")},
                                  {FVec("staple1yt06_22RE"), FVec("staple1yt06_22IM")} } },
                                        { { {FVec("staple1yt07_00RE"), FVec("staple1yt07_00IM")},
                                  {FVec("staple1yt07_01RE"), FVec("staple1yt07_01IM")},
                                  {FVec("staple1yt07_02RE"), FVec("staple1yt07_02IM")} },
                                { {FVec("staple1yt07_10RE"), FVec("staple1yt07_10IM")},
                                  {FVec("staple1yt07_11RE"), FVec("staple1yt07_11IM")},
                                  {FVec("staple1yt07_12RE"), FVec("staple1yt07_12IM")} },
                                { {FVec("staple1yt07_20RE"), FVec("staple1yt07_20IM")},
                                  {FVec("staple1yt07_21RE"), FVec("staple1yt07_21IM")},
                                  {FVec("staple1yt07_22RE"), FVec("staple1yt07_22IM")} } },
                                        { { {FVec("staple1yt10_00RE"), FVec("staple1yt10_00IM")},
                                  {FVec("staple1yt10_01RE"), FVec("staple1yt10_01IM")},
                                  {FVec("staple1yt10_02RE"), FVec("staple1yt10_02IM")} },
                                { {FVec("staple1yt10_10RE"), FVec("staple1yt10_10IM")},
                                  {FVec("staple1yt10_11RE"), FVec("staple1yt10_11IM")},
                                  {FVec("staple1yt10_12RE"), FVec("staple1yt10_12IM")} },
                                { {FVec("staple1yt10_20RE"), FVec("staple1yt10_20IM")},
                                  {FVec("staple1yt10_21RE"), FVec("staple1yt10_21IM")},
                                  {FVec("staple1yt10_22RE"), FVec("staple1yt10_22IM")} } },
                                        { { {FVec("staple1yt11_00RE"), FVec("staple1yt11_00IM")},
                                  {FVec("staple1yt11_01RE"), FVec("staple1yt11_01IM")},
                                  {FVec("staple1yt11_02RE"), FVec("staple1yt11_02IM")} },
                                { {FVec("staple1yt11_10RE"), FVec("staple1yt11_10IM")},
                                  {FVec("staple1yt11_11RE"), FVec("staple1yt11_11IM")},
                                  {FVec("staple1yt11_12RE"), FVec("staple1yt11_12IM")} },
                                { {FVec("staple1yt11_20RE"), FVec("staple1yt11_20IM")},
                                  {FVec("staple1yt11_21RE"), FVec("staple1yt11_21IM")},
                                  {FVec("staple1yt11_22RE"), FVec("staple1yt11_22IM")} } },
                                        { { {FVec("staple1yt12_00RE"), FVec("staple1yt12_00IM")},
                                  {FVec("staple1yt12_01RE"), FVec("staple1yt12_01IM")},
                                  {FVec("staple1yt12_02RE"), FVec("staple1yt12_02IM")} },
                                { {FVec("staple1yt12_10RE"), FVec("staple1yt12_10IM")},
                                  {FVec("staple1yt12_11RE"), FVec("staple1yt12_11IM")},
                                  {FVec("staple1yt12_12RE"), FVec("staple1yt12_12IM")} },
                                { {FVec("staple1yt12_20RE"), FVec("staple1yt12_20IM")},
                                  {FVec("staple1yt12_21RE"), FVec("staple1yt12_21IM")},
                                  {FVec("staple1yt12_22RE"), FVec("staple1yt12_22IM")} } },
                                        { { {FVec("staple1yt13_00RE"), FVec("staple1yt13_00IM")},
                                  {FVec("staple1yt13_01RE"), FVec("staple1yt13_01IM")},
                                  {FVec("staple1yt13_02RE"), FVec("staple1yt13_02IM")} },
                                { {FVec("staple1yt13_10RE"), FVec("staple1yt13_10IM")},
                                  {FVec("staple1yt13_11RE"), FVec("staple1yt13_11IM")},
                                  {FVec("staple1yt13_12RE"), FVec("staple1yt13_12IM")} },
                                { {FVec("staple1yt13_20RE"), FVec("staple1yt13_20IM")},
                                  {FVec("staple1yt13_21RE"), FVec("staple1yt13_21IM")},
                                  {FVec("staple1yt13_22RE"), FVec("staple1yt13_22IM")} } },
                                        { { {FVec("staple1yt14_00RE"), FVec("staple1yt14_00IM")},
                                  {FVec("staple1yt14_01RE"), FVec("staple1yt14_01IM")},
                                  {FVec("staple1yt14_02RE"), FVec("staple1yt14_02IM")} },
                                { {FVec("staple1yt14_10RE"), FVec("staple1yt14_10IM")},
                                  {FVec("staple1yt14_11RE"), FVec("staple1yt14_11IM")},
                                  {FVec("staple1yt14_12RE"), FVec("staple1yt14_12IM")} },
                                { {FVec("staple1yt14_20RE"), FVec("staple1yt14_20IM")},
                                  {FVec("staple1yt14_21RE"), FVec("staple1yt14_21IM")},
                                  {FVec("staple1yt14_22RE"), FVec("staple1yt14_22IM")} } },
                                        { { {FVec("staple1yt15_00RE"), FVec("staple1yt15_00IM")},
                                  {FVec("staple1yt15_01RE"), FVec("staple1yt15_01IM")},
                                  {FVec("staple1yt15_02RE"), FVec("staple1yt15_02IM")} },
                                { {FVec("staple1yt15_10RE"), FVec("staple1yt15_10IM")},
                                  {FVec("staple1yt15_11RE"), FVec("staple1yt15_11IM")},
                                  {FVec("staple1yt15_12RE"), FVec("staple1yt15_12IM")} },
                                { {FVec("staple1yt15_20RE"), FVec("staple1yt15_20IM")},
                                  {FVec("staple1yt15_21RE"), FVec("staple1yt15_21IM")},
                                  {FVec("staple1yt15_22RE"), FVec("staple1yt15_22IM")} } },
                                        { { {FVec("staple1yt16_00RE"), FVec("staple1yt16_00IM")},
                                  {FVec("staple1yt16_01RE"), FVec("staple1yt16_01IM")},
                                  {FVec("staple1yt16_02RE"), FVec("staple1yt16_02IM")} },
                                { {FVec("staple1yt16_10RE"), FVec("staple1yt16_10IM")},
                                  {FVec("staple1yt16_11RE"), FVec("staple1yt16_11IM")},
                                  {FVec("staple1yt16_12RE"), FVec("staple1yt16_12IM")} },
                                { {FVec("staple1yt16_20RE"), FVec("staple1yt16_20IM")},
                                  {FVec("staple1yt16_21RE"), FVec("staple1yt16_21IM")},
                                  {FVec("staple1yt16_22RE"), FVec("staple1yt16_22IM")} } },
                                        { { {FVec("staple1yt17_00RE"), FVec("staple1yt17_00IM")},
                                  {FVec("staple1yt17_01RE"), FVec("staple1yt17_01IM")},
                                  {FVec("staple1yt17_02RE"), FVec("staple1yt17_02IM")} },
                                { {FVec("staple1yt17_10RE"), FVec("staple1yt17_10IM")},
                                  {FVec("staple1yt17_11RE"), FVec("staple1yt17_11IM")},
                                  {FVec("staple1yt17_12RE"), FVec("staple1yt17_12IM")} },
                                { {FVec("staple1yt17_20RE"), FVec("staple1yt17_20IM")},
                                  {FVec("staple1yt17_21RE"), FVec("staple1yt17_21IM")},
                                  {FVec("staple1yt17_22RE"), FVec("staple1yt17_22IM")} } } }, 

					{ { { {FVec("staple1zt00_00RE"), FVec("staple1zt00_00IM")},
                                  {FVec("staple1zt00_01RE"), FVec("staple1zt00_01IM")},
                                  {FVec("staple1zt00_02RE"), FVec("staple1zt00_02IM")} },
                                { {FVec("staple1zt00_10RE"), FVec("staple1zt00_10IM")},
                                  {FVec("staple1zt00_11RE"), FVec("staple1zt00_11IM")},
                                  {FVec("staple1zt00_12RE"), FVec("staple1zt00_12IM")} },
                                { {FVec("staple1zt00_20RE"), FVec("staple1zt00_20IM")},
                                  {FVec("staple1zt00_21RE"), FVec("staple1zt00_21IM")},
                                  {FVec("staple1zt00_22RE"), FVec("staple1zt00_22IM")} } },
                                        { { {FVec("staple1zt01_00RE"), FVec("staple1zt01_00IM")},
                                  {FVec("staple1zt01_01RE"), FVec("staple1zt01_01IM")},
                                  {FVec("staple1zt01_02RE"), FVec("staple1zt01_02IM")} },
                                { {FVec("staple1zt01_10RE"), FVec("staple1zt01_10IM")},
                                  {FVec("staple1zt01_11RE"), FVec("staple1zt01_11IM")},
                                  {FVec("staple1zt01_12RE"), FVec("staple1zt01_12IM")} },
                                { {FVec("staple1zt01_20RE"), FVec("staple1zt01_20IM")},
                                  {FVec("staple1zt01_21RE"), FVec("staple1zt01_21IM")},
                                  {FVec("staple1zt01_22RE"), FVec("staple1zt01_22IM")} } },
                                        { { {FVec("staple1zt02_00RE"), FVec("staple1zt02_00IM")},
                                  {FVec("staple1zt02_01RE"), FVec("staple1zt02_01IM")},
                                  {FVec("staple1zt02_02RE"), FVec("staple1zt02_02IM")} },
                                { {FVec("staple1zt02_10RE"), FVec("staple1zt02_10IM")},
                                  {FVec("staple1zt02_11RE"), FVec("staple1zt02_11IM")},
                                  {FVec("staple1zt02_12RE"), FVec("staple1zt02_12IM")} },
                                { {FVec("staple1zt02_20RE"), FVec("staple1zt02_20IM")},
                                  {FVec("staple1zt02_21RE"), FVec("staple1zt02_21IM")},
                                  {FVec("staple1zt02_22RE"), FVec("staple1zt02_22IM")} } },
                                        { { {FVec("staple1zt03_00RE"), FVec("staple1zt03_00IM")},
                                  {FVec("staple1zt03_01RE"), FVec("staple1zt03_01IM")},
                                  {FVec("staple1zt03_02RE"), FVec("staple1zt03_02IM")} },
                                { {FVec("staple1zt03_10RE"), FVec("staple1zt03_10IM")},
                                  {FVec("staple1zt03_11RE"), FVec("staple1zt03_11IM")},
                                  {FVec("staple1zt03_12RE"), FVec("staple1zt03_12IM")} },
                                { {FVec("staple1zt03_20RE"), FVec("staple1zt03_20IM")},
                                  {FVec("staple1zt03_21RE"), FVec("staple1zt03_21IM")},
                                  {FVec("staple1zt03_22RE"), FVec("staple1zt03_22IM")} } },
                                        { { {FVec("staple1zt04_00RE"), FVec("staple1zt04_00IM")},
                                  {FVec("staple1zt04_01RE"), FVec("staple1zt04_01IM")},
                                  {FVec("staple1zt04_02RE"), FVec("staple1zt04_02IM")} },
                                { {FVec("staple1zt04_10RE"), FVec("staple1zt04_10IM")},
                                  {FVec("staple1zt04_11RE"), FVec("staple1zt04_11IM")},
                                  {FVec("staple1zt04_12RE"), FVec("staple1zt04_12IM")} },
                                { {FVec("staple1zt04_20RE"), FVec("staple1zt04_20IM")},
                                  {FVec("staple1zt04_21RE"), FVec("staple1zt04_21IM")},
                                  {FVec("staple1zt04_22RE"), FVec("staple1zt04_22IM")} } },
                                        { { {FVec("staple1zt05_00RE"), FVec("staple1zt05_00IM")},
                                  {FVec("staple1zt05_01RE"), FVec("staple1zt05_01IM")},
                                  {FVec("staple1zt05_02RE"), FVec("staple1zt05_02IM")} },
                                { {FVec("staple1zt05_10RE"), FVec("staple1zt05_10IM")},
                                  {FVec("staple1zt05_11RE"), FVec("staple1zt05_11IM")},
                                  {FVec("staple1zt05_12RE"), FVec("staple1zt05_12IM")} },
                                { {FVec("staple1zt05_20RE"), FVec("staple1zt05_20IM")},
                                  {FVec("staple1zt05_21RE"), FVec("staple1zt05_21IM")},
                                  {FVec("staple1zt05_22RE"), FVec("staple1zt05_22IM")} } },
                                        { { {FVec("staple1zt06_00RE"), FVec("staple1zt06_00IM")},
                                  {FVec("staple1zt06_01RE"), FVec("staple1zt06_01IM")},
                                  {FVec("staple1zt06_02RE"), FVec("staple1zt06_02IM")} },
                                { {FVec("staple1zt06_10RE"), FVec("staple1zt06_10IM")},
                                  {FVec("staple1zt06_11RE"), FVec("staple1zt06_11IM")},
                                  {FVec("staple1zt06_12RE"), FVec("staple1zt06_12IM")} },
                                { {FVec("staple1zt06_20RE"), FVec("staple1zt06_20IM")},
                                  {FVec("staple1zt06_21RE"), FVec("staple1zt06_21IM")},
                                  {FVec("staple1zt06_22RE"), FVec("staple1zt06_22IM")} } },
                                        { { {FVec("staple1zt07_00RE"), FVec("staple1zt07_00IM")},
                                  {FVec("staple1zt07_01RE"), FVec("staple1zt07_01IM")},
                                  {FVec("staple1zt07_02RE"), FVec("staple1zt07_02IM")} },
                                { {FVec("staple1zt07_10RE"), FVec("staple1zt07_10IM")},
                                  {FVec("staple1zt07_11RE"), FVec("staple1zt07_11IM")},
                                  {FVec("staple1zt07_12RE"), FVec("staple1zt07_12IM")} },
                                { {FVec("staple1zt07_20RE"), FVec("staple1zt07_20IM")},
                                  {FVec("staple1zt07_21RE"), FVec("staple1zt07_21IM")},
                                  {FVec("staple1zt07_22RE"), FVec("staple1zt07_22IM")} } },
                                        { { {FVec("staple1zt10_00RE"), FVec("staple1zt10_00IM")},
                                  {FVec("staple1zt10_01RE"), FVec("staple1zt10_01IM")},
                                  {FVec("staple1zt10_02RE"), FVec("staple1zt10_02IM")} },
                                { {FVec("staple1zt10_10RE"), FVec("staple1zt10_10IM")},
                                  {FVec("staple1zt10_11RE"), FVec("staple1zt10_11IM")},
                                  {FVec("staple1zt10_12RE"), FVec("staple1zt10_12IM")} },
                                { {FVec("staple1zt10_20RE"), FVec("staple1zt10_20IM")},
                                  {FVec("staple1zt10_21RE"), FVec("staple1zt10_21IM")},
                                  {FVec("staple1zt10_22RE"), FVec("staple1zt10_22IM")} } },
                                        { { {FVec("staple1zt11_00RE"), FVec("staple1zt11_00IM")},
                                  {FVec("staple1zt11_01RE"), FVec("staple1zt11_01IM")},
                                  {FVec("staple1zt11_02RE"), FVec("staple1zt11_02IM")} },
                                { {FVec("staple1zt11_10RE"), FVec("staple1zt11_10IM")},
                                  {FVec("staple1zt11_11RE"), FVec("staple1zt11_11IM")},
                                  {FVec("staple1zt11_12RE"), FVec("staple1zt11_12IM")} },
                                { {FVec("staple1zt11_20RE"), FVec("staple1zt11_20IM")},
                                  {FVec("staple1zt11_21RE"), FVec("staple1zt11_21IM")},
                                  {FVec("staple1zt11_22RE"), FVec("staple1zt11_22IM")} } },
                                        { { {FVec("staple1zt12_00RE"), FVec("staple1zt12_00IM")},
                                  {FVec("staple1zt12_01RE"), FVec("staple1zt12_01IM")},
                                  {FVec("staple1zt12_02RE"), FVec("staple1zt12_02IM")} },
                                { {FVec("staple1zt12_10RE"), FVec("staple1zt12_10IM")},
                                  {FVec("staple1zt12_11RE"), FVec("staple1zt12_11IM")},
                                  {FVec("staple1zt12_12RE"), FVec("staple1zt12_12IM")} },
                                { {FVec("staple1zt12_20RE"), FVec("staple1zt12_20IM")},
                                  {FVec("staple1zt12_21RE"), FVec("staple1zt12_21IM")},
                                  {FVec("staple1zt12_22RE"), FVec("staple1zt12_22IM")} } },
                                        { { {FVec("staple1zt13_00RE"), FVec("staple1zt13_00IM")},
                                  {FVec("staple1zt13_01RE"), FVec("staple1zt13_01IM")},
                                  {FVec("staple1zt13_02RE"), FVec("staple1zt13_02IM")} },
                                { {FVec("staple1zt13_10RE"), FVec("staple1zt13_10IM")},
                                  {FVec("staple1zt13_11RE"), FVec("staple1zt13_11IM")},
                                  {FVec("staple1zt13_12RE"), FVec("staple1zt13_12IM")} },
                                { {FVec("staple1zt13_20RE"), FVec("staple1zt13_20IM")},
                                  {FVec("staple1zt13_21RE"), FVec("staple1zt13_21IM")},
                                  {FVec("staple1zt13_22RE"), FVec("staple1zt13_22IM")} } },
                                        { { {FVec("staple1zt14_00RE"), FVec("staple1zt14_00IM")},
                                  {FVec("staple1zt14_01RE"), FVec("staple1zt14_01IM")},
                                  {FVec("staple1zt14_02RE"), FVec("staple1zt14_02IM")} },
                                { {FVec("staple1zt14_10RE"), FVec("staple1zt14_10IM")},
                                  {FVec("staple1zt14_11RE"), FVec("staple1zt14_11IM")},
                                  {FVec("staple1zt14_12RE"), FVec("staple1zt14_12IM")} },
                                { {FVec("staple1zt14_20RE"), FVec("staple1zt14_20IM")},
                                  {FVec("staple1zt14_21RE"), FVec("staple1zt14_21IM")},
                                  {FVec("staple1zt14_22RE"), FVec("staple1zt14_22IM")} } },
                                        { { {FVec("staple1zt15_00RE"), FVec("staple1zt15_00IM")},
                                  {FVec("staple1zt15_01RE"), FVec("staple1zt15_01IM")},
                                  {FVec("staple1zt15_02RE"), FVec("staple1zt15_02IM")} },
                                { {FVec("staple1zt15_10RE"), FVec("staple1zt15_10IM")},
                                  {FVec("staple1zt15_11RE"), FVec("staple1zt15_11IM")},
                                  {FVec("staple1zt15_12RE"), FVec("staple1zt15_12IM")} },
                                { {FVec("staple1zt15_20RE"), FVec("staple1zt15_20IM")},
                                  {FVec("staple1zt15_21RE"), FVec("staple1zt15_21IM")},
                                  {FVec("staple1zt15_22RE"), FVec("staple1zt15_22IM")} } },
                                        { { {FVec("staple1zt16_00RE"), FVec("staple1zt16_00IM")},
                                  {FVec("staple1zt16_01RE"), FVec("staple1zt16_01IM")},
                                  {FVec("staple1zt16_02RE"), FVec("staple1zt16_02IM")} },
                                { {FVec("staple1zt16_10RE"), FVec("staple1zt16_10IM")},
                                  {FVec("staple1zt16_11RE"), FVec("staple1zt16_11IM")},
                                  {FVec("staple1zt16_12RE"), FVec("staple1zt16_12IM")} },
                                { {FVec("staple1zt16_20RE"), FVec("staple1zt16_20IM")},
                                  {FVec("staple1zt16_21RE"), FVec("staple1zt16_21IM")},
                                  {FVec("staple1zt16_22RE"), FVec("staple1zt16_22IM")} } },
                                        { { {FVec("staple1zt17_00RE"), FVec("staple1zt17_00IM")},
                                  {FVec("staple1zt17_01RE"), FVec("staple1zt17_01IM")},
                                  {FVec("staple1zt17_02RE"), FVec("staple1zt17_02IM")} },
                                { {FVec("staple1zt17_10RE"), FVec("staple1zt17_10IM")},
                                  {FVec("staple1zt17_11RE"), FVec("staple1zt17_11IM")},
                                  {FVec("staple1zt17_12RE"), FVec("staple1zt17_12IM")} },
                                { {FVec("staple1zt17_20RE"), FVec("staple1zt17_20IM")},
                                  {FVec("staple1zt17_21RE"), FVec("staple1zt17_21IM")},
                                  {FVec("staple1zt17_22RE"), FVec("staple1zt17_22IM")} } } } 

 					};

static FVec Herm1[9] = { FVec("Herm1_0"), FVec("Herm1_1"), FVec("Herm1_2"), FVec("Herm1_3"), FVec("Herm1_4"), FVec("Herm1_5"), FVec("Herm1_6"), FVec("Herm1_7"), FVec("Herm1_8") }; 
static FVec Herm2[9] = { FVec("Herm2_0"), FVec("Herm2_1"), FVec("Herm2_2"), FVec("Herm2_3"), FVec("Herm2_4"), FVec("Herm2_5"), FVec("Herm2_6"), FVec("Herm2_7"), FVec("Herm2_8") };
static FVec Herm3[9] = { FVec("Herm3_0"), FVec("Herm3_1"), FVec("Herm3_2"), FVec("Herm3_3"), FVec("Herm3_4"), FVec("Herm3_5"), FVec("Herm3_6"), FVec("Herm3_7"), FVec("Herm3_8") };
static FVec Herm17W[7][9] = { { FVec("Herm10_0"), FVec("Herm10_1"), FVec("Herm10_2"), FVec("Herm10_3"), FVec("Herm10_4"), FVec("Herm10_5"), FVec("Herm10_6"), FVec("Herm10_7"), FVec("Herm10_8") },
				{ FVec("Herm11_0"), FVec("Herm11_1"), FVec("Herm11_2"), FVec("Herm11_3"), FVec("Herm11_4"), FVec("Herm11_5"), FVec("Herm11_6"), FVec("Herm11_7"), FVec("Herm11_8") },
				{ FVec("Herm12_0"), FVec("Herm12_1"), FVec("Herm12_2"), FVec("Herm12_3"), FVec("Herm12_4"), FVec("Herm12_5"), FVec("Herm12_6"), FVec("Herm12_7"), FVec("Herm12_8") },
				{ FVec("Herm13_0"), FVec("Herm13_1"), FVec("Herm13_2"), FVec("Herm13_3"), FVec("Herm13_4"), FVec("Herm13_5"), FVec("Herm13_6"), FVec("Herm13_7"), FVec("Herm13_8") },
				{ FVec("Herm14_0"), FVec("Herm14_1"), FVec("Herm14_2"), FVec("Herm14_3"), FVec("Herm14_4"), FVec("Herm14_5"), FVec("Herm14_6"), FVec("Herm14_7"), FVec("Herm14_8") },
				{ FVec("Herm15_0"), FVec("Herm15_1"), FVec("Herm15_2"), FVec("Herm15_3"), FVec("Herm15_4"), FVec("Herm15_5"), FVec("Herm15_6"), FVec("Herm15_7"), FVec("Herm15_8") },
				{ FVec("Herm16_0"), FVec("Herm16_1"), FVec("Herm16_2"), FVec("Herm16_3"), FVec("Herm16_4"), FVec("Herm16_5"), FVec("Herm16_6"), FVec("Herm16_7"), FVec("Herm16_8") } };

static FVec Herm27W[7][9] = { { FVec("Herm20_0"), FVec("Herm20_1"), FVec("Herm20_2"), FVec("Herm20_3"), FVec("Herm20_4"), FVec("Herm20_5"), FVec("Herm20_6"), FVec("Herm20_7"), FVec("Herm20_8") },                                { FVec("Herm21_0"), FVec("Herm21_1"), FVec("Herm21_2"), FVec("Herm21_3"), FVec("Herm21_4"), FVec("Herm21_5"), FVec("Herm21_6"), FVec("Herm21_7"), FVec("Herm21_8") },
                                { FVec("Herm22_0"), FVec("Herm22_1"), FVec("Herm22_2"), FVec("Herm22_3"), FVec("Herm22_4"), FVec("Herm22_5"), FVec("Herm22_6"), FVec("Herm22_7"), FVec("Herm22_8") },
                                { FVec("Herm23_0"), FVec("Herm23_1"), FVec("Herm23_2"), FVec("Herm23_3"), FVec("Herm23_4"), FVec("Herm23_5"), FVec("Herm23_6"), FVec("Herm23_7"), FVec("Herm23_8") },
                                { FVec("Herm24_0"), FVec("Herm24_1"), FVec("Herm24_2"), FVec("Herm24_3"), FVec("Herm24_4"), FVec("Herm24_5"), FVec("Herm24_6"), FVec("Herm24_7"), FVec("Herm24_8") },
                                { FVec("Herm25_0"), FVec("Herm25_1"), FVec("Herm25_2"), FVec("Herm25_3"), FVec("Herm25_4"), FVec("Herm25_5"), FVec("Herm25_6"), FVec("Herm25_7"), FVec("Herm25_8") },
                                { FVec("Herm26_0"), FVec("Herm26_1"), FVec("Herm26_2"), FVec("Herm26_3"), FVec("Herm26_4"), FVec("Herm26_5"), FVec("Herm26_6"), FVec("Herm26_7"), FVec("Herm26_8") } };

static
int index(int dir, int dir2)
{
  int dt = (dir<dir2?dir:dir2);
  if(dt==0) return (dir+dir2 - 1);
  else return (dir+dir2);
}

static inline 
int upDir(int dir)
{
  return (2*dir+1);
}

static inline 
int dnDir(int dir)
{
  return (2*dir);
}

static void zeroGResult(InstVector& ivector, vMatrix out, bool compress12)
{
    for(int r=0; r<3-compress12; ++r) for(int c=0; c < 3; c++) for(int ir=0; ir<2; ++ir) {
        setZero(ivector, out[r][c][ir]);
    }
}

static void loadGaugeDirXYZT(InstVector& ivector, vMatrix matr, string inGbase, int dirg, int dir, bool compress12, string mask)
{
  loadGaugeDir(ivector, matr, inGbase, dir, compress12, mask);
  int d = ((1 << dirg));
  string isBoundary("isBoundary");
  stringstream n_str;
  n_str << d;
  ifStatement(ivector, isBoundary + " & " + n_str.str());
  {
    for(int r=0; r<3-compress12; ++r) for(int c=0; c<3; ++c) for(int ir=0; ir<2; ++ir)
      permuteXYZTFVec(ivector, matr[r][c][ir], matr[r][c][ir], dirg);
  }
  endScope(ivector);
}

static void loadGauge7WayDirXYZT(InstVector& ivector, vMatrix *matr, string inGbase, int dir, bool compress12, string mask)
{
  for(int i=0, n=0; i<8; ++i) if(i!=OPPDIR(dir)) {
    loadGaugeDir(ivector, matr[n], inGbase, i, compress12, mask);
    n++;
  }
  string isBoundary("isBoundary");
  stringstream n_str;
  n_str << (1 << dir);
  ifStatement(ivector, isBoundary + " & " + n_str.str());
  {
    for(int i=0; i<7; ++i) for(int r=0; r<3-compress12; ++r) for(int c=0; c<3; ++c) for(int ir=0; ir<2; ++ir)
	permuteXYZTFVec(ivector, matr[i][r][c][ir], matr[i][r][c][ir], dir);
  }
  endScope(ivector);

}

static void loadHermitDirXYZT(InstVector& ivector, FVec *herm, string inHbase, int dirh, int dir, string mask)
{
  loadHermitDir(ivector, herm, inHbase, dir, mask);
  int d = ((1 << dirh));
  string isBoundary("isBoundary");
  stringstream n_str;
  n_str << d;
  ifStatement(ivector, isBoundary + " & " + n_str.str());
  {
    for(int k=0; k<8; ++k)
	permuteXYZTFVec(ivector, herm[k], herm[k], dirh);
  }
  endScope(ivector);
  naddFVec(ivector, herm[8], herm[6], herm[7], mask);
}

static void loadHermit7WayDirXYZT(InstVector& ivector, FVec *herm[9], string inHbase, int dir, string mask)
{
  int d0=0;
  for(int d=0; d<8; ++d) if(d!=OPPDIR(dir))
    loadHermitDir(ivector, herm[d0++], inHbase, d, mask);
  int d = (1 << dir);
  string isBoundary("isBoundary");
  stringstream n_str;
  n_str << d;
  ifStatement(ivector, isBoundary + " & " + n_str.str());
  {
    for(int dd=0; dd<7; ++dd) for(int k=0; k<8; ++k)
	permuteXYZTFVec(ivector, herm[dd][k], herm[dd][k], dir);
  }
  endScope(ivector);
  for(int dd=0; dd<7; ++dd)
	naddFVec(ivector, herm[dd][8], herm[dd][6], herm[dd][7], mask);
}

static void storeHermitDirXYZT(InstVector& ivector, string outHbase, FVec *herm, int dirh, int dir, string mask, bool isstream)
{
  int d = ((1 << dirh));
  string isBoundary("isBoundary");
  stringstream n_str;
  n_str << d;
  ifStatement(ivector, isBoundary + " & " + n_str.str());
  {
    for(int k=0; k<8; ++k)
	aPermuteXYZTFVec(ivector, herm[k], herm[k], dirh);
  }
  endScope(ivector);
  storeHermitDir(ivector, outHbase, herm, dir, mask, isstream);
}

static int OPPDIR(int dir)
{
  if(dir%2==1)
    return dir-1;
  else
    return dir+1;
}

/*
static void loadHermitHelperDirXYZT(InstVector& ivector, vHerm herm, int dir1, int dir2, string isBoundary, string mask)
{
  loadHermitHelperDir(ivector, herm, inHHBase[dir1], OPPDIR(dir1), dir2, mask);
    int d = ((1 << dir1));
  stringstream n_str;
  n_str << d;
  ifStatement(ivector, isBoundary + " & " + n_str.str());
  {
    for(int k=0; k<9; ++k) 
      permuteXYZTFVec(ivector, herm[k], herm[k], dir1);
  }
  endScope(ivector);
}
*/

static void aPermuteXYZTGauge(InstVector& ivector, vMatrix matr, int dir, string isBoundary, bool compress12, string mask)
{
  stringstream n_str;
  n_str << ( (1 << dir) );
  ifStatement(ivector, isBoundary + " & " + n_str.str());
  {
    for(int r=0; r<3-compress12; ++r) for(int c=0; c<3; ++c) for(int ir=0; ir<2; ++ir)
      aPermuteXYZTFVec(ivector, matr[r][c][ir], matr[r][c][ir], dir);
  }
  endScope(ivector);
}

/*
 * Out : (0,0,0,0), (2,0,0,0), (1,-1,0,0), (1,1,0,0), (1,0,-1,0), (1,0,1,0), (1,0,0,-1), (1,0,0,1).
 */

/* 
 * dir :     X          X              Y          Y             Z           Z             T           T
 * Out : (0,0,0,0), (2,0,0,0), or (1,-1,0,0), (1,1,0,0), or (1,0,-1,0), (1,0,1,0), or (1,0,0,-1), (1,0,0,1)
 *                    ____       __   
 * Out1 ([2*dir]) += |  __| + |____| ,
 *                      ____     __
 * Out2 ([2*dir+1]) += |__  | + |____| .
 */
static
void g_force_1( InstVector& ivector, vMatrix Out1, vMatrix Out2, string &kappa, int dir, bool compress12, int add=1, bool load=1 )
{
/*                ____
 * Rectangular : |    | + |____| .
 */
  //string accumulate("accumulate");
  if(load)
  {
  loadGaugeDirXYZT(ivector, tmpMtr2, inGBase[2*dir+1], 2*dir+1, 2*dir, compress12, mask);
  loadGaugeDirXYZT(ivector, tmpMtr3, inGBase[2*dir], 2*dir, 2*dir+1, compress12, mask);
  }
  zeroGResult(ivector, tmpMtr, compress12);
  for(int dir2 = 0; dir2<LDim; ++dir2) if(dir2 != dir)
  {
    /* Up */
    stringstream d_str;
    d_str << ((1 << (2*dir2+1)) );
    //ifStatement(ivector, accumulate + " & " + d_str.str());
    {
      rectangle0(ivector, tmpMtr, staple1[index(dir,dir2)][8*(dir2<dir)+2], staple1[index(dir,dir2)][8*(dir2<dir)+3], 1, mask);
    }
    //endScope(ivector);
    /* Down */
    stringstream d2_str;
    d2_str << ((1 << (2*dir2)) );
    //ifStatement(ivector, accumulate + " & " + d2_str.str());
    {
      rectangle0(ivector, tmpMtr, staple1[index(dir,dir2)][8*(dir2<dir)+4], staple1[index(dir,dir2)][8*(dir2<dir)+5], 1, mask); 
    }
    //endScope(ivector);
  }
  //rectangle1(ivector, Out1, Out2, tmpMtr, In[2*dir+1][2*dir], In[2*dir][2*dir], add, compress12, mask);
  if(load) rectangle1(ivector, Out1, Out2, tmpMtr, tmpMtr2, tmpMtr3, kappa, add, mask);
  else rectangle1(ivector, Out1, Out2, tmpMtr, In[2*dir+1][2*dir], In[2*dir][2*dir], kappa, add, mask);
}

/*
 * dir :    	  !X			    !Y			      !Z			!T
 * Out : (0,0,0,0), (2,0,0,0), or (1,-1,0,0), (1,1,0,0), or (1,0,-1,0), (1,0,1,0), or (1,0,0,-1), (1,0,0,1)
 */ 
static 
void g_force_2( InstVector& ivector, vMatrix Out1, vMatrix Out2, string &kappar, string &kappab, int dir, int dir2, int updown, bool compress12, int add=1, bool load=1 )
{
/*
 * Rectangulars with only up or down staples
 */
    //string accumulate("accumulate");
    stringstream d1_str, d2_str, d3_str;
    int d = ((1 << (2*dir)) );
    d1_str << d;
    d = ((1 << (2*dir+1)) );
    d2_str << d;
    if(updown==UP) d = ((1 << (2*dir2)) );
    else d = ((1 << (2*dir2+1)) );
    d3_str << d;
    //ifStatement(ivector, accumulate + " & " + d1_str.str());
    //ifStatement(ivector, accumulate + " & " + d2_str.str());
    {
      /* Up staple */
      if(updown==UP)
      {
	if(load) {
            loadGaugeDirXYZT(ivector, tmpMtr2, inGBase[2*dir2], 2*dir2, upDir(dir), compress12, mask);
            loadGaugeDirXYZT(ivector, tmpMtr3, inGBase[2*dir2], 2*dir2, dnDir(dir), compress12, mask);
        }
	rectangle0(ivector, tmpMtr, staple1[index(dir,dir2)][8*(dir2<dir)], staple1[index(dir,dir2)][8*(dir2<dir)+1], 0, mask);
        //rectangle1(ivector, Out1, Out2, tmpMtr, In[2*dir2][2*dir+1-(dir>dir2)], In[2*dir2][2*dir-(dir>dir2)], add, compress12, mask);
        if(load) rectangle1(ivector, Out1, Out2, tmpMtr, tmpMtr2, tmpMtr3, kappar, add, mask);
	else rectangle1(ivector, Out1, Out2, tmpMtr, In[2*dir2][upDir(dir)-(dir>dir2)], In[2*dir2][dnDir(dir)-(dir>dir2)], kappar, add, mask);
#if 0
        loadGaugeDir(ivector, tmpMtr2, inGBase[2*dir], dnDir(dir2), compress12, mask);
        loadGaugeDir(ivector, tmpMtr3, inGBase[2*dir2], upDir(dir2), compress12, mask);
        //rectangle2(ivector, Out1, In[2*dir][2*dir2-(dir2>dir)], staple1[index(dir,dir2)][8*(dir2<dir)+2], In[2*dir2][2*dir2], UP, compress12, mask);
        rectangle2(ivector, Out1, tmpMtr2, staple1[index(dir,dir2)][8*(dir2<dir)+2], tmpMtr3, tmpMtr, UP, mask);
        loadGaugeDir(ivector, tmpMtr2, inGBase[2*dir+1], dnDir(dir2),  compress12, mask);
        //rectangle2(ivector, Out2, In[2*dir2][2*dir2], staple1[index(dir,dir2)][8*(dir2<dir)+3], In[2*dir+1][2*dir2-(dir2>dir)], UP, compress12, mask);
        rectangle2(ivector, Out2, tmpMtr3, staple1[index(dir,dir2)][8*(dir2<dir)+3], tmpMtr2, tmpMtr, UP, mask);
        add=1;
#endif
      }
      /* Down staple */
      else
      {
	if(load) {
            loadGaugeDirXYZT(ivector, tmpMtr2, inGBase[2*dir2+1], 2*dir2+1, upDir(dir), compress12, mask);
            loadGaugeDirXYZT(ivector, tmpMtr3, inGBase[2*dir2+1], 2*dir2+1, dnDir(dir), compress12, mask);
        }
	rectangle0(ivector, tmpMtr, staple1[index(dir,dir2)][8*(dir2<dir)+6], staple1[index(dir,dir2)][8*(dir2<dir)+7], 0, mask);
        //rectangle1(ivector, Out1, Out2, tmpMtr, In[2*dir2+1][2*dir+1-(dir>dir2)], In[2*dir2+1][2*dir-(dir>dir2)], add, compress12, mask);
        if(load) rectangle1(ivector, Out1, Out2, tmpMtr, tmpMtr2, tmpMtr3, kappar, add, mask);
	else rectangle1(ivector, Out1, Out2, tmpMtr, In[2*dir2+1][upDir(dir)-(dir>dir2)], In[2*dir2+1][dnDir(dir)-(dir>dir2)], kappar, add, mask);
#if 0
        loadGaugeDir(ivector, tmpMtr2, inGBase[2*dir], upDir(dir2), compress12, mask);
        loadGaugeDir(ivector, tmpMtr3, inGBase[2*dir2+1], dnDir(dir2), compress12, mask);
      //rectangle2(ivector, Out1, In[2*dir][2*dir2+1-(dir2>dir)], staple1[index(dir,dir2)][8*(dir2<dir)+4], In[2*dir2+1][2*dir2], DN, compress12, mask);
        rectangle2(ivector, Out1, tmpMtr2, staple1[index(dir,dir2)][8*(dir2<dir)+4], tmpMtr3, tmpMtr, DN, mask);
        loadGaugeDir(ivector, tmpMtr2, inGBase[2*dir+1], upDir(dir2), compress12, mask);
      //rectangle2(ivector, Out2, In[2*dir2+1][2*dir2], staple1[index(dir,dir2)][8*(dir2<dir)+5], In[2*dir+1][2*dir2+1-(dir2>dir)], DN, compress12, mask);
        rectangle2(ivector, Out2, tmpMtr3, staple1[index(dir,dir2)][8*(dir2<dir)+5], tmpMtr2, tmpMtr, DN, mask);
        add=1;
#endif
      }
    }
    //endScope(ivector);
    //endScope(ivector);

    stringstream d4_str, d5_str, d6_str;
    d = ((1 << (2*dir2)) );
    d4_str << d;
    d = ((1 << (2*dir2+1)) );
    d5_str << d;
    //ifStatement(ivector, accumulate + " & " + d4_str.str());
    //ifStatement(ivector, accumulate + " & " + d5_str.str());
    {
      //ifStatement(ivector, accumulate + " & " + d1_str.str());
      {
        //loadGaugeDirXYZT(ivector, tmpMtr2, inGBase[2*dir], dnDir(dir2), compress12, mask);
        if(updown==UP) {
	  if(load) {
	      loadGaugeDirXYZT(ivector, tmpMtr2, inGBase[2*dir], 2*dir, dnDir(dir2), compress12, mask);
	      loadGaugeDirXYZT(ivector, tmpMtr3, inGBase[2*dir2], 2*dir2, upDir(dir2), compress12, mask);
              rectangle2(ivector, Out1, tmpMtr2, staple1[index(dir,dir2)][8*(dir2<dir)+2], tmpMtr3, tmpMtr, kappar, UP, mask);
	  }
	  else rectangle2(ivector, Out1, In[2*dir][dnDir(dir2)-(dir2>dir)], staple1[index(dir,dir2)][8*(dir2<dir)+2], In[2*dir2][dnDir(dir2)], tmpMtr, kappar, UP, mask);
	}
	else {
	  if(load) {
	      loadGaugeDirXYZT(ivector, tmpMtr2, inGBase[2*dir], 2*dir, upDir(dir2), compress12, mask);
	      loadGaugeDirXYZT(ivector, tmpMtr3, inGBase[2*dir2+1], 2*dir2+1, dnDir(dir2), compress12, mask);
	      rectangle2(ivector, Out1, tmpMtr2, staple1[index(dir,dir2)][8*(dir2<dir)+4], tmpMtr3, tmpMtr, kappar, DN, mask);
	  }
	  else rectangle2(ivector, Out1, In[2*dir][upDir(dir2)-(dir2>dir)], staple1[index(dir,dir2)][8*(dir2<dir)+4], In[2*dir2+1][dnDir(dir2)], tmpMtr, kappar, DN, mask);
	}
        //loadGaugeDir(ivector, tmpMtr2, inGBase[2*dir+1], dnDir(dir2),  compress12, mask);
        //rectangle2(ivector, Out2, tmpMtr3, staple1[index(dir,dir2)][8*(dir2<dir)+3], tmpMtr2, tmpMtr, UP, mask);
      }
      //endScope(ivector);
      //ifStatement(ivector, accumulate + " & " + d2_str.str());
      {
	if(updown==UP) {
	  if(load) {
	      loadGaugeDirXYZT(ivector, tmpMtr3, inGBase[2*dir2], 2*dir2, upDir(dir2), compress12, mask);
	      loadGaugeDirXYZT(ivector, tmpMtr2, inGBase[2*dir+1], 2*dir+1, dnDir(dir2),  compress12, mask);
	      rectangle2(ivector, Out2, tmpMtr3, staple1[index(dir,dir2)][8*(dir2<dir)+3], tmpMtr2, tmpMtr, kappar, UP, mask);
	  }
	  else rectangle2(ivector, Out2, In[2*dir2][dnDir(dir2)], staple1[index(dir,dir2)][8*(dir2<dir)+3], In[2*dir+1][dnDir(dir2)-(dir2>dir)], tmpMtr, kappar, UP, mask);
        }
	else {
	  if(load) {
	      loadGaugeDirXYZT(ivector, tmpMtr3, inGBase[2*dir2+1], 2*dir2+1, dnDir(dir2), compress12, mask);
	      loadGaugeDirXYZT(ivector, tmpMtr2, inGBase[2*dir+1], 2*dir+1, upDir(dir2),  compress12, mask);
	      rectangle2(ivector, Out2, tmpMtr3, staple1[index(dir,dir2)][8*(dir2<dir)+5], tmpMtr2, tmpMtr, kappar, DN, mask);
	  }
	  else rectangle2(ivector, Out2, In[2*dir2+1][dnDir(dir2)], staple1[index(dir,dir2)][8*(dir2<dir)+5], In[2*dir+1][upDir(dir2)-(dir2>dir)], tmpMtr, kappar, DN, mask);
	}
      }
      //endScope(ivector);
    }
    //endScope(ivector);
    //endScope(ivector);

/*
 *  Bentchairs
 */    
    //int ff=0;
    /* Out1 */
    //ifStatement(ivector, accumulate + " & " + d1_str.str());
    {
      zeroGResult(ivector, tmpMtr, compress12);
      for(int dir3=0; dir3<LDim; ++dir3) if(dir3!=dir && dir3!=dir2) {
	stringstream d8_str, d9_str;
        int d = ((1 << (2*dir3+1)) );
        d8_str << d;
	d = ((1 << (2*dir3)) );
	d9_str << d;
	//ifStatement(ivector, accumulate + " & " + d8_str.str());
	{
	  bentchair0(ivector, tmpMtr, staple1[index(dir,dir3)][8*(dir>dir3)+2], staple1[index(dir2,dir3)][8*(dir2>dir3)+2+(updown==DN)], false, (updown==UP), 1, mask);
	}
	//endScope(ivector);
	//ifStatement(ivector, accumulate + " & " + d9_str.str());
	{
	  bentchair0(ivector, tmpMtr, staple1[index(dir,dir3)][8*(dir>dir3)+4], staple1[index(dir2,dir3)][8*(dir2>dir3)+4+(updown==DN)], false, (updown==UP), 1, mask);
	}
	//endScope(ivector);
      //bentchair0(ivector, tmpMtr2, staple1[index(dir2,dir3)][8*(dir2>dir3)+2+(updown==UP)], staple1[index(dir,dir3)][8*(dir>dir3)+3], (updown==DN), false, ff, mask);
      //bentchair0(ivector, tmpMtr2, staple1[index(dir2,dir3)][8*(dir2>dir3)+4+(updown==UP)], staple1[index(dir,dir3)][8*(dir>dir3)+5], (updown==DN), false, 1, mask);
      //ff=1;
      } 
      if(load) {
      	loadGaugeDirXYZT(ivector, tmpMtr3, inGBase[2*dir], 2*dir, 2*dir2+(updown==DN), compress12, mask);
    //addAdjMatMultMat(ivector, Out1, Out1, In[2*dir][2*dir2+(updown==UP)-(dir2>dir)], tmpMtr, (updown==UP), compress12, mask);
      	addSMulAdjMatMultMat(ivector, Out1, Out1, tmpMtr3, tmpMtr, kappab, (updown==DN), mask);
      }
      else addSMulAdjMatMultMat(ivector, Out1, Out1, In[2*dir][2*dir2+(updown==DN)-(dir2>dir)], tmpMtr, kappab, (updown==DN), mask);
    }
    //endScope(ivector);
    //ff=0;
    /* Out2 */
    //ifStatement(ivector, accumulate + " & " + d2_str.str());
    {
      zeroGResult(ivector, tmpMtr, compress12);
      for(int dir3=0; dir3<LDim; ++dir3) if(dir3!=dir && dir3!=dir2) {
        stringstream d8_str, d9_str;
        int d = ((1 << (2*dir3+1)) );
        d8_str << d;
        d = ((1 << (2*dir3)) );
        d9_str << d;
        //ifStatement(ivector, accumulate + " & " + d8_str.str());
        {
	  bentchair0(ivector, tmpMtr, staple1[index(dir2,dir3)][8*(dir2>dir3)+2+(updown==DN)], staple1[index(dir,dir3)][8*(dir>dir3)+3], (updown==DN), false, 1, mask);
	}
	//endScope(ivector);
	//ifStatement(ivector, accumulate + " & " + d9_str.str());
	{
	  bentchair0(ivector, tmpMtr, staple1[index(dir2,dir3)][8*(dir2>dir3)+4+(updown==DN)], staple1[index(dir,dir3)][8*(dir>dir3)+5], (updown==DN), false, 1, mask);
      //ff=1;
	}
	//endScope(ivector);
      }
      if(load) {
      	loadGaugeDirXYZT(ivector, tmpMtr3, inGBase[2*dir+1], 2*dir+1, 2*dir2+(updown==DN), compress12, mask);
    //addMatMultMat(ivector, Out2, Out2, tmpMtr2, In[2*dir+1][2*dir2+(updown==UP)-(dir2>dir)], (updown==DN), compress12, mask);
      	addSMulMatMultMat(ivector, Out2, Out2, tmpMtr, tmpMtr3, kappab, (updown==UP), mask);
      }
      else addSMulMatMultMat(ivector, Out2, Out2, tmpMtr, In[2*dir+1][2*dir2+(updown==UP)-(dir2>dir)], kappab, (updown==UP), mask);
    }
    //endScope(ivector);

} 

/* 
 * Calculate gauge forces,
 * and update H. No updates on U.
 *
 * For leapfrog, 2G1F, 3G1F. 
 * //Output and input U are stored separately.
 *
 */
void g_force_imp_body(InstVector& ivector, bool compress12)
{
  string face[6] = {"xy", "xz", "xt", "yz", "yt", "zt"};
  string loop[16] = { "00", "01", "02", "03", "04", "05", "06", "07", "10", "11", "12", "13", "14", "15", "16", "17"};
  for(int i=0; i<6; ++i) for(int j=0; j<16; ++j)
  {
    declare_Gauge(ivector, staple1[i][j]);
  }
  declare_Gauge(ivector, tmpMtr);
  declare_Gauge(ivector, tmpMtr2);
  declare_Gauge(ivector, tmpMtr3);
  declare_Gauge(ivector, Out1);
  declare_Gauge(ivector, Out2);
  declare_Hermit(ivector, Herm1);

  /* Load gauges from inGBase to In */
#if 0
  for(int dir1=0; dir1<2*LDim; ++dir1) for(int dir2=0; dir2<2*LDim; ++dir2) {
    loadGaugeDir(ivector, In[dir1][dir2-(dir2/2>dir1/2)], inGBase[dir1], ((dir1/2==dir2/2) ? dir2+(dir1+1)%2 : dir2), compress12, mask);
    if(dir1/2==dir2/2) dir2++;
  }
#endif

  //string accumulate("accumulate");
  string isBoundary("isBoundary");

  /***********/
  /* Staples */
  /***********/
  for(int dir1=0; dir1<LDim; ++dir1) for(int dir2=dir1+1; dir2<LDim; ++dir2) 
  {
    for(int n=0; n<4; ++n) 
    {
      stringstream d1_str, d2_str;
      int d = ((1 << (2*dir1+n%2)) );
      d1_str << d;
      d = ((1 << (2*dir2+n/2)) );
      d2_str << d;
      //ifStatement(ivector, accumulate + " & " + d1_str.str());
      //ifStatement(ivector, accumulate + " & " + d2_str.str());
      {
	loadGaugeDir(ivector, tmpMtr, inGBase[2*dir1+n%2], 2*dir1+(n+1)%2, compress12, mask);
	loadGaugeDir(ivector, tmpMtr2, inGBase[2*dir1+n%2], 2*dir2+n/2, compress12, mask);
	ifStatement(ivector, isBoundary + " & " + d1_str.str());
	{
	  for(int r=0; r<3-compress12; ++r) for(int c=0; c<3; ++c) for(int ir=0; ir<2; ++ir)
	    permuteXYZTFVec(ivector, tmpMtr[r][c][ir], tmpMtr[r][c][ir], 2*dir1+n%2);
	  for(int r=0; r<3-compress12; ++r) for(int c=0; c<3; ++c) for(int ir=0; ir<2; ++ir)
	    permuteXYZTFVec(ivector, tmpMtr2[r][c][ir], tmpMtr2[r][c][ir], 2*dir1+n%2);
	}
	endScope(ivector);
	loadGaugeDir(ivector, Out1, inGBase[2*dir2+n/2], 2*dir1+n%2, compress12, mask);
	loadGaugeDir(ivector, Out2, inGBase[2*dir2+n/2], 2*dir2+1-n/2, compress12, mask);
	ifStatement(ivector, isBoundary + " & " + d2_str.str());
	{
	  for(int r=0; r<3-compress12; ++r) for(int c=0; c<3; ++c) for(int ir=0; ir<2; ++ir)
	    permuteXYZTFVec(ivector, Out1[r][c][ir], Out1[r][c][ir], 2*dir2+n/2);
	  for(int r=0; r<3-compress12; ++r) for(int c=0; c<3; ++c) for(int ir=0; ir<2; ++ir)
	    permuteXYZTFVec(ivector, Out2[r][c][ir], Out2[r][c][ir], 2*dir2+n/2);
	}
	endScope(ivector);
        if(n==0)
        {
          staple(ivector, staple1[index(dir1,dir2)][n], tmpMtr2, tmpMtr, Out2, UP, mask);
	  staple(ivector, staple1[index(dir1,dir2)][n+4], tmpMtr2, Out1, Out2, DN, mask);
	  staple(ivector, staple1[index(dir1,dir2)][n+8], Out1, Out2, tmpMtr, UP, mask);
	  staple(ivector, staple1[index(dir1,dir2)][n+12], Out1, tmpMtr2, tmpMtr, DN, mask);
        }
        else if(n==1)
        {
	  staple(ivector, staple1[index(dir1,dir2)][n], Out2, tmpMtr, tmpMtr2, UP, mask);
	  staple(ivector, staple1[index(dir1,dir2)][n+4], Out2, Out1, tmpMtr2, DN, mask);
	  staple(ivector, staple1[index(dir1,dir2)][n+9], Out1, tmpMtr2, tmpMtr, UP, mask);
	  staple(ivector, staple1[index(dir1,dir2)][n+13], Out1, Out2, tmpMtr, DN, mask);
        }
	else if(n==2)
	{
	  staple(ivector, staple1[index(dir1,dir2)][n], tmpMtr2, Out1, Out2, UP, mask);
	  staple(ivector, staple1[index(dir1,dir2)][n+4], tmpMtr2, tmpMtr, Out2, DN, mask);
	  staple(ivector, staple1[index(dir1,dir2)][n+7], tmpMtr, Out2, Out1, UP, mask);
	  staple(ivector, staple1[index(dir1,dir2)][n+11], tmpMtr, tmpMtr2, Out1, DN, mask);
	}
	else
	{
	  staple(ivector, staple1[index(dir1,dir2)][n], Out2, Out1, tmpMtr2, UP, mask);
	  staple(ivector, staple1[index(dir1,dir2)][n+4], Out2, tmpMtr, tmpMtr2, DN, mask);
	  staple(ivector, staple1[index(dir1,dir2)][n+8], tmpMtr, tmpMtr2, Out1, UP, mask);
	  staple(ivector, staple1[index(dir1,dir2)][n+12], tmpMtr, Out2, Out1, DN, mask);
	}
      }
      //endScope(ivector);
      //endScope(ivector);
    } /* n loop */
    
  }

  /* Loop over directions of gauge offsets */
  for(int dir1=0; dir1<LDim; ++dir1) {
    /* Loop over dirction of gauge links */
    for(int dir2=0; dir2<LDim; ++dir2) {
      /* direction of gauge offset is the same as gauge direction */
      if(dir2==dir1) {
        zeroGResult(ivector, Out1, compress12);
        zeroGResult(ivector, Out2, compress12);
	/**************************/
	/* Calculate gauge forces */
	/**************************/
        /* From staples */
        stringstream d1_str, d2_str;
	int d = ((1 << (2*dir1)) );
	d1_str << d;
	d = ((1 << (2*dir1+1)) );
	d2_str << d;
	//ifStatement(ivector, accumulate + " & " + d1_str.str());
	{
          for(int dir3=0; dir3<LDim; ++dir3) if(dir3!=dir2) {
	    stringstream d3_str;
	    d = ((1 << (2*dir3+1)) );
	    d3_str << d;
	    //ifStatement(ivector, accumulate + " & " + d3_str.str());
	    {
              /* Up */
              addSMulMat(ivector, Out1, Out1, staple1[index(dir2,dir3)][8*(dir3<dir2)+2], kappaS, mask);
	    }
	    //endScope(ivector);
	    d = ((1 << (2*dir3)) );
	    stringstream d4_str;
	    d4_str << d;
            //ifStatement(ivector, accumulate + " & " + d4_str.str());
	    {
	      /* Down */
              addSMulMat(ivector, Out1, Out1, staple1[index(dir2,dir3)][8*(dir3<dir2)+4], kappaS, mask);
	    }
	    //endScope(ivector);
	  }
	}
	//endScope(ivector);
	//ifStatement(ivector, accumulate + " & " + d2_str.str());
	{
	  for(int dir3=0; dir3<LDim; ++dir3) if(dir3!=dir2) {
            stringstream d3_str;
            d = ((1 << (2*dir3+1)) );
            d3_str << d;
	    //ifStatement(ivector, accumulate + " & " + d3_str.str());
            {
              /* Up */
	      addSMulMat(ivector, Out2, Out2, staple1[index(dir2,dir3)][8*(dir3<dir2)+3], kappaS, mask);
	    }
	    //endScope(ivector);
            d = ((1 << (2*dir3)) );
            stringstream d4_str;
            d4_str << d;
            //ifStatement(ivector, accumulate + " & " + d4_str.str());
            {
              /* Down */
	      addSMulMat(ivector, Out2, Out2, staple1[index(dir2,dir3)][8*(dir3<dir2)+5], kappaS, mask);
            }
            //endScope(ivector);
          }
        }
	//endScope(ivector);

        /* From rectangulars */
	//ifStatement(ivector, accumulate + " & " + d1_str.str());
	//ifStatement(ivector, accumulate + " & " + d2_str.str());
	{
          g_force_1(ivector, Out1, Out2, kappaR, dir2, compress12);
	}
	//endScope(ivector);
	//endScope(ivector);

	/********************/
	/* Update momentum  */
	/********************/
	/* Backward */
	//ifStatement(ivector, accumulate + " & " + d1_str.str());
	{
	loadGaugeDirXYZT(ivector, tmpMtr2, inGBase[2*dir1], 2*dir1, upDir(dir2), compress12, mask);
	if(dir1==0) {
		loadHermitDirXYZT(ivector, Herm1, outHBase[dnDir(dir1)], dnDir(dir1), upDir(dir2), mask);
	}
	//if(dir1==0) loadHermitDir(ivector, Herm1, outHBase[dnDir(dir1)], upDir(dir2), mask);
	//else loadHermitDirXYZT(ivector, Herm1, outHBase[dnDir(dir1)], dnDir(dir1), dnDir(dir2), mask);
	matMultMatSconj(ivector, tmpMtr, tmpMtr2, Out1, mask);
	matAntiHermTrls(ivector, &tmpMtr2[0][0][0], tmpMtr, mask);
	if(dir1==0) {
		addSMultAntiHermIc(ivector, Herm1, Herm1, half_epsilonH, &tmpMtr2[0][0][0], mask);
		storeHermitDirXYZT(ivector, outHBase[dnDir(dir1)], Herm1, dnDir(dir1), upDir(dir2), mask, 1);
	}
	else {
		sMultAntiHermIc(ivector, &Out1[0][0][0], half_epsilonH, &tmpMtr2[0][0][0], mask);
		storeHermitDirXYZT(ivector, outHBase[dnDir(dir1)], &Out1[0][0][0], dnDir(dir1), dnDir(dir2), mask, 1);
	}
	}
	//endScope(ivector);
	/* Forward */
	//ifStatement(ivector, accumulate + " & " + d2_str.str());
	{
	loadGaugeDirXYZT(ivector, tmpMtr2, inGBase[2*dir1+1], 2*dir1+1, dnDir(dir2), compress12, mask);
	if(dir1==0) loadHermitDirXYZT(ivector, Herm1, outHBase[upDir(dir1)], upDir(dir1), dnDir(dir2), mask);
	//else loadHermitDirXYZT(ivector, Herm1, outHBase[upDir(dir1)], upDir(dir1), dnDir(dir2), mask);
	matMultMatSconj(ivector, tmpMtr, tmpMtr2, Out2, mask);
	matAntiHermTrls(ivector, &tmpMtr2[0][0][0], tmpMtr, mask);
	if(dir1==0) {
		addSMultAntiHermIc(ivector, Herm1, Herm1, half_epsilonH, &tmpMtr2[0][0][0], mask);
		storeHermitDirXYZT(ivector, outHBase[upDir(dir1)], Herm1, upDir(dir1), dnDir(dir2), mask, 1);
	}
	else {
		sMultAntiHermIc(ivector, &Out2[0][0][0], half_epsilonH, &tmpMtr2[0][0][0], mask);
		storeHermitDirXYZT(ivector, outHBase[upDir(dir1)], &Out2[0][0][0], upDir(dir1), dnDir(dir2), mask, 1);
	}
	}
	//endScope(ivector);
      }
      /* otherwise */
      else for(int ud=0; ud<2; ++ud) {
        zeroGResult(ivector, Out1, compress12);
        zeroGResult(ivector, Out2, compress12);
        /**************************/
        /* Calculate gauge forces */
        /**************************/
	/* From staples NOT needed */
	/* From rectangulars & bentchairs */
	stringstream d_str;
        int d = ((1 << (2*dir1+ud)) );
        d_str << d;
	//ifStatement(ivector, accumulate + " & " + d_str.str());
	{
	  g_force_2(ivector, Out1, Out2, kappaR, kappaB, dir2, dir1, (ud==0? UP : DN), compress12);
          /********************/
          /* Update momentum  */
          /* No updates on gauge fields */
          /********************/
	  /* Forward */
	  loadGaugeDirXYZT(ivector, tmpMtr2, inGBase[2*dir1+ud], 2*dir1+ud, upDir(dir2), compress12, mask);
	  if(dir1==0) loadHermitDirXYZT(ivector, Herm1, outHBase[ud==0?dnDir(dir1):upDir(dir1)], ud==0?dnDir(dir1):upDir(dir1), upDir(dir2), mask);
	  //else loadHermitDirXYZT(ivector, Herm1, outHBase[ud==0?dnDir(dir1):upDir(dir1)], ud==0?dnDir(dir1):upDir(dir1), upDir(dir2)-(dir2>dir1), mask);
	  matMultMatSconj(ivector, tmpMtr, tmpMtr2, Out2, mask);
	  matAntiHermTrls(ivector, &tmpMtr2[0][0][0], tmpMtr, mask);
	  if(dir1==0) {
	  	addSMultAntiHermIc(ivector, Herm1, Herm1, half_epsilonH, &tmpMtr2[0][0][0], mask);
	  	storeHermitDirXYZT(ivector, outHBase[ud==0?dnDir(dir1):upDir(dir1)], Herm1, ud==0?dnDir(dir1):upDir(dir1), upDir(dir2), mask, 1);
	  }
	  else {
		sMultAntiHermIc(ivector, &Out2[0][0][0], half_epsilonH, &tmpMtr2[0][0][0], mask);
		storeHermitDirXYZT(ivector, outHBase[ud==0?dnDir(dir1):upDir(dir1)], &Out2[0][0][0], ud==0?dnDir(dir1):upDir(dir1), upDir(dir2)-(dir2>dir1), mask, 1);
	  }
	  /* Backward */
	  loadGaugeDirXYZT(ivector, tmpMtr2, inGBase[2*dir1+ud], 2*dir1+ud, dnDir(dir2), compress12, mask);
	  if(dir1==0) loadHermitDirXYZT(ivector, Herm1, outHBase[ud==0?dnDir(dir1):upDir(dir1)], ud==0?dnDir(dir1):upDir(dir1), dnDir(dir2), mask);
	  //else loadHermitDirXYZT(ivector, Herm1, outHBase[ud==0?dnDir(dir1):upDir(dir1)], ud==0?dnDir(dir1):upDir(dir1), dnDir(dir2)-(dir2>dir1), mask);
	  matMultMatSconj(ivector, tmpMtr, tmpMtr2, Out1, mask);
	  matAntiHermTrls(ivector, &tmpMtr2[0][0][0], tmpMtr, mask);
	  if(dir1==0) {
	  	addSMultAntiHermIc(ivector, Herm1, Herm1, half_epsilonH, &tmpMtr2[0][0][0], mask);
	  	storeHermitDirXYZT(ivector, outHBase[ud==0?dnDir(dir1):upDir(dir1)], Herm1, ud==0?dnDir(dir1):upDir(dir1), dnDir(dir2), mask, 1);
	  }
	  else {
		sMultAntiHermIc(ivector, &Out1[0][0][0], half_epsilonH, &tmpMtr2[0][0][0], mask);
		storeHermitDirXYZT(ivector, outHBase[ud==0?dnDir(dir1):upDir(dir1)], &Out1[0][0][0], ud==0?dnDir(dir1):upDir(dir1), dnDir(dir2)-(dir2>dir1), mask, 1);
	  }
	}
	//endScope(ivector);
      } /* dir1 != dir2 */

    } /* gauge directions loop */

  } /* gauge offsets loop */


}

static string gBase("gBase");
static string hBase("hBase"); 
static string hYZTBase("hYZTBase");
static FVec tmp2[4] = { FVec("tmp20RE"), FVec("tmp20IM"), FVec("tmp21RE"), FVec("tmp21IM") };

/* Update momentum and gauge field body part */
void update_mg_body(InstVector& ivector, bool compress12)
{
  declare_Hermit(ivector, Herm1);
  declare_Hermit(ivector, Herm2);
  declare_Gauge(ivector, tmpMtr);
/*
  for(int i=0; i<4; ++i)
    declareFVecFromFVec(ivector, tmp2[i]);
*/
  //string accumulate("accumulate");
  for(int dir=0; dir<8; ++dir)
  {
    loadHermitDir(ivector, Herm1, hBase, dir, mask);
    loadGaugeDir(ivector, tmpMtr, gBase, dir, compress12, mask);
    for(int dir2=0; dir2<6; ++dir2)
    {
      int d=(1<<(dir2+2));
      stringstream n_str;
      n_str << d;
      //ifStatement(ivector, accumulate + " & " + n_str.str());
      if( dir2+2 != OPPDIR(dir) ) {
        loadHermitYZTDir(ivector, Herm2, hYZTBase, dir2, dir, mask);
        addHerm(ivector, Herm1, Herm1, Herm2, mask);
      }
      //endScope(ivector);
    }
    storeHermitDir(ivector, hBase, Herm1, dir, mask, 1);
    addSMultHermMultMatIc(ivector, tmpMtr, epsilonU, Herm1, tmpMtr, /*tmp2,*/ mask);
    storeGaugeDir(ivector, gBase, tmpMtr, dir, mask, compress12, 1);
  }
}

static FVec tmpMtr7w[7*2][3][2] = { { {FVec("tmpMtr_00RE"), FVec("tmpMtr_00IM")},
                                  {FVec("tmpMtr_01RE"), FVec("tmpMtr_01IM")},
                                  {FVec("tmpMtr_02RE"), FVec("tmpMtr_02IM")} },
                                { {FVec("tmpMtr_10RE"), FVec("tmpMtr_10IM")},
                                  {FVec("tmpMtr_11RE"), FVec("tmpMtr_11IM")},
                                  {FVec("tmpMtr_12RE"), FVec("tmpMtr_12IM")} },
                                { {FVec("tmpMtr_20RE"), FVec("tmpMtr_20IM")},
                                  {FVec("tmpMtr_21RE"), FVec("tmpMtr_21IM")},
                                  {FVec("tmpMtr_22RE"), FVec("tmpMtr_22IM")} },
				{ {FVec("tmpMtr_30RE"), FVec("tmpMtr_30IM")},
                                  {FVec("tmpMtr_31RE"), FVec("tmpMtr_31IM")},
                                  {FVec("tmpMtr_32RE"), FVec("tmpMtr_32IM")} },
                                { {FVec("tmpMtr_40RE"), FVec("tmpMtr_40IM")},
                                  {FVec("tmpMtr_41RE"), FVec("tmpMtr_41IM")},
                                  {FVec("tmpMtr_42RE"), FVec("tmpMtr_42IM")} },
                                { {FVec("tmpMtr_50RE"), FVec("tmpMtr_50IM")},
                                  {FVec("tmpMtr_51RE"), FVec("tmpMtr_1IM")},
                                  {FVec("tmpMtr_52RE"), FVec("tmpMtr_52IM")} },
				{ {FVec("tmpMtr_60RE"), FVec("tmpMtr_60IM")},
                                  {FVec("tmpMtr_61RE"), FVec("tmpMtr_61IM")},
                                  {FVec("tmpMtr_62RE"), FVec("tmpMtr_62IM")} },
                                { {FVec("tmpMtr_70RE"), FVec("tmpMtr_70IM")},
                                  {FVec("tmpMtr_71RE"), FVec("tmpMtr_71IM")},
                                  {FVec("tmpMtr_72RE"), FVec("tmpMtr_72IM")} },
                                { {FVec("tmpMtr_80RE"), FVec("tmpMtr_80IM")},
                                  {FVec("tmpMtr_81RE"), FVec("tmpMtr_81IM")},
                                  {FVec("tmpMtr_82RE"), FVec("tmpMtr_82IM")} }, 
				{ {FVec("tmpMtr_90RE"), FVec("tmpMtr_90IM")},
                                  {FVec("tmpMtr_91RE"), FVec("tmpMtr_91IM")},
                                  {FVec("tmpMtr_92RE"), FVec("tmpMtr_92IM")} },
                                { {FVec("tmpMtr_A0RE"), FVec("tmpMtr_A0IM")},
                                  {FVec("tmpMtr_A1RE"), FVec("tmpMtr_A1IM")},
                                  {FVec("tmpMtr_A2RE"), FVec("tmpMtr_A2IM")} },
                                { {FVec("tmpMtr_B0RE"), FVec("tmpMtr_B0IM")},
                                  {FVec("tmpMtr_B1RE"), FVec("tmpMtr_B1IM")},
                                  {FVec("tmpMtr_B2RE"), FVec("tmpMtr_B2IM")} },
				{ {FVec("tmpMtr_C0RE"), FVec("tmpMtr_C0IM")},
                                  {FVec("tmpMtr_C1RE"), FVec("tmpMtr_C1IM")},
                                  {FVec("tmpMtr_C2RE"), FVec("tmpMtr_C2IM")} },
                                { {FVec("tmpMtr_D0RE"), FVec("tmpMtr_D0IM")},
                                  {FVec("tmpMtr_D1RE"), FVec("tmpMtr_D1IM")},
                                  {FVec("tmpMtr_D2RE"), FVec("tmpMtr_D2IM")} } };

void gf_pack_face(InstVector& ivector, int dir)
{
    std::string l_out("lBuf");
    std::string r_out("rBuf");
    //PrefetchL1GaugeDir(ivector, out, dir, true, 2 /*Exclusive*/);
    /* This will write it to outbufs */
    declare_Gauge7Way(ivector, tmpMtr7w, 1);
    loadGauge12Dir(ivector, tmpMtr7w, "giBase", dir, mask);
    PackGauge7WayDir(ivector, tmpMtr7w, l_out, r_out, dir, 1);
}

/* Receive off-node Gauge fields order : 7 -> 0 */
void gf_unpack_gauge_pack_momentum_face(InstVector& ivector, int Dir, bool compress12)
{
    std::string l_in[8] = { "lBuf[0]", "lBuf[1]", "lBuf[2]", "lBuf[3]", "lBuf[4]", "lBuf[5]", "lBuf[6]", "lBuf[7]" };
    std::string r_in[8] = { "rBuf[0]", "rBuf[1]", "rBuf[2]", "rBuf[3]", "rBuf[4]", "rBuf[5]", "rBuf[6]", "rBuf[7]" };
    std::string h_out[8] = { "hBuf[0]", "hBuf[1]", "hBuf[2]", "hBuf[3]", "hBuf[4]", "hBuf[5]", "hBuf[6]", "hBuf[7]" };
    std::string inGdir[7] = { "inGdir1", "inGdir2", "inGdir3", "inGdir4", "inGdir5", "inGdir6", "inGdir7" };
/*
    std::string l_in0("lBuf[0]"), l_in1("lBuf[1]"), l_in2("lBuf[2]"), l_in3("lBuf[3]"), l_in4("lBuf[4]"), l_in5("lBuf[5]"), l_in6("lBuf[6]"), l_in7("lBuf[7]");
    std::string r_in0("rBuf[0]"), r_in1("rBuf[1]"), r_in2("rBuf[2]"), r_in3("rBuf[3]"), r_in4("rBuf[4]"), r_in5("rBuf[5]"), r_in6("rBuf[6]"), r_in7("rBuf[7]");
    std::string h_out0("hBuf[0]"), h_out1("hBuf[1]"), h_out2("hBuf[2]"), h_out3("hBuf[3]"), h_out4("hBuf[4]"), h_out5("hBuf[5]"), h_out6("hBuf[6]"), h_out7("hBuf[7]");
    std::string inGdir1("inGdir1");
    std::string inGdir2("inGdir2");
    std::string inGdir3("inGdir3");
    std::string inGdir4("inGdir4");
    std::string inGdir5("inGdir5");
    std::string inGdir6("inGdir6");
    std::string inGdir7("inGdir7");
*/
    string face[6] = {"xy", "xz", "xt", "yz", "yt", "zt"};
    string loop[16] = { "00", "01", "02", "03", "04", "05", "06", "07", "10", "11", "12", "13", "14", "15", "16", "17"};
    for(int i=0; i<6; ++i) for(int j=0; j<16; ++j)
    {
        declare_Gauge(ivector, staple1[i][j]);
    }
  declare_Gauge(ivector, tmpMtr);
  declare_Gauge(ivector, tmpMtr2);
  declare_Gauge(ivector, tmpMtr3);
  declare_Gauge(ivector, Out1);
  declare_Gauge(ivector, Out2);
  declare_Hermit7W(ivector, Herm17W); 
  declare_Hermit7W(ivector, Herm27W);   
  string isBoundary("isBoundary");
  declareMask(ivector, "gf_mask");

  /* Load gauges from buffers or local memory */
  for(int dir1=0, i=0; dir1<2*LDim; ++dir1) {
    if(dir1==Dir) UnpackGauge7WayDir(ivector, In[dir1], l_in[dir1], r_in[dir1], dir1, compress12);
    else if(dir1>Dir) {
        ifStatement(ivector, r_in[dir1] + " != 0x00");
        {
	    UnpackGauge7WayDir(ivector, In[dir1], l_in[dir1], r_in[dir1], dir1, compress12);
    	}
    	elseStatement(ivector);
    	{
	    loadGauge7WayDirXYZT(ivector, In[dir1], inGdir[i], dir1, compress12, mask);
    	}
    	endScope(ivector);
    }
    else loadGauge7WayDirXYZT(ivector, In[dir1], inGdir[i], dir1, compress12, mask);
    if(dir1!=Dir) i++;
  }

  /***********/
  /* Staples */
  /***********/
  for(int dir1=0; dir1<LDim; ++dir1) for(int dir2=dir1+1; dir2<LDim; ++dir2) 
  {
    for(int n=0; n<4; ++n)
    {
      if(n==0) {
	staple(ivector, staple1[index(dir1,dir2)][n], /* tmpMtr2 */In[2*dir1+n%2][2*dir2+n/2-1], /* tmpMtr */In[2*dir1+n%2][2*dir1], /* Out2 */In[2*dir2+n/2][2*dir2], UP, mask);
	staple(ivector, staple1[index(dir1,dir2)][n+4], In[2*dir1+n%2][2*dir2+n/2-1], /* Out1 */In[2*dir2+n/2][2*dir1+n%2], In[2*dir2+n/2][2*dir2], DN, mask);
	staple(ivector, staple1[index(dir1,dir2)][n+8], In[2*dir2+n/2][2*dir1+n%2], In[2*dir2+n/2][2*dir2], In[2*dir1+n%2][2*dir1], UP, mask);
	staple(ivector, staple1[index(dir1,dir2)][n+12], In[2*dir2+n/2][2*dir1+n%2], In[2*dir1+n%2][2*dir2+n/2-1], In[2*dir1+n%2][2*dir1], DN, mask);
      }
      else if(n==1)
      {
	staple(ivector, staple1[index(dir1,dir2)][n], In[2*dir2+n/2][2*dir2], In[2*dir1+n%2][2*dir1], In[2*dir1+n%2][2*dir2+n/2-1], UP, mask);
	staple(ivector, staple1[index(dir1,dir2)][n+4], In[2*dir2+n/2][2*dir2], In[2*dir2+n/2][2*dir1+n%2], In[2*dir1+n%2][2*dir2+n/2-1], DN, mask);
	staple(ivector, staple1[index(dir1,dir2)][n+9], In[2*dir2+n/2][2*dir1+n%2], In[2*dir1+n%2][2*dir2+n/2-1], In[2*dir1+n%2][2*dir1], UP, mask);
	staple(ivector, staple1[index(dir1,dir2)][n+13], In[2*dir2+n/2][2*dir1+n%2], In[2*dir2+n/2][2*dir2], In[2*dir1+n%2][2*dir1], DN, mask);
      }
      else if(n==2)
      {
	staple(ivector, staple1[index(dir1,dir2)][n], In[2*dir1+n%2][2*dir2+n/2-1], In[2*dir2+n/2][2*dir1+n%2], In[2*dir2+n/2][2*dir2], UP, mask);
	staple(ivector, staple1[index(dir1,dir2)][n+4], In[2*dir1+n%2][2*dir2+n/2-1], In[2*dir1+n%2][2*dir1], In[2*dir2+n/2][2*dir2], DN, mask);
	staple(ivector, staple1[index(dir1,dir2)][n+7], In[2*dir1+n%2][2*dir1], In[2*dir2+n/2][2*dir2], In[2*dir2+n/2][2*dir1+n%2], UP, mask);
	staple(ivector, staple1[index(dir1,dir2)][n+11], In[2*dir1+n%2][2*dir1], In[2*dir1+n%2][2*dir2+n/2-1], In[2*dir2+n/2][2*dir1+n%2], DN, mask);
      }
      else {
	staple(ivector, staple1[index(dir1,dir2)][n], In[2*dir2+n/2][2*dir2], In[2*dir2+n/2][2*dir1+n%2], In[2*dir1+n%2][2*dir2+n/2-1], UP, mask);
	staple(ivector, staple1[index(dir1,dir2)][n+4], In[2*dir2+n/2][2*dir2], In[2*dir1+n%2][2*dir1], In[2*dir1+n%2][2*dir2+n/2-1], DN, mask);
	staple(ivector, staple1[index(dir1,dir2)][n+8], In[2*dir1+n%2][2*dir1], In[2*dir1+n%2][2*dir2+n/2-1], In[2*dir2+n/2][2*dir1+n%2], UP, mask);
	staple(ivector, staple1[index(dir1,dir2)][n+12], In[2*dir1+n%2][2*dir1], In[2*dir2+n/2][2*dir2], In[2*dir2+n/2][2*dir1+n%2], DN, mask);
      }
    }
  }

  /* Loop over directions of gauge offsets */
  for(int dir1=0; dir1<LDim; ++dir1) { 
    if(dir1==0) {
#if VECLEN==1
	string mvb("0x1"), mvf("0x1");
	string fullMask("0x1");
#elif VECLEN==2
	string mvb("0x1"), mvf("0x2");
	string fullMask("0x3");
#elif VECLEN==4
	string mvb("0x5"), mvf("0xA");
	string fullMask("0xF");
#if VECLEN==8
	string mvb("0x55"), mvf("0xAA");
	string fullMask("0xFF");
#elif VECLEN==16
	string mvb("0x5555"), mvf("0xAAAA");
	string fullMask("0xFFFF");
#endif
	if(dir1>=Dir/2) {
	ifStatement(ivector, h_out[dnDir(dir1)] + " != 0x00");
	{
	IntToMask(ivector, "gf_mask", mvb);
	}
	elseStatement(ivector);
	{
	IntToMask(ivector, "gf_mask", "");
	}
	endScope(ivector);
	loadHermit7WayDirXYZT(ivector, Herm17W, outHBase[dnDir(dir1)], dnDir(dir1), "gf_mask");
	ifStatement(ivector, h_out[upDir(dir1)] + " != 0x00");
	{
	IntToMask(ivector, "gf_mask", mvf);
	}
	elseStatement(ivector);
	{
	IntToMask(ivector, "gf_mask", "");
	}
	endScope(ivector);
	loadHermit7WayDirXYZT(ivector, Herm27W, outHBase[upDir(dir1)], upDir(dir1), "gf_mask");
	}
	else {
	loadHermit7WayDirXYZT(ivector, Herm17W, outHBase[dnDir(dir1)], dnDir(dir1), mask);
	loadHermit7WayDirXYZT(ivector, Herm27W, outHBase[upDir(dir1)], upDir(dir1), mask);
	}
    }

    /* Loop over dirction of gauge links */
    for(int dir2=0; dir2<LDim; ++dir2) {
      /* direction of gauge offset is the same as gauge direction */
      if(dir2==dir1) {
        zeroGResult(ivector, Out1, compress12);
        zeroGResult(ivector, Out2, compress12);
        /**************************/
        /* Calculate gauge forces */
        /**************************/
        /* From staples */
	for(int dir3=0; dir3<LDim; ++dir3) if(dir3!=dir2) {
	    /* Up */
	    addSMulMat(ivector, Out1, Out1, staple1[index(dir2,dir3)][8*(dir3<dir2)+2], kappaS, mask);
	    /* Down */
	    addSMulMat(ivector, Out1, Out1, staple1[index(dir2,dir3)][8*(dir3<dir2)+4], kappaS, mask);
	    /* Up */
	    addSMulMat(ivector, Out2, Out2, staple1[index(dir2,dir3)][8*(dir3<dir2)+3], kappaS, mask);
	    /* Down */
	    addSMulMat(ivector, Out2, Out2, staple1[index(dir2,dir3)][8*(dir3<dir2)+5], kappaS, mask);
	}
	/* From rectangulars */
	g_force_1(ivector, Out1, Out2, kappaR, dir2, compress12, 1, 0);
	/********************/
	/* Update momentum  */
	/********************/
	matMultMatSconj(ivector, tmpMtr, In[2*dir1][upDir(dir2)], Out1, mask);
	matAntiHermTrls(ivector, &tmpMtr2[0][0][0], tmpMtr, mask);
	if(dir1==0) {
	    addSMultAntiHermIc(ivector, Herm17W[dnDir(dir2)], Herm17W[dnDir(dir2)], half_epsilonH, &tmpMtr2[0][0][0], mask);
        }
        else {
	    sMultAntiHermIc(ivector, Herm17W[dnDir(dir2)], half_epsilonH, &tmpMtr2[0][0][0], mask);
        }
	matMultMatSconj(ivector, tmpMtr, In[2*dir1+1][dnDir(dir2)], Out2, mask);
	matAntiHermTrls(ivector, &tmpMtr2[0][0][0], tmpMtr, mask);
	if(dir1==0) 
	    addSMultAntiHermIc(ivector, Herm27W[dnDir(dir2)], Herm27W[dnDir(dir2)], half_epsilonH, &tmpMtr2[0][0][0], mask);
	else
	    sMultAntiHermIc(ivector, Herm27W[dnDir(dir2)], half_epsilonH, &tmpMtr2[0][0][0], mask);
      }

      /* otherwise */
      else {
	FVec Herm7W[2][7][9] = {Herm17W, Herm27W};
	for(int ud=0; ud<2; ++ud) {
	    zeroGResult(ivector, Out1, compress12);
	    zeroGResult(ivector, Out2, compress12);
	    /**************************/
	    /* Calculate gauge forces */
	    /**************************/
	    g_force_2(ivector, Out1, Out2, kappaR, kappaB, dir2, dir1, (ud==0? UP : DN), compress12, 1, 0);

	    /********************/
	    /* Update momentum  */
	    /********************/
	    /* Forward */
	    matMultMatSconj(ivector, tmpMtr, In[2*dir1+ud][upDir(dir2)-(dir2>dir1)], Out2, mask);
	    matAntiHermTrls(ivector, &tmpMtr2[0][0][0], tmpMtr, mask);
	    if(dir1==0) 
	    	addSMultAntiHermIc(ivector, Herm7W[ud][upDir(dir2)-1], Herm7W[ud][upDir(dir2)-1], half_epsilonH, &tmpMtr2[0][0][0], mask);
	    else 
	    	sMultAntiHermIc(ivector, Herm7W[ud][upDir(dir2)-(dir2>dir1)], half_epsilonH, &tmpMtr2[0][0][0], mask);
	    /* Backward */
 	    matMultMatSconj(ivector, tmpMtr, In[2*dir1+ud][dnDir(dir2)-(dir2>dir1)], Out1, mask);
	    matAntiHermTrls(ivector, &tmpMtr2[0][0][0], tmpMtr, mask);
	    if(dir1==0)
	    	addSMultAntiHermIc(ivector, Herm7W[ud][dnDir(dir2)-1], Herm7W[ud][dnDir(dir2)-1], half_epsilonH, &tmpMtr2[0][0][0], mask);
	    else
	    	sMultAntiHermIc(ivector, Herm7W[ud][dnDir(dir2)-(dir2>dir1)], half_epsilonH, &tmpMtr2[0][0][0], mask);
        }
      }
    }

#if VECLEN==1
        string mvb[4] = {"0x1", "0x1", "0x1", "0x1"};
	string mvf[4] = mvb;
        string fullMask("0x1");
#elif VECLEN==2
        string mvf[4] = {"0x1", "0x3", "0x3", "0x3"};
	string mvb[4] = {"0x2", "0x3", "0x3", "0x3"};
#elif VECLEN==4
        string mvf[4] = {"0x5", "0x3", "0xF", "0xF"};
	string mvb[4] = {"0xA", "0xC", "0xF", "0xF"};
#if VECLEN==8
        string mvf[4] = {"0x55", "0x33", "0x0F", "0xFF"};
	string mvb[4] = {"0xAA", "0xCC", "0xF0", "0xFF"};
#elif VECLEN==16
        string mvf[4] = {"0x5555", "0x3333", "0x0F0F", "0x00FF"};
	string mvb[4] = {"0xAAAA", "0xCCCC", "0xF0F0", "0xFF00"};
#endif
    if(2*dir1>Dir) {
    ifStatement(ivector, h_out[dnDir(dir1)] + " != 0x00");
    {
    PackHermit7WayDir(ivector, Herm17W, h_out[dnDir(dir1)], mask);
    IntToMask(ivector, "gf_mask", mvb[dir1]);
    }
    elseStatement(ivector);
    {
    IntToMask(ivector, "gf_mask", "");
    }
    endScope(ivector);
    storeHermit7WayDirXYZT(ivector, outHBase[dnDir(dir1)], Herm17W, dnDir(dir1), "gf_mask", 1);
    }
    else if(2*dir1==Dir) {
	IntToMask(ivector, "gf_mask", mvb[dir1]);
	storeHermit7WayDirXYZT(ivector, outHBase[dnDir(dir1)], Herm17W, dnDir(dir1), "gf_mask", 1);
    }
    else storeHermit7WayDirXYZT(ivector, outHBase[dnDir(dir1)], Herm17W, dnDir(dir1), mask, 1);
    if(2*dir1+1>Dir) {
    ifStatement(ivector, h_out[upDir(dir1)] + " != 0x00");
    {
    PackHermit7WayDir(ivector, Herm27W, h_out[upDir(dir1)], mask);
    IntToMask(ivector, "gf_mask", mvf[dir1]);
    }
    elseStatement(ivector);
    {
    IntToMask(ivector, "gf_mask", "");
    }
    endScope(ivector);
    storeHermit7WayDirXYZT(ivector, outHBase[upDir(dir1)], Herm27W, upDir(dir1), "gf_mask", 1);
    endScope(ivector);
    }
    else if(2*dir1+1==Dir) {
	IntToMask(ivector, "gf_mask", mvf[dir1]);
	storeHermit7WayDirXYZT(ivector, outHBase[upDir(dir1)], Herm27W, upDir(dir1), "gf_mask", 1);
    }
    else storeHermit7WayDirXYZT(ivector, outHBase[upDir(dir1)], Herm27W, upDir(dir1), mask, 1);
  }

}

void generateGFImpL2Prefetches(InstVector& ivector, bool compress12)
{
  for(int i=0; i<4; ++i)
    PrefetchL2FullGaugeIn(ivector, prefGBD[i][0], prefGBD[i][1], compress12);
/*
  for(int i=0; i<4; ++i)
    PrefetchL2FullHermitIn(ivector, prefHBD[i][0], prefHBD[i][1]);
*/
}

void generateGFImpAddL2Prefetches(InstVector& ivector, bool compress12)
{
  PrefetchL2FullHermitIn(ivector, inHBase, iprefHBD);
  for(int i=0; i<4; ++i)
    PrefetchL2FullHermitHelperIn(ivector, prefHBD[i][0], prefHBD[i][1]);
}

int main(void)
{
    InstVector ivector;
    InstVector l2prefs;
    bool compress12;

    const string HermitTypeName = getTypeName(sizeof(HermitBaseType));
    const string GaugeTypeName = getTypeName(sizeof(GaugeBaseType));

#if 0
#ifdef NO_MASKS
    for(int i = 0; i < 4; i++) requireAllOneCheck[i] = false;
#endif
#endif

    for(int num_components=12; num_components <=18; num_components+=6) {
        compress12 = ( num_components==12 );

	/* GF body part calculations */
        std::ostringstream filename;
        filename << "./"<<ARCH_NAME<<"/" << "gf_imp_body_" << HermitTypeName << "_" << GaugeTypeName << "_v"<< VECLEN <<"_s"<<"_"<< num_components;
        cout << "GENERATING " << filename.str() << endl;
        l2prefs.resize(0);
        ivector.resize(0);
        generateGFImpL2Prefetches(l2prefs,compress12);
        g_force_imp_body(ivector,compress12);
        mergeIvectorWithL2Prefetches(ivector, l2prefs);
        dumpIVector(ivector,filename.str());
        filename.str("");
        filename.clear();
#if 0
        filename << "./"<<ARCH_NAME<<"/" << "gf_imp_body_add_" << HermitTypeName << "_" << GaugeTypeName << "_v"<< VECLEN <<"_s"<<"_"<< num_components;
        cout << "GENERATING " << filename.str() << endl;
        l2prefs.resize(0);
        ivector.resize(0);
        generateGFImpAddL2Prefetches(l2prefs,compress12);
        g_force_imp_body_add(ivector,compress12);
        mergeIvectorWithL2Prefetches(ivector, l2prefs);
        dumpIVector(ivector,filename.str());
        filename.str("");
        filename.clear();
#endif

	filename << "./"<<ARCH_NAME<<"/" << "update_mg_body_" << HermitTypeName << "_" << GaugeTypeName << "_v"<< VECLEN <<"_s"<<"_"<< num_components;
        cout << "GENERATING " << filename.str() << endl;
        l2prefs.resize(0);
        ivector.resize(0);
	update_mg_body(ivector,compress12);
	mergeIvectorWithL2Prefetches(ivector, l2prefs);
	dumpIVector(ivector,filename.str());
	filename.str("");
	filename.clear();

    }

    std::ostringstream bf[2];
    bf[0] << "back";
    bf[1] << "forw";
    char Dir[4] = { 'X', 'Y', 'Z', 'T' };

    /* Pack gauge fields */
    for(int dir=0; dir<8; ++dir)
    {
    	std::ostringstream filename;
    	filename << "./"<<ARCH_NAME<<"/" << "gauge_pack_to_" << bf[dir&1].str() << "_" << Dir[dir/2] << "_" << GaugeTypeName << "_v"<< VECLEN <<"_s"<<"_12";
    	cout << "GENERATING " << filename.str() << endl;
    	l2prefs.resize(0);
    	ivector.resize(0);
    	gf_pack_face(ivector, int dir);
    	mergeIvectorWithL2Prefetches(ivector, l2prefs);
    	dumpIVector(ivector,filename.str());
    	filename.str("");
    	filename.clear();
    }

}
