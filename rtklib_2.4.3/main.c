#include "src\rtklib.h"

#define TRACEFILE   ""

void main()
{
	traceopen(TRACEFILE);
	tracelevel(3);

	gtime_t ts={0},te={0};
	prcopt_t prcopt=prcopt_default;
	solopt_t solopt=solopt_default;
	filopt_t filopt={0};
	char* infile[]={{"C:\\Users\\LW\\Documents\\rtklib\\vs0\\rtklib_2.4.3\\rtklib_2.4.3\\data\\1.obs"},
					{"C:\\Users\\LW\\Documents\\rtklib\\vs0\\rtklib_2.4.3\\rtklib_2.4.3\\data\\2.obs"},
					{"C:\\Users\\LW\\Documents\\rtklib\\vs0\\rtklib_2.4.3\\rtklib_2.4.3\\data\\1.nav"},
	};
	char* outfile="C:\\Users\\LW\\Documents\\rtklib\\vs0\\rtklib_2.4.3\\rtklib_2.4.3\\data\\pos.pos";

	prcopt.mode= PMODE_STATIC;
	prcopt.navsys = SYS_GPS;
	solopt.posf=SOLF_ENU;
	postpos(ts,te,0,0,&prcopt,&solopt,&filopt,infile,3,outfile,NULL,NULL);
	
	/*double opos[3], rr[3],pos[3], r[3], enu[3];
	rr[0] = -2191128.608;
	rr[1] = 5198843.810;
	r[2] = 2965371.956;
	for (int i = 0; i < 3; i++) r[i] = rr[i] - opos[i];
	ecef2pos(opos, pos);
	ecef2enu(pos, r, enu);*/
}