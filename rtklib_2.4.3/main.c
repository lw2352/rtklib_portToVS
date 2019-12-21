#include "src\rtklib.h"


void main1()
{

	gtime_t ts={0},te={0};
	prcopt_t prcopt=prcopt_default;
	solopt_t solopt=solopt_default;
	filopt_t filopt={0};
	char* infile[]={{"D:\\rtklib\\testdata\\07590920.05o"},
					{"D:\\rtklib\\testdata\\30400920.05o"},
					{"D:\\rtklib\\testdata\\07590920.05n"},
	};
	char* outfile="D:\\rtklib\\testdata\\pos.pos";

	prcopt.mode=PMODE_MOVEB;
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