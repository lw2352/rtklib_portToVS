#include "src/rtklib.h"

#define TRACEFILE   "../rtklib_2.4.3/trace_%Y%m%d%h%M.txt"
#define OPTSDIR     "."                 /* default config directory */
#define OPTSFILE    "../rtklib_2.4.3/rtk.conf"       /* default config file */
#define MAXSTR      1024                /* max length of a stream */

static char passwd[MAXSTR] = "admin";     /* login password */
static int timetype = 0;             /* time format (0:gpst,1:utc,2:jst,3:tow) */
static int soltype = 0;             /* sol format (0:dms,1:deg,2:xyz,3:enu,4:pyl) */
static int solflag = 2;             /* sol flag (1:std+2:age/ratio/ns) */
static int svrcycle = 10;            /* server cycle (ms) */
static int timeout = 10000;         /* timeout time (ms) */
static int reconnect = 10000;         /* reconnect interval (ms) */
static int nmeacycle = 5000;          /* nmea request cycle (ms) */
static int buffsize = 32768;         /* input buffer size (bytes) */
static int navmsgsel = 0;             /* navigation mesaage select */
static char proxyaddr[256] = "";          /* http/ntrip proxy */
static int nmeareq = 0;             /* nmea request type (0:off,1:lat/lon,2:single) */
static double nmeapos[] = { 0,0,0 };       /* nmea position (lat/lon/height) (deg,m) */
static char rcvcmds[3][MAXSTR] = { "" };    /* receiver commands files */
static char startcmd[MAXSTR] = "";        /* start command */
static char stopcmd[MAXSTR] = "";        /* stop command */
static int modflgr[256] = { 0 };           /* modified flags of receiver options */
static int modflgs[256] = { 0 };           /* modified flags of system options */
static int moniport = 0;             /* monitor port */
static int keepalive = 0;             /* keep alive flag */
static int fswapmargin = 30;            /* file swap margin (s) */
static char sta_name[256] = "";           /* station name */

static int strtype[] = {                  /* stream types */
    STR_SERIAL,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE
};
static char strpath[8][MAXSTR] = { "","","","","","","","" }; /* stream paths */
static int strfmt[] = {                   /* stream formats */
    STRFMT_UBX,STRFMT_RTCM3,STRFMT_SP3,SOLF_LLH,SOLF_NMEA
};
/* receiver options table ----------------------------------------------------*/
#define TIMOPT  "0:gpst,1:utc,2:jst,3:tow"
#define CONOPT  "0:dms,1:deg,2:xyz,3:enu,4:pyl"
#define FLGOPT  "0:off,1:std+2:age/ratio/ns"
#define ISTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,7:ntripcli,8:ftp,9:http"
#define OSTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,6:ntripsvr,11:ntripc_c"
#define FMTOPT  "0:rtcm2,1:rtcm3,2:oem4,3:oem3,4:ubx,5:ss2,6:hemis,7:skytraq,8:gw10,9:javad,10:nvs,11:binex,12:rt17,13:sbf,14:cmr,15:tersus,18:sp3"
#define NMEOPT  "0:off,1:latlon,2:single"
#define SOLOPT  "0:llh,1:xyz,2:enu,3:nmea,4:stat"
#define MSGOPT  "0:all,1:rover,2:base,3:corr"

static opt_t rcvopts[] = {
    {"console-passwd",  2,  (void*)passwd,              ""     },
    {"console-timetype",3,  (void*)&timetype,           TIMOPT },
    {"console-soltype", 3,  (void*)&soltype,            CONOPT },
    {"console-solflag", 0,  (void*)&solflag,            FLGOPT },

    {"inpstr1-type",    3,  (void*)&strtype[0],         ISTOPT },
    {"inpstr2-type",    3,  (void*)&strtype[1],         ISTOPT },
    {"inpstr3-type",    3,  (void*)&strtype[2],         ISTOPT },
    {"inpstr1-path",    2,  (void*)strpath[0],         ""     },
    {"inpstr2-path",    2,  (void*)strpath[1],         ""     },
    {"inpstr3-path",    2,  (void*)strpath[2],         ""     },
    {"inpstr1-format",  3,  (void*)&strfmt[0],         FMTOPT },
    {"inpstr2-format",  3,  (void*)&strfmt[1],         FMTOPT },
    {"inpstr3-format",  3,  (void*)&strfmt[2],         FMTOPT },
    {"inpstr2-nmeareq", 3,  (void*)&nmeareq,            NMEOPT },
    {"inpstr2-nmealat", 1,  (void*)&nmeapos[0],         "deg"  },
    {"inpstr2-nmealon", 1,  (void*)&nmeapos[1],         "deg"  },
    {"inpstr2-nmeahgt", 1,  (void*)&nmeapos[2],         "m"    },
    {"outstr1-type",    3,  (void*)&strtype[3],         OSTOPT },
    {"outstr2-type",    3,  (void*)&strtype[4],         OSTOPT },
    {"outstr1-path",    2,  (void*)strpath[3],         ""     },
    {"outstr2-path",    2,  (void*)strpath[4],         ""     },
    {"outstr1-format",  3,  (void*)&strfmt[3],         SOLOPT },
    {"outstr2-format",  3,  (void*)&strfmt[4],         SOLOPT },
    {"logstr1-type",    3,  (void*)&strtype[5],         OSTOPT },
    {"logstr2-type",    3,  (void*)&strtype[6],         OSTOPT },
    {"logstr3-type",    3,  (void*)&strtype[7],         OSTOPT },
    {"logstr1-path",    2,  (void*)strpath[5],         ""     },
    {"logstr2-path",    2,  (void*)strpath[6],         ""     },
    {"logstr3-path",    2,  (void*)strpath[7],         ""     },

    {"misc-svrcycle",   0,  (void*)&svrcycle,           "ms"   },
    {"misc-timeout",    0,  (void*)&timeout,            "ms"   },
    {"misc-reconnect",  0,  (void*)&reconnect,          "ms"   },
    {"misc-nmeacycle",  0,  (void*)&nmeacycle,          "ms"   },
    {"misc-buffsize",   0,  (void*)&buffsize,           "bytes"},
    {"misc-navmsgsel",  3,  (void*)&navmsgsel,          MSGOPT },
    {"misc-proxyaddr",  2,  (void*)proxyaddr,           ""     },
    {"misc-fswapmargin",0,  (void*)&fswapmargin,        "s"    },

    {"misc-startcmd",   2,  (void*)startcmd,            ""     },
    {"misc-stopcmd",    2,  (void*)stopcmd,             ""     },

    {"file-cmdfile1",   2,  (void*)rcvcmds[0],          ""     },
    {"file-cmdfile2",   2,  (void*)rcvcmds[1],          ""     },
    {"file-cmdfile3",   2,  (void*)rcvcmds[2],          ""     },

    {"",0,NULL,""}
};

void main()
{
    int level = 1;
#if 1
    //test();
	rtkrcv(level);
#elif 0
    //my test
    test();
    return 0;
#else
	traceopen(TRACEFILE);
	tracelevel(level);
	char file[MAXSTR] = "";
	gtime_t ts={0},te={0};
    double es[] = { 2021,6,8,11,22,00 }, ee[] = { 2021,6,8,11,26,30 };
    //ts = epoch2time(es);te = epoch2time(ee);
    
	prcopt_t prcopt=prcopt_default;
	solopt_t solopt=solopt_default;
	filopt_t filopt={""};
	char* infile[]={
        {"D:\\data\\6-8\\rover_202106081120.obs"},
        {"D:\\data\\6-8\\base_202106081120.obs"},
	    {"D:\\data\\6-8\\rover_202106081120.nav"}
	};
	char* outfile="D:\\data\\6-8\\sol.txt";
	
	//prcopt.mode= PMODE_STATIC;
	//prcopt.navsys = SYS_GPS;
	//solopt.posf=SOLF_ENU;

	/* load options file */
	if (!*file) sprintf(file, "%s/%s", OPTSDIR, OPTSFILE);

	resetsysopts();
	int a = loadopts(file, rcvopts);
	int b = loadopts(file, sysopts);

	getsysopts(&prcopt, &solopt, &filopt);

	postpos(ts,te,0,0,&prcopt,&solopt,&filopt,infile,3,outfile,NULL,NULL);
#endif
    return 0;
}