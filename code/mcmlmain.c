#pragma once
/***********************************************************
 *  Copyright Univ. of Texas M.D. Anderson Cancer Center
 *  1992.
 *
 *	main program for Monte Carlo simulation of photon
 *	distribution in multi-layered turbid media.
 *
 ****/

 /****
  *	THINKCPROFILER is defined to generate profiler calls in
  *	Think C. If 1, remember to turn on "Generate profiler
  *	calls" in the options menu.
  ****/
#define THINKCPROFILER 0	

  /* GNU cc does not support difftime() and CLOCKS_PER_SEC.*/
#define GNUCC 0

#if THINKCPROFILER
#include <profile.h>
#include <console.h>
#endif

#include "mcml.h"

/*	Declare before they are used in main(). */
FILE* GetFile(char*);
short ReadNumRuns(FILE*);
void ReadParm(FILE*, InputStruct*);
void CheckParm(FILE*, InputStruct*);
void InitOutputData(InputStruct, OutStruct*);
void InitEmit(circle_ini*);
void FreeData(InputStruct, OutStruct*,circle_ini*);
double Rspecular(LayerStruct*);
void LaunchPhoton(double ,LayerStruct* , PhotonStruct* ,
    circle_ini* , InputStruct* );
void HopDropSpin(InputStruct*, PhotonStruct*, OutStruct*);
void SumScaleResult(InputStruct, OutStruct*);
void WriteResult(InputStruct, OutStruct, char*);
void WriteMatrix(InputStruct, OutStruct);
void ini_xy(circle_ini*, InputStruct, PhotonStruct*);

int raman_count = 0;
int axis_time = 0;
int matrix = 1;
double radius;
int mode_shift = 0;
int photon_count=0;
int nums_photon_global=0;
double laser_r = 0.2;
int version;
double depth;
char name[64] = "";





/***********************************************************
 *	If F = 0, reset the clock and return 0.
 *
 *	If F = 1, pass the user time to Msg and print Msg on
 *	screen, return the real time since F=0.
 *
 *	If F = 2, same as F=1 except no printing.
 *
 *	Note that clock() and time() return user time and real
 *	time respectively.
 *	User time is whatever the system allocates to the
 *	running of the program;
 *	real time is wall-clock time.  In a time-shared system,
 *	they need not be the same.
 *
 *	clock() only hold 16 bit integer, which is about 32768
 *	clock ticks.
 ****/
time_t PunchTime(char F, char* Msg)
{
#if GNUCC
    return(0);
#else
    static clock_t ut0;	/* user time reference. */
    static time_t  rt0;	/* real time reference. */
    double secs;
    char s[STRLEN];

    if (F == 0) {
        ut0 = clock();
        rt0 = time(NULL);
        return(0);
    }
    else if (F == 1) {
        secs = (clock() - ut0) / (double)CLOCKS_PER_SEC;
        if (secs < 0) secs = 0;	/* clock() can overflow. */
        sprintf(s, "User time: %8.0lf sec = %8.2lf hr.  %s\n",
            secs, secs / 3600.0, Msg);
        puts(s);
        strcpy(Msg, s);
        return(difftime(time(NULL), rt0));
    }
    else if (F == 2) return(difftime(time(NULL), rt0));
    else return(0);
#endif
}

/***********************************************************
 *	Print the current time and the estimated finishing time.
 *
 *	P1 is the number of computed photon packets.
 *	Pt is the total number of photon packets.
 ****/
void PredictDoneTime(long P1, long Pt)
{
    time_t now, done_time;
    struct tm* date;
    char s[80];

    now = time(NULL);
    date = localtime(&now);
    strftime(s, 80, "%H:%M %x", date);
    printf("Now %s, ", s);

    done_time = now +
        (time_t)(PunchTime(2, "") / (double)P1 * (Pt - P1));
    date = localtime(&done_time);
    strftime(s, 80, "%H:%M %x", date);
    printf("End %s\n", s);
}

/***********************************************************
 *	Report time and write results.
 ****/
void ReportResult(InputStruct In_Parm, OutStruct Out_Parm)
{
    char time_report[STRLEN];

    strcpy(time_report, " Simulation time of this run.");
    PunchTime(1, time_report);

    SumScaleResult(In_Parm, &Out_Parm);
    if (matrix == 1)
    {
        WriteMatrix(In_Parm, Out_Parm);
    }
    printf("~~~~~~~~~~~~\n");
    WriteResult(In_Parm, Out_Parm, time_report);
    printf("~~~~~~~~~~~~");
}

/***********************************************************
 *	Get the file name of the input data file from the
 *	argument to the command line.
 ****/
void GetFnameFromArgv(int argc,
    char* argv[],
    char* input_filename)
{
    if (argc >= 2) {			/* filename in command line */
        strcpy(input_filename, argv[1]);
    }
    else
        input_filename[0] = '\0';
}


/***********************************************************
 *	Execute Monte Carlo simulation for one independent run.
 ****/
void DoOneRun(short NumRuns, InputStruct* In_Ptr)
{
    register long i_photon;
    /* index to photon. register for speed.*/
    OutStruct out_parm;		/* distribution of photons.*/
    PhotonStruct photon;
    circle_ini ini_ptr;
    long num_photons = In_Ptr->num_photons, photon_rep = 10;
    nums_photon_global = In_Ptr->num_photons;
    

#if THINKCPROFILER
    InitProfile(200, 200); cecho2file("prof.rpt", 0, stdout);
#endif
    InitOutputData(*In_Ptr, &out_parm);
    InitEmit(&ini_ptr);
    ini_xy(&ini_ptr, *In_Ptr, &photon);
    printf("\n******************************初始化完成*********************************");
    out_parm.Rsp = Rspecular(In_Ptr->layerspecs);
    i_photon = num_photons;

    PunchTime(0, "");

    do {
        printf("******NO.:%d******", i_photon);
        LaunchPhoton(out_parm.Rsp, In_Ptr->layerspecs, &photon, &ini_ptr, In_Ptr->num_layers, In_Ptr);
        printf("\n******************************坐标输出完成*********************************\n");
        do  HopDropSpin(In_Ptr, &photon, &out_parm);
        while (!photon.dead==1);
    } while (--i_photon);

#if THINKCPROFILER
    exit(0);
#endif

    ReportResult(*In_Ptr, out_parm);
    FreeData(*In_Ptr, &out_parm,&ini_ptr);
}



void Write_raman_givexy(double a,double b)
{


    FILE* file;
    char        fname[STRLEN];
    
    sprintf(fname, "ini_xy-%s.txt",name);
     
    
    file = fopen(fname, "a+");
    if (file == NULL)
        return;

    fprintf(file, "%-12.4E\t%-12.4E\n", a, b);

    
    /*
    for (int i = 0; i < nums_photon_global; i++)
    {
        xs=ini_ptr->x_r[i];
        ys=ini_ptr->y_r[i];  
        
    }*/
    fclose(file);
}
void ini_xy(circle_ini* ini_ptr, InputStruct In_Parm,PhotonStruct* Photon_Ptr)
{
    
    int requeire;//光子模式选择
    ini_ptr->give_count = 0;
    printf("\n******************************请输入你想要的******************************\n");
    printf("#1:随机均匀分布圆\n");
    printf("#2:MC-F立体高斯\n");
    printf("#3:平面高斯圆（r）\n");
    printf("#4:原始（zr）\n");
    printf("#5:特定layer\n");
    printf("你快输入:\n");
    scanf_s("%d", &requeire);
    switch (requeire)
    {
    case(1):
    {
        version = 1;
        double rr;
        printf("\n******************************请输入半径******************************\n");
        scanf_s("%lf", &rr);
        rr = rr / 10000;
        printf("你输入的r是：%f", rr);
        for (int i = 0;i < nums_photon_global;i++)
        {
            double u = sqrt(RAND) * rr;
            double v = RAND * 2 * PI;
            ini_ptr->x_r[i] = u * cos(v);
            ini_ptr->y_r[i] = u * sin(v);
        } 
        //Write_raman_givexy(ini_ptr->x_r, ini_ptr->y_r);
        break;
    }
    case(2):
    {
        version = 2;
        printf("直接pass后面在初始化\n");
        /*
        scanf_s("%lf", &shift_r);
        printf("平移几次\n");
        scanf_s("%d", &shift_role);
        laser_num= nums_photon_global / shift_role;
        shift_r = shift_r / 10000;
        for (int i = 0; i < shift_role; i++)
        {
            for (int k = 0; k < laser_num; k++)
            {
                ini_ptr->x_r[flag_end+k] = i * shift_r;
                if(k==laser_num-1 )
                {
                    flag_end = (i+1)*laser_num-1;
                }
            }
        }*/
        break;
    }

    case(3):
    {
        version = 3;
        //double gau_mean = 0;
        //double gau_variance = 0;
        char input_xy[512];
        char buf[1024];

        FILE* XY;
        int c = 0;
        printf("输入配置文件：");
        XY = GetFile(input_xy);
        while (fgets(buf, 1024, XY) != NULL)
        {
            sscanf(buf, "%lf%lf", &(ini_ptr->x_r[c]), &(ini_ptr->y_r[c]));
            c++;
        }
        //Write_raman_givexy(ini_ptr->x_r, ini_ptr->y_r);
        fclose(XY);
        break;

    }

    case(4):
    {
        printf("raw");
        version = 4;
        break;
    }
    case(5):
    {
        version = 5;
        double l_x;
        double l_y;
        double l_z;
        printf("\n******************************请输入X******************************\n");
        scanf_s("%lf", &l_x);
        printf("\n******************************请输入Y******************************\n");
        scanf_s("%lf", &l_y);

        char input_xy[1024];
        char buf[1024];

        FILE* XY;
        int c = 0;
        printf("输入配置文件：");
        XY = GetFile(input_xy);
        while (fgets(buf, 1024, XY) != NULL)
        {
            sscanf(buf, "%lf%lf", &(ini_ptr->x_r[c]), &(ini_ptr->y_r[c]));
            c++;
        }
        //Write_raman_givexy(ini_ptr->x_r, ini_ptr->y_r);
        fclose(XY);
        
        
        printf("\n******************************请输入Z******************************\n");
        //Photon_Ptr->z = 0.0;
        scanf_s("%lf", &Photon_Ptr->z);
        for (int j = 1; j <= In_Parm.num_layers; j++)
        {
            if (Photon_Ptr->z >= In_Parm.layerspecs[j].z0 && Photon_Ptr->z < In_Parm.layerspecs[j].z1)
            {
                Photon_Ptr->layer = j;
            }
        }
        printf("\n******************************选择方向******************************\n");
        //Photon_Ptr->z = 0.0;
        scanf_s("%lf", &Photon_Ptr->uz);
        break;
    }
    /*
    for (int i = 0; i < (int)sqrt_n; i++)
    {
        a = i * rr / sqrt_n;
        b = i * (360 / sqrt_n);
        r_ini[i] = a;
        theta_ini[i] = b;
    }
    if (res != 0)
    {
       for (int resss = 0; resss < res; resss++)
        {
           ini_ptr->x_r[resss] = 0;
           ini_ptr->y_r[resss] = 0;
        }
       k = res;
    }

    for (int i = 0; i < (int)sqrt_n; i++)
    {
        for (int j = 0; j < (int)sqrt_n; j++)
        {
            c = r_ini[i] * cos(theta_ini[j] * PI / 180);
            ini_ptr->x_r[k]=c;
            d = r_ini[i] * sin(theta_ini[j] * PI / 180);
            ini_ptr->y_r[k]=d;
            k++;
        }

    }*/
    }
   

}
/***********************************************************
 *	The argument to the command line is filename, if any.
 *	Macintosh does not support command line.
 ****/

char main(int argc, char* argv[])
{
    do{
        char input_filename[STRLEN];
        FILE* input_file_ptr;
        short num_runs;	/* number of independent runs. */
        InputStruct in_parm;
        ShowVersion("Version 1.2, 1993");
        GetFnameFromArgv(argc, argv, input_filename);
        input_file_ptr = GetFile(input_filename);
        CheckParm(input_file_ptr, &in_parm);
        num_runs = ReadNumRuns(input_file_ptr);
    
        while (num_runs--) {
            ReadParm(input_file_ptr, &in_parm);
            printf("\n输入文件名\n");
            scanf("%s", &name);
            DoOneRun(num_runs, &in_parm);
            printf("\n************************总共有%d拉曼光子************************\n", raman_count);
        }
    
        fclose(input_file_ptr);
        getchar();
        system("pause");
    } while (1);
    return 0;
   
}