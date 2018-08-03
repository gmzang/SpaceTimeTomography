#include "GLUTSceneContainer.h" // 
#include "ST_Tomography.h"
#include "VolumeData.h" 
#include <argtable2.h>

extern "C" { FILE __iob_func[3] = { *stdin,*stdout,*stderr }; }

int main(int argc,char** argv)
{



	struct arg_dbl *d_sigma=arg_dbl0("s","Sigma",        "<double>"     ,     "Sigma in volume density reconstruction (default as 0.2)");
	struct arg_dbl *d_tau=arg_dbl0("t","Tau",  "<double>"       ,         "Lambda in volume density reconstruction (default as 0.2)");
	struct arg_dbl *d_lambda=arg_dbl0("l","Lambda",  "<double>"       ,         "TV prior weight in volume density reconstruction ");
	struct arg_dbl *d_tp = arg_dbl0("k", "TemporalPrior", "<double>", "temporal prior weight in volume density reconstruction  ");
	struct arg_dbl *d_bc = arg_dbl0("e", "Bright constancy prior", "<double>", "Bright constancy weight in volume density reconstruction ");
	struct arg_dbl *d_huber = arg_dbl0(NULL, "huber", "<double>", "trade-off weight for huber norm epsilon ");
	struct arg_dbl *alp=arg_dbl0("a","Alpha",  "<double>"       ,         "stepsize for SART algorithm  (default as 0.3) ");
	
	struct arg_file *outfile = arg_file0("o",NULL,"<output>",           "output 3d image (default as \"hello.tif\")");

	struct arg_int *nframe = arg_int0("v", "nframes", "<int>", "Number of frames");
	struct arg_int *nrounds = arg_int0("r", "nframes", "<int>", "Number of rounds for each frames");

	struct arg_int  *vols = arg_intn("f", "XYZ", NULL, 1, 3, "takes an integer value (defaults to 9)");

	struct arg_int *imgWH = arg_intn("u", "imgUV", NULL,1,2, "w and h of projection image default:512 ");


	struct arg_dbl *vs=arg_dbl0("c","voxelspacing",  "<double>"       ,         "voxel spacing (default as 1)");
	struct arg_dbl *ds=arg_dbl0("d","detector spacing",  "<double>"       ,         "detector spacing (default as 1)");
	struct arg_file *projfiles = arg_filen("i","projsvolume", NULL,1,20, "Input Dataset file(.mha\.tif\.tiff\ 2d bmp and other raw file are supported)");

	struct arg_dbl *sdd=arg_dbl0("j","sdd",  "<double>"       ,         " Source to Detector Distance (default as 1000mm)");
	struct arg_dbl *sid=arg_dbl0("g","sid",  "<double>"       ,         " Source to Iso-Object Distance  (default as 600mm)");
	struct arg_dbl *offsetXYZ = arg_dbln(NULL, "oxyz", NULL, 1, 3, " offset x y z of the center of volume (default as 0)");
	struct arg_dbl *startdegree = arg_dbln(NULL, "sdg", NULL, 1, 20, " start degree for each proj sequence");
	struct arg_int *AlgoIter=arg_int0("b","AlgorithmIter",  "<int>"       ,         "Number of algorithm iterations, default as 20");
	struct arg_int *SartIter=arg_int0("p","SartIter",  "<int>"       ,         "Number of SART nested iterations, default as 1");
	// add for prior options
	struct arg_str  *prior    = arg_str0(NULL,"prior", "{STV,SAD,ATV}", "specify the prior you are applying from {STV,SAD,TV} , default as TV");
    struct arg_str  *bp    = arg_str0(NULL,"bp", "{Voxelbased,Raybased}", "specify the backprojection methods you are applying from {Voxelbased,Raybased}, default as voxelbased");
  
	struct arg_lit  *help    = arg_lit0("h","help",                    "print this help and exit");
	struct arg_lit  *version = arg_lit0(NULL,"version",                 "print version information and exit");
	struct arg_end  *end     = arg_end(20);


	void* argtable[] = {d_sigma,d_tau,d_lambda,d_tp,d_bc,d_huber,alp,outfile,projfiles,vols,
		imgWH,vs,ds,sid,sdd,offsetXYZ,startdegree,nframe,nrounds,AlgoIter,SartIter,prior, bp,help,version,end};
	const char* progname = "Test";
	int nerrors;
	int exitcode=0;

	/* verify the argtable[] entries were allocated sucessfully */
	if (arg_nullcheck(argtable) != 0)
	{
		/* NULL entries were detected, some allocations m_must have failed */
		printf("%s: insufficient memory\n",progname);
		exitcode=1;
		goto exit;
	}

	/* set any command line default values prior to parsing */
	
	outfile->filename[0]="Test/test";
	d_sigma->dval[0]=0.3;
	d_tau->dval[0]=0.3;
	d_lambda->dval[0]=0.1;
	d_tp->dval[0] = 0.0;
	d_bc->dval[0] = 0.0;
	d_huber->dval[0] = 0.1;
	AlgoIter->ival[0]=2;
	SartIter->ival[0]=2;
	alp->dval[0]=0.1;
	vs->dval[0]=1.0;
	sdd->dval[0]=200.0f;
	vols->ival[0] = 63;
	vols->ival[1] = 63;
	vols->ival[2] = 63;
	offsetXYZ->count = 3;
	offsetXYZ->dval[0] = 0.0f;
	offsetXYZ->dval[1] = 0.0f;
	offsetXYZ->dval[2] = 0.0f;

	// add prior
	 prior->sval[0]  = "ATV";     /* --prior={STV,SAD,ATV} */
	 bp->sval[0]  = "Voxelbased";
	sid->dval[0]=400.0;
	projfiles->filename[0]="proj";
	ds->dval[0]=1.0;
	nframe->ival[0] = 1;
	nrounds->ival[0] = 1;


	nerrors = arg_parse(argc,argv,argtable);

	/* special case: '--help' takes precedence over error reporting */
	if (help->count > 0)
	{
		printf("Usage: %s", progname);
		arg_print_syntax(stdout,argtable,"\n");
		printf("This program demonstrates the use of the argtable2 library\n");
		printf("for parsing command line arguments. Argtable accepts integers\n");
		printf("in decimal (123), hexadecimal (0xff), octal (0o123) and binary\n");
		printf("(0b101101) formats. Suffixes KB, MB and GB are also accepted.\n");
		arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		exitcode=0;
		goto exit;
	}

	/* special case: '--version' takes precedence error reporting */
	if (version->count > 0)
	{
		printf("'%s' example program for the \"argtable\" command line argument parser.\n",progname);
		printf("September 2003, Stewart Heitmann\n");
		exitcode=0;
		goto exit;
	}

	/* If the parser returned any errors then display them and exit */
	if (nerrors > 0)
	{
		/* Display the error details contained in the arg_end struct.*/
		arg_print_errors(stdout,end,progname);
		printf("Try '%s --help' for more information.\n",progname);
		exitcode=1;
		goto exit;
	}

	/* special case: uname with no command line options induces brief help */
	if (argc==1)
	{
		printf("Try '%s --help' for more information.\n",progname);
		exitcode=0;
		goto exit;
	}

	glutInit(&argc,argv);

	VolumeData* data=new VolumeData;
	data->setParas(vols->ival, vs->dval[0], nframe->ival[0]);

	data->ReadStructuredFile(); 
	
	// For each round, 10 projection images are applied.
	int proNo = 10 * nrounds->ival[0];
	ST_Tomography *scene = new ST_Tomography(d_sigma->dval[0], d_tau->dval[0], d_lambda->dval[0], d_tp->dval[0], d_bc->dval[0], d_huber->dval[0], outfile->filename[0],
		imgWH->ival, proNo, AlgoIter->ival[0], SartIter->ival[0], alp->dval[0],
		sid->dval[0], projfiles->filename, offsetXYZ->dval, sdd->dval[0], ds->dval[0], prior->sval[0], 
		bp->sval[0], vols->ival, nframe->ival[0], startdegree->dval,nrounds->ival[0]);
	scene->SetData(data);
	GLUTSceneContainer *mainWindow=GLUTSceneContainer::GetMainWindow(imgWH->ival,ds->dval[0]);
	mainWindow->SetTitle("S&T-tomography Framework");
	mainWindow->SetScene(scene);
	
	mainWindow->Create(10,850);	

exit:
	/* deallocate each non-null entry in argtable[] */
	arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
	//system("pause");
	return 0;
}

