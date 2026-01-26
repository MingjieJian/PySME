#ifndef SME_DLL
#ifdef BUILDING_SME_WIN_DLL
#define SME_DLL __declspec(dllexport)
#else
#define SME_DLL
#endif
#endif

// The SME library version (and compilation date)
#ifndef VERSION
#define VERSION "6.13, June 2025"
#endif

/* Datafile locations */
// DATA_DIR is defined in platform.h

#define DATAFILE_FE "Fe1_Bautista2017.dat.INTEL"
#define DATAFILE_NH "NH_Stancil2018.dat.INTEL"
#define DATAFILE_STEHLE "stehle_long.dat.INTEL"
#define DATAFILE_BPO "bpo_self.grid.INTEL"
#define DATAFILE_VCS "vcsbalmer.dat"

#ifndef IDL_DEFINE
#define IDL_DEFINE
// Define IDL String
typedef int IDL_STRING_SLEN_T;
#define IDL_STRING_MAX_SLEN 2147483647

typedef struct
{                         /* Define string descriptor */
  IDL_STRING_SLEN_T slen; /* Length of string, 0 for null */
  short stype;            /* type of string, static or dynamic */
  char *s;                /* Addr of string */
} IDL_STRING;
#endif
// define global parameter access
extern "C" int SME_DLL GetNLINES(void);
extern "C" short SME_DLL GetNRHOX(void);
extern "C" char *SME_DLL GetSPNAME(void);

// define the external methods
extern "C" const char *SME_DLL SMELibraryVersion(int n, void *arg[]); /* Return SME library version */
extern "C" const char *SME_DLL GetDataFiles(int n, void *arg[]);      /* Return the required data files */
extern "C" const char *SME_DLL GetLibraryPath(int n, void *arg[]);    /* Return the current data file directory */
extern "C" const char *SME_DLL SetLibraryPath(int n, void *arg[]);    /* Set the data file directory */
extern "C" const char *SME_DLL InputWaveRange(int n, void *arg[]);    /* Read in Wavelength range */
extern "C" const char *SME_DLL SetVWscale(int n, void *arg[]);        /* Set van der Waals scaling factor */
extern "C" const char *SME_DLL SetH2broad(int n, void *arg[]);        /* Set flag for H2 molecule */
extern "C" const char *SME_DLL ClearH2broad(int n, void *arg[]);      /* Clear flag for H2 molecule */
extern "C" const char *SME_DLL InputLineList(int n, void *arg[]);     /* Read in line list */
extern "C" const char *SME_DLL OutputLineList(int n, void *arg[]);    /* Return line list */
extern "C" const char *SME_DLL UpdateLineList(int n, void *arg[]);    /* Change line list parameters */
extern "C" const char *SME_DLL InputModel(int n, void *arg[]);        /* Read in model atmosphere */
extern "C" const char *SME_DLL InputDepartureCoefficients(int n, void *arg[]);
extern "C" const char *SME_DLL GetDepartureCoefficients(int n, void *arg[]);   /* Get NLTE b's for specific line */
extern "C" const char *SME_DLL GetNLTEflags(int n, void *arg[]);               /* Get line list NLTE flags */
extern "C" const char *SME_DLL ResetDepartureCoefficients(int n, void *arg[]); /* Reset LTE */
extern "C" const char *SME_DLL InputAbund(int n, void *arg[]);                 /* Read in abundances */
extern "C" const char *SME_DLL Opacity(int n, void *arg[]);                    /* Calculate opacities */
extern "C" const char *SME_DLL GetOpacity(int n, void *arg[]);                 /* Returns specific cont. opacity */
extern "C" const char *SME_DLL GetFraction(int n, void *arg[]);                /* Returns fraction or partition function */
extern "C" const char *SME_DLL Ionization(int n, void *arg[]);                 /* Perfrom EOS calculations */
extern "C" const char *SME_DLL GetDensity(int n, void *arg[]);                 /* Returns density in g/cm^3 */
extern "C" const char *SME_DLL GetNatom(int n, void *arg[]);                   /* Returns atomic number density */
extern "C" const char *SME_DLL GetNelec(int n, void *arg[]);                   /* Returns electron number density */
extern "C" const char *SME_DLL Transf(int n, void *arg[]);                     /* Computes spectral synthesis */
extern "C" const char *SME_DLL CentralDepth(int n, void *arg[]);               /* Computes line central depths */
extern "C" const char *SME_DLL GetLineOpacity(int n, void *arg[]);             /* Returns specific line opacity */
extern "C" const char *SME_DLL GetLineRange(int n, void *arg[]);               /* Get validity range for every line */
extern "C" const char *SME_DLL Contribution_functions(int n, void *arg[]);     /* Get contribution functions */