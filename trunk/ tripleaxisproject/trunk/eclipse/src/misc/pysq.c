/*
 * pysq.c --
 *
 *    This module implements "structure factor calculation" Python extension.
 *
 */

/*
  include data structures
  and locate Python.h and other build dependencies
*/

static char progid[] = "# pysq extension for python" ;

#include "pysq.h"


/*-----------------------------------------------------------------------
 *
 * This initial hack has most variables as following GLOBALS
 *
 * ----------------------------------------------------------------------
 */

/*
 * global dictionary for atomDEFs
 * and global storage ATOMSlist and Qlist
 */
static PyDictObject *atomDEFdict = NULL ;
static PyDictObject *spaceGRPdict = NULL ;
static PyDictObject *spacegrpALIASdict = NULL ;
static PyObject *pysqdict = NULL ;
static PyObject *scriptsdict = NULL ;

static ATOMSlist atomslist = {0, 0, NULL} ;
static Qlist qlist = {0, 0, NULL} ;

/* init the counters for atoms and Qs */

static int nAtomDef = 0 ;


/* put atom flags */
static Flag putflags[] = {
  {1,  "spacegroupActive"},
  {0,  "resetspacechange"},
  {0,  "autoprune"},
  {-3, "positionTolerance"},
  {0,  "multicell"},
  {0,  "amin"},
  {1,  "amax"},
  {0,  "bmin"},
  {1,  "bmax"},
  {0,  "cmin"},
  {1,  "cmax"}
} ;
static int nputflags = 11 ;

static int *usespacegroup = &(putflags[0].i) ;
static int *spacereset = &(putflags[1].i) ;
static int *autoprune = &(putflags[2].i) ;
static int *atomtolerance = &(putflags[3].i) ;
static int *multicell = &(putflags[4].i) ;
static int *amin = &(putflags[5].i) ;
static int *amax = &(putflags[6].i) ;
static int *bmin = &(putflags[7].i) ;
static int *bmax = &(putflags[8].i) ;
static int *cmin = &(putflags[9].i) ;
static int *cmax = &(putflags[10].i) ;


static int allpositions = 1 ;

/*
 allpositions flag will generate all general positions for putAtom and prune
*/
//static double atomDiffTolerance = 0.001 ;
static int readAtomic = 0 ;
static int readP = 0 ;
static int putAtomic = 0 ;
static int writeAtomic = 0 ;


static Flag calcflags[] = {
  {0, "nuconly"},
  {0, "magonly"},
  {0, "powder"},
  {1, "reset"},
  {1, "print"},
  {1, "atomsPrint"},
  {1, "setupPrint"},
  {1, "perQ"},
  {0, "pointers"},
  {0, "selectQ"},
  {1, "sortQ"},
  {0, "selectAtoms"},
} ;
static int ncalcflags = 12 ;

static Flag lorentzflags[] = {
  {0, "applyTocalc"},
  {1, "Qsteps"},
  {1, "twotheta"},
  {0, "powder"},
} ;
static int nlorentzflags = 4 ;

static Flag calccolumns[] = {
  {1, "Hcol"},
  {2, "Kcol"},
  {3, "Lcol"},
  {4, "Qcol"},
  {5, "QPcol"},
  {0, "SQPcol"},
  {0, "SQMcol"},
  {6, "SQPPcol"},
  {7, "SQPMcol"},
  {0, "SQMMcol"},
  {0, "SQMPcol"},
  {8, "SQcol"},
  {0, "TTcol"},
  {0, "OMcol"},
  {0, "LORcol"},
} ;
static int ncalccolumns = 15 ;

static Flag unknowncolumns[] = {
  {2, "atomname"},
  {3, "X"},
  {4, "Y"},
  {5, "Z"},
  {0, "Biso"},
  {0, "occ"},
  {0, "mom"},
  {0, "sx"},
  {0, "sy"},
  {0, "sz"},
} ;
static int nunknowncolumns = 10 ;

static int *nuconly = &(calcflags[0].i) ;
static int *magonly = &(calcflags[1].i) ;
static int *powderavg = &(calcflags[2].i) ;
static int *resetcalc = &(calcflags[3].i) ;
static int *calcprint = &(calcflags[4].i) ;
static int *printAtoms = &(calcflags[5].i) ;
static int *printSetup = &(calcflags[6].i) ;
static int *perQ = &(calcflags[7].i) ;
static int *calcpointers = &(calcflags[8].i) ;
static int *selectQ = &(calcflags[9].i) ;
static int *sortQ = &(calcflags[10].i) ;
static int *selectAtoms = &(calcflags[11].i) ;

static int *lorentzcalc = &(lorentzflags[0].i) ;
static int *Qlorentz = &(lorentzflags[1].i) ;
static int *TTlorentz = &(lorentzflags[2].i) ;
static int *Plorentz = &(lorentzflags[3].i) ;

/*
dont need calccolumns pointers
static int *Hcol = &(calccolumns[0].i) ;
static int *Kcol = &(calccolumns[1].i) ;
static int *Lcol = &(calccolumns[2].i) ;
static int *Qcol = &(calccolumns[3].i) ;
static int *QPcol = &(calccolumns[4].i) ;
static int *SQPcol = &(calccolumns[5].i) ;
static int *SQMcol = &(calccolumns[6].i) ;
static int *SQPPcol = &(calccolumns[7].i) ;
static int *SQPMcol = &(calccolumns[8].i) ;
static int *SQMMcol = &(calccolumns[9].i) ;
static int *SQMPcol = &(calccolumns[10].i) ;
static int *SQcol = &(calccolumns[11].i) ;
static int *TTcol = &(calccolumns[12].i) ;
static int *OMcol = &(calccolumns[13].i) ;
*/

/*
  output storage
*/
static double *SQh = NULL ;
static double *SQk = NULL ;
static double *SQl = NULL ;
static double *SQq = NULL ;
static double *SQp = NULL ;
static double *SQcs = NULL ;
static double *SQcsp = NULL ;
static double *SQcsm = NULL ;
static double *SQpp = NULL ;
static double *SQmm = NULL ;
static double *SQpm = NULL ;
static double *SQmp = NULL ;
static double *SQtt = NULL ;
static double *SQom = NULL ;
static double *SQlor = NULL ;

/* associate output columns and labels */
typedef struct {
  char *l ;
  double **d ;
} outCol ;
static outCol outcols[] = {
  {"      H     ", &SQh},
  {"      K     ", &SQk},
  {"      L     ", &SQl},
  {"      Q     ", &SQq},
  {"  Qoutplane ", &SQp},
  {"    SQ+     ", &SQcsp},
  {"    SQ-     ", &SQcsm},
  {"    ++      ", &SQpp},
  {"    +-      ", &SQpm},
  {"    --      ", &SQmm},
  {"    -+      ", &SQmp},
  {"  SQtotal   ", &SQcs},
  {"  twotheta  ", &SQtt},
  {"   omega    ", &SQom},
  {" Lorentz-fac", &SQlor}
} ;

static Int_List selectQlist = {0, NULL} ;
static Int_List selectAlist = {0, NULL} ;
/* wordlist for reading external files */
static WordList wordlist = {NULL, NULL, NULL, 0, 0, 0} ;

static char flagfile[] = { "pysqflags.py" } ;
static int calcOK = 0 ;    /* set 1 after calc, set 0 when params change */

static char titlebuf[512] = { "" } ;

/* init the spacegroup selection */
static char curSpaceGrp[64] = { "" } ;
static int  curSpaceGrpN = -1 ;
static int  NspaceGrpsLoaded = 0 ;
static char spacegroupsFile[512] ;
static char spacegroupsAliasFile[512] ;
static Spcgrp *spcgrps = NULL ;

static double lambda = 2.359 ;
static double ki ;
static double Ei ;

static double pbx[3] = XVEC ;
static double pby[3] = YVEC ;
static double pbz[3] = ZVEC ;

static int pbq = 0 ;

static double trx[3][3], try[3][3], trz[3][3] ;

static double avec[3] = { TWOPI, 0., 0. } ;
static double bvec[3] = { 0., TWOPI, 0. } ;
static double cvec[3] = { 0., 0., TWOPI } ;
static double ast[3] = XVEC ;
static double bst[3] = YVEC ;
static double cst[3] = ZVEC ;

static double latt[3] = { TWOPI, TWOPI, TWOPI } ;
static double angl[3] = { 90., 90., 90. } ;
static double rlat[3] = { 1., 1., 1. } ;
static double rang[3] = { 90., 90., 90. } ;

static double auni[3] = XVEC ;
static double buni[3] = YVEC ;
static double cuni[3] = ZVEC ;

static double orig[3] = { 0., 0., 0. } ;

#define  smallestSQ     0.0001
#define  smallestQ      0.00001
static double Bohrmagtob = 0.539 ;

static Q ip1 = { XVEC, XVEC } ;
static Q ip2 = { YVEC, YVEC } ;
static Q ip3 = { ZVEC, ZVEC } ;


static int NQ = 0 ;
static int NQalloc = 0 ;

static int Nelem = 103 ;
static char *ElemSymbols[] = {
  "H D T",                                                                                              "He",
  "Li", "Be",                                                             "B",  "C",  "N",  "O",  "F",  "Ne",
  "Na", "Mg",                                                             "Al", "Si", "P",  "S",  "Cl", "Ar",
  "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
  "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
  "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
  "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
  "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lw"
} ;


/* forward declarations for this modules functions */

typedef PyObject* (PyCmd) (PyObject *self, PyObject *args);

static PyCmd title;
static PyCmd lattice;
static PyCmd orient;
static PyCmd polarize;
static PyCmd spacegroup;
static PyCmd alias;
static PyCmd wavelength;
static PyCmd energy;
static PyCmd setupFile;

static PyCmd atomDef;
static PyCmd atomDefFF;
static PyCmd atomDefDW;
static PyCmd atomDefCopy;
static PyCmd atomDefRead;
static PyCmd atomDefFFRead;
static PyCmd atomDefWrite;
static PyCmd bEdep;

static PyCmd atomPut;
static PyCmd atomPutFlags;
//static PyCmd origin;
static PyCmd atomDel;
static PyCmd atomPrune;
static PyCmd atomList;
static PyCmd atomPrint;
static PyCmd atomFile;
static PyCmd unknownColumns;

static PyCmd calcFlags;
static PyCmd lorentzFlags;
static PyCmd calcColumns;
static PyCmd calc;
static PyCmd calcFile;

static PyCmd editAtom;
static PyCmd editMom;
static PyCmd editMag;
static PyCmd editFF;
static PyCmd editDW;
static PyCmd editOcc;
static PyCmd editDir;
static PyCmd editPos;
static PyCmd edit;
static PyCmd editByPos;

static PyCmd Qadd;
static PyCmd Qi;
static PyCmd Qindices;
static PyCmd Qdelete;
static PyCmd Qsort;
static PyCmd Qlisting;
static PyCmd QFile;

static PyCmd readFile;
static PyCmd saveFlags;
static PyCmd allFlags;

static struct PyMethodDef sq_methods[] = {
  {"title", title, 1},
  {"lattice", lattice, 1}, 
  {"orient", orient, 1},
  {"polarize", polarize, 1},
  {"spacegroup", spacegroup, 1},
  {"alias", alias, 1},
  {"wavelength", wavelength, 1},
  {"energy", energy, 1},
  {"setupSave", setupFile, 1},
  {"setupsave", setupFile, 1},
  {"atomdef", atomDef, 1},
  {"bEdep", bEdep, 1},
  {"atomdefFF", atomDefFF, 1},
  {"atomdefDW", atomDefDW, 1},
  {"atomdefCopy", atomDefCopy, 1},
  {"atomdefRead", atomDefRead, 1},
  {"atomdefFFRead", atomDefFFRead, 1},
  {"atomdefWrite", atomDefWrite, 1},
  {"atomdefFile", atomDefWrite, 1},
  {"atomPut", atomPut, 1},
  {"atomput", atomPut, 1},
  {"put", atomPut, 1},
  {"putflags", atomPutFlags, 1},
  {"atomDel", atomDel, 1},
  {"atomdel", atomDel, 1},
  {"delete", atomDel, 1},
  {"prune", atomPrune, 1},
  {"atomList", atomList, 1},
  {"atomlist", atomList, 1},
  {"editAtom", editAtom, 1},
  {"edit", edit, 1},
  {"editPos", editPos, 1},
  {"editDir", editDir, 1},
  {"editMom", editMom, 1},
  {"editMag", editMag, 1},
  {"editFF", editFF, 1},
  {"editDW", editDW, 1},
  {"mag", editMag, 1},
  {"editOcc", editOcc, 1},
  {"editByPos", editByPos, 1},
  {"atomPrint", atomPrint, 1},
  {"atomprint", atomPrint, 1},
  {"atomSave", atomFile, 1},
  {"atomsave", atomFile, 1},
  {"al", atomPrint, 1},
  {"calc", calc, 1},
  {"calcSave", calcFile, 1},
  {"calcsave", calcFile, 1},
  {"calcflags", calcFlags, 1},
  {"lorentz", lorentzFlags, 1},
  {"calccolumns", calcColumns, 1},
  {"Qadd", Qadd, 1},
  {"Q", Qadd, 1},
  {"Qi", Qi, 1},
  {"Qgen", Qindices, 1},
  {"indices", Qindices, 1},
  {"Qdelete", Qdelete, 1},
  {"Qsort", Qsort, 1},
  {"Qlist", Qlisting, 1},
  {"Ql", Qlisting, 1},
  {"QSave", QFile, 1},
  {"Qsave", QFile, 1},
  {"readFile", readFile, 1},
  {"readfile", readFile, 1},
  {"unknowncolumns", unknownColumns, 1},
  {"flagsave", saveFlags, 1},
  {"flagSave", saveFlags, 1},
  {"flags", allFlags, 1},
  {NULL, NULL}
} ;

static void pbcalc() ;
static int readAtomDefFile(char *fil) ;
static int readAtomDefFileFF(char *fil) ;
static int readAliases(char *fil, PyDictObject *d, PyDictObject *g) ;
static int readPYSQ(char *fil, int ph) ;

/* module init */
void
initpysq()
{
  static char location[] = "initpysq" ;
  PyObject *pysqmod ;
  PyObject *pmod, *pdict, *pathlist, *item ;
  PyObject *scriptsmod, *mdict, *pymod ;
  PyObject *value, *pyflags ;

  int i, j, nread ;

  char atomdefFile[512] ;
  char *spacegroupsPtr, *atomFilePtr ;

  atomDEFdict = (PyDictObject*) PyDict_New() ;
  spaceGRPdict = (PyDictObject*) PyDict_New() ;
  spacegrpALIASdict = (PyDictObject*) PyDict_New() ;

  pysqmod = Py_InitModule3("pysq", sq_methods, "pysq module") ;
  pysqdict = PyModule_GetDict(pysqmod) ;

  /* import the pyfitScripts which are mostly sys and os calls */
  scriptsmod = PyImport_ImportModule("pysqScripts") ;
  if( scriptsmod == NULL ) {
    printf("failed pysqScripts import\n") ;
    exit ;
  }
  scriptsdict = PyModule_GetDict(scriptsmod) ;
  /*
    now with PyRun_String(command, Py_eval_input, scriptsdict, scriptsdict)
    should be able to run scripts
  */

  /* check if the sys.modules dictionary contains pysq */
  pmod = PyImport_ImportModule("sys") ;
  if( pmod == NULL ) {
    printf("failed sys import\n") ;
    exit ;
  }
  pdict = PyModule_GetDict(pmod) ;
  /* now get the modules dictionary from the sys dictionary */
  mdict = PyDict_GetItemString(pdict, "modules") ;
  /* check if pyfit is in this dictionary */
  if ( (pymod = PyDict_GetItemString(mdict, "pysq")) == NULL ) {
    printf("couldnt find pysq in sys.modules\n") ;
  }

  /*
    in principle we could now extract the path to pyfit from
    the returned value pymod
    for now we are going to assume modules need to be loaded from
    cwd so insert cwd at head of sys.path
  */

  item = PyRun_String("cwd()", Py_eval_input, scriptsdict, scriptsdict) ;

  pathlist = PyDict_GetItemString(pdict, "path") ;
  //item = Py_BuildValue("s", "./") ;
  PyList_Insert(pathlist, 0, item) ;
  PyDict_SetItemString(pdict, "path", pathlist) ;

  /* find the spacegroups.dat file */
  item = PyRun_String("findsysfile('spacegroups.dat')",
		      Py_eval_input, scriptsdict, scriptsdict) ;
  spacegroupsPtr = PyString_AsString(item) ;
  if( strlen(spacegroupsPtr) < 1 ) {
    Py_DECREF(item) ;
    printf("Cant find spacegroups.dat file\n") ;
    exit(0) ;
  }
  strcpy(spacegroupsFile, spacegroupsPtr) ;
  Py_DECREF(item) ;

  if( NspaceGrpsLoaded < 230 ) {
    spcgrps = (Spcgrp *)calloc(NSPACEGROUPS, sizeof(Spcgrp)) ;
    for( i=0 ; i<NSPACEGROUPS ; i++ ) {
      for( j=0 ; j<4 ; j++ ) spcgrps[i].trans[j].ptstr = NULL ;
    }
    if( memExc(isNULL(spcgrps), location) ) {
      printf("failed to allocate spacegroups\n") ;
      exit(0) ;
      //return NULL ;
    }
    nread = readWyck(spacegroupsFile, spcgrps, NSPACEGROUPS) ;
    NspaceGrpsLoaded = nread ;
    if( nread < 230 ) {
      printf("failed to read spacegroups\n") ;
      exit(0) ;
      //PyErr_SetString(PyExc_ValueError,"failed to read all spacegroups") ;
      //return NULL ;
    }

    /* set deflt spacegroup */
    strcpy(curSpaceGrp,spcgrps[0].symbol) ;
    curSpaceGrpN = 0 ;

    for( i=0 ; i<nread ; i++ ) {
      if ((value=PyDict_GetItemString((PyObject *)spaceGRPdict,
				      spcgrps[i].symbol))
	  == NULL) {
	value = PyCObject_FromVoidPtr((void*)(spcgrps+i), NULL) ;
	if ( PyDict_SetItemString((PyObject *)spaceGRPdict,
				  spcgrps[i].symbol, value) ) {
	  printf("failed to load spacegroups into dictionary\n") ;
	  exit(0) ;
	  //PyErr_SetString(PyExc_ValueError,
	  //	  "failed to put spacegroups in dictionary") ;
	  //return NULL ;
	}
      }
    }
  }

  /* find the spacegroups.aliases file */
  item = PyRun_String("findsysfile('spacegroups.aliases')",
		      Py_eval_input, scriptsdict, scriptsdict) ;
  spacegroupsPtr = PyString_AsString(item) ;
  if( strlen(spacegroupsPtr) < 1 ) {
    printf("Cant find spacegroups.aliases file\n") ;
    Py_DECREF(item) ;
    exit(0) ;
  }
  strcpy(spacegroupsAliasFile, spacegroupsPtr) ;
  Py_DECREF(item) ;
  nread = readAliases(spacegroupsAliasFile, spacegrpALIASdict, spaceGRPdict) ;

  /* now read in the atom definitions */
  /* find the atomDef.dat file */
  item = PyRun_String("findsysfile('atomDef.dat')",
		      Py_eval_input, scriptsdict, scriptsdict) ;
  atomFilePtr = PyString_AsString(item) ;
  if( strlen(atomFilePtr) < 1 ) {
    printf("Cant find atomDef.dat file\n") ;
    Py_DECREF(item) ;
    exit(0) ;
  } else {
    //printf("about to load atomDefFile=%s\n", atomFilePtr) ;
  }
  strcpy(atomdefFile, atomFilePtr) ;
  Py_DECREF(item) ;
  if( readAtomDefFile(atomdefFile) < 0 ) {
    printf("Failed to load atomDef.dat file\n") ;
    exit(0) ;
  }

  /* now read in the magnetic form factors */
  /* find the magff.dat file */
  item = PyRun_String("findsysfile('magff.dat')",
		      Py_eval_input, scriptsdict, scriptsdict) ;
  atomFilePtr = PyString_AsString(item) ;
  if( strlen(atomFilePtr) < 1 ) {
    printf("Cant find magff.dat file\n") ;
    Py_DECREF(item) ;
    exit(0) ;
  } else {
    //printf("about to load magff=%s\n", atomFilePtr) ;
  }
  strcpy(atomdefFile, atomFilePtr) ;
  Py_DECREF(item) ;
  if( readAtomDefFileFF(atomdefFile) < 0 ) {
    printf("Failed to load magff.dat file\n") ;
    exit(0) ;
  }

  /* try to read a local pysqflags.py file */
  readPYSQ("pysqflags.py",0) ;
  if( PyErr_Occurred() ) PyErr_Clear() ;

  ki = TWOPI/lambda ;
  Ei = DNEUT * ki*ki ;
  pbcalc() ;
}




/* SETUP commands */

static PyObject *
title(PyObject *self, PyObject *args)
{
  int narg ;
  char *str ;
  PyObject *value ;

  /*
   * title(titlestring)
   */

  narg = PyTuple_Size(args) ;
  if( narg > 0 ) {
    if ( (value = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
    if ( ! PyString_Check(value) ) {
      PyErr_SetString(PyExc_ValueError, "title arg must be string") ;
      return NULL ;
    } 
    str = PyString_AsString(value) ;
    if ( str == NULL ) {
      PyErr_SetString(PyExc_ValueError, "title arg converted to NULL") ;
      return NULL ;
    }
    if( strlen(str) > 511 ) str[511] = '\0' ;
    strcpy(titlebuf, str) ;
  }
  return Py_BuildValue("s", titlebuf) ;
}


static int recip() ;
static int checkLattice(int fix) ;

static PyObject *
lattice(PyObject *self, PyObject *args)
{
  int i, narg, ok ;
  double dbl ;
  PyObject *Llist, *Alist ;
  PyObject *value, *src ;
  static char *lattcheck[] = {
    "lattice and spacegroup incompatible",
    "lattice and spacegroup compatible"
  } ;

  /*
   * lattice( a ?b c ?bcAng acAng abAng?? )
   */

  narg = PyTuple_Size(args) ;
  src = args ;
  if( narg > 0 ) {
    if ( (value = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
    if ( PyTuple_Check(value) ) {
      src = value ;
      narg = PyTuple_Size(value) ;
    } 
  }

  for( i=0 ; i<(narg<=6?narg:6) ; i++ ) {
    if ( (value = PyTuple_GetItem(src, i)) == NULL ) return NULL ;
    if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError, "args to lattice must be numeric") ;
      return NULL ;
    }

    dbl = PyFloat_AsDouble(value) ;
    if( i<3 ) latt[i] = dbl ;
    else angl[i-3] = dbl ;
  }
  if( narg == 1 ) { latt[1] = latt[0] ; latt[2] = latt[0] ; }
  if( narg > 0 && narg < 4 ) { angl[0] = angl[1] = angl[2] = 90. ; }
  if( narg == 4 ) { angl[1] = angl[0] ; angl[2] = angl[0] ; }


  if( !recip() ) {
    PyErr_SetString(PyExc_ValueError, "lattice failed recip lattice setup") ;
    return NULL ;
  }
  if( narg > 0 ) calcOK = 0 ;
  //Llist = MakeDblList(3, latt) ;
  //Alist = MakeDblList(3, angl) ;
  ok = checkLattice(0) ;
  return Py_BuildValue("ddddddss",
		       latt[0],latt[1],latt[2],
		       angl[0],angl[1],angl[2],
		       "lattice params in Angstroms and angles in degrees",
		       lattcheck[ok]) ;
}

static void rhsetup(double hkl1[3], double hkl2[3], double hkl3[3]) ;

static PyObject *
orient(PyObject *self, PyObject *args)
{
  int i, narg ;
  double dbl ;
  PyObject *list1, *list2 ;
  PyObject *value, *src ;

  /*
   * orient( h,k,l,h,k,l )
   */

  narg = PyTuple_Size(args) ;
  src = args ;
  if( narg > 0 ) {
    if ( (value = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
    if ( PyTuple_Check(value) ) {
      src = value ;
      narg = PyTuple_Size(value) ;
    } 
  }

  if( narg > 0 && narg < 6 ) {
    PyErr_SetString(PyExc_ValueError, "orient requires hkl1 hkl2 = 6 args") ;
    return NULL ;
  }

  for( i=0 ; i<(narg<=6?narg:6) ; i++ ) {
    if ( (value = PyTuple_GetItem(src, i)) == NULL ) return NULL ;
    if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError, "args to orient must be numeric") ;
      return NULL ;
    }

    dbl = PyFloat_AsDouble(value) ;
    if( i<3 ) ip1.hkl[i] = dbl ;
    else ip2.hkl[i-3] = dbl ;
  }

  rhsetup(ip1.hkl, ip2.hkl, ip2.hkl) ;
  /* whenever orient changes must recalc polarization transform to xtal coord*/
  pbcalc() ;
  if( narg > 0 ) calcOK = 0 ;
  //list1 = MakeDblList(3, ip1.hkl) ;
  //list2 = MakeDblList(3, ip2.hkl) ;
  return Py_BuildValue("dddddds",
		       ip1.hkl[0],ip1.hkl[1],ip1.hkl[2],
		       ip2.hkl[0],ip2.hkl[1],ip2.hkl[2],
		       "orientation recip lattice vectors hkl1 and hkl2") ;
}

static PyObject *
polarize(PyObject *self, PyObject *args)
{
  int i, narg ;
  PyObject *list ;
  PyObject *value, *src ;

  /*
   * polarize( alongQ, pbzX, pbzY, pbzZ)
   */

  narg = PyTuple_Size(args) ;
  src = args ;
  if( narg > 0 ) {
    if ( (value = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
    if ( PyTuple_Check(value) ) {
      src = value ;
      narg = PyTuple_Size(value) ;
    } 
  }

  if( narg > 0 && narg < 4 ) {
    PyErr_SetString(PyExc_ValueError,
		    "polarize requires alongQ X Y Z = 4 args") ;
    return NULL ;
  }

  if ( narg > 0 ) {
    if ( (value = PyTuple_GetItem(src, 0)) == NULL ) return NULL ;
    if ( ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError, "polarize first arg is alongQ flag") ;
      return NULL ;
    }
    pbq = PyInt_AsLong(value) ;
    for( i=1 ; i<4 ; i++ ) {
      if ( (value = PyTuple_GetItem(src, i)) == NULL ) return NULL ;
      if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError, "args to polarize must be numeric") ;
	return NULL ;
      }

      pbz[i-1] = PyFloat_AsDouble(value) ;
    }
    pbcalc() ;
    calcOK = 0 ;
  }

  //list = MakeDblList(3, pbz) ;
  return Py_BuildValue("iddds", pbq, pbz[0], pbz[1], pbz[2],
		       "alongQFlag pbzX pbzY pbzZ") ;
}

static PyObject *
wavelength(PyObject *self, PyObject *args)
{
  int narg ;
  PyObject *Llist, *Alist ;
  PyObject *value, *src ;

  /*
   * wavelength(lambda)
   */

  narg = PyTuple_Size(args) ;
  src = args ;
  if( narg > 0 ) {
    if ( (value = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
    if ( PyTuple_Check(value) ) {
      src = value ;
      narg = PyTuple_Size(value) ;
    } 
  }

  if( narg > 0 ) {
    if ( (value = PyTuple_GetItem(src, 0)) == NULL ) return NULL ;
    if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError, "arg to wavelength must be numeric") ;
      return NULL ;
    }
    lambda = PyFloat_AsDouble(value) ;
  }
  if( lambda > 0. ) { ki = TWOPI/lambda ; Ei = DNEUT * ki*ki ; }
  if( narg > 0 ) calcOK = 0 ;
  return Py_BuildValue("ds", lambda, "wavelength in Angstroms") ;
}
static PyObject *
energy(PyObject *self, PyObject *args)
{
  int narg ;
  PyObject *Llist, *Alist ;
  PyObject *value, *src ;

  /*
   * energy(Emev)
   */

  narg = PyTuple_Size(args) ;
  src = args ;
  if( narg > 0 ) {
    if ( (value = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
    if ( PyTuple_Check(value) ) {
      src = value ;
      narg = PyTuple_Size(value) ;
    } 
  }

  if( narg > 0 ) {
    if ( (value = PyTuple_GetItem(src, 0)) == NULL ) return NULL ;
    if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError, "arg to energy must be numeric") ;
      return NULL ;
    }
    Ei = PyFloat_AsDouble(value) ;
    lambda = PyFloat_AsDouble(value) ;
  }
  if( Ei > 0. ) { ki = sqrt(Ei/ DNEUT ) ; lambda = TWOPI/ki ; }
  if( narg > 0 ) calcOK = 0 ;
  return Py_BuildValue("ds", Ei, "energy in meV") ;
}

static DCMPLX bGSAS(DCMPLX b, double E, int nc, double *c) ;
static PyObject *
bEdep(PyObject *self, PyObject *args)
{
  /*
    args: atomName, Estart, Efinal, Estep
  */
  int i, narg, nE ;
  char *name ;
  double Es[3], Etemp, En ;
  DCMPLX b ;
  ATOMdef *atomDefPtr ;
  PyObject *value, *atom ;

  narg = PyTuple_Size(args) ;
  if( narg < 3 ) {
    PyErr_SetString(PyExc_ValueError, "bEdep requires atomName, Ei, Ef args") ;
    return NULL ;
  }

  if ( (value = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;

  if( ! PyString_Check(value) ) {
    PyErr_SetString(PyExc_ValueError, "bEdep requires atomName string") ;
    return NULL ;
  }

  if( (atom = PyDict_GetItem((PyObject *)atomDEFdict, value)) == NULL ) {
    PyErr_SetString(PyExc_ValueError,
		    "that ATOM is not in ATOMdef dictionary") ;
    return NULL ;
  }
  atomDefPtr = (ATOMdef *) PyCObject_AsVoidPtr(atom) ;
  if ( atomDefPtr == NULL ) {
    if ( PyErr_Occurred() ) return NULL ;
    PyErr_SetString(PyExc_ValueError,
		    "failed to convert to ATOMdef pointer") ;
    return NULL ;
  }
  name = PyString_AsString(value) ;
  
  Es[2] = 0. ;
  for( i=1 ; i<(narg < 4 ? narg : 4) ; i++ ) {
    if ( (value = PyTuple_GetItem(args, i)) == NULL ) return NULL ;
    if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError,
		      "bEdep args after atomName must be numeric") ;
      return NULL ;
    }
    Es[i-1] = PyFloat_AsDouble(value) ;
  }
  if( Es[0] > Es[1] ) { Etemp = Es[0] ; Es[0] = Es[1] ; Es[1] = Etemp ; }
  if( Es[2] <= 0. ) {
    nE = 21 ;
    Es[2] = (Es[1] - Es[0])/20. ;
  } else {
    nE = (int)((Es[1] - Es[0])/Es[2]) + 1 ;
    if( nE > 1 ) Es[2] = (Es[1] - Es[0])/(nE - 1) ;
  }
  printf("b(E) for atom=%s\n    E      bReal      bImag\n", name) ;
  for( i=0 ; i<nE ; i++ ) {
    En = Es[0] + i*Es[2] ;
    b = bGSAS(atomDefPtr->b, En, atomDefPtr->nbcoef, atomDefPtr->bcoef) ;
    printf("%7f %10.4g %10.4g\n", En, b.r, b.i) ;
  }
  return Py_BuildValue("") ;
}

static void resetAtomsList(ATOMSlist *atomslst) ;
static int getSpaceGroupIndex(int ispc, int iset) ;
static char *stripstring(char *s) ;

static PyObject *
spacegroup(PyObject *self, PyObject *args)
{
  static char location[] = "spacegroup" ;
  /*
    spacegroup(ispace,[iset] OR name)
  */

  PyObject *arg, *slist, *sublist, *ptlist, *value, *tlist ;

  static char one[2] = "1" ;
  char *sitesym ;
  int ok, narg, i, j, ptlen, ispc, iset ;
  int nread ;
  char coords[128] ;
  char buf[128] ;
  char *str, *ptPtr ;
  static char *lattcheck[] = {
    "lattice repaired due to spacegroup incompatibility",
    "lattice is compatible with spacegroup"
  } ;
  Spcgrp *spcgrpPtr ;
  Subgrp *subgrpPtr ;

  narg = PyTuple_Size(args) ;

  if( narg < 1 ) {

    spcgrpPtr = spcgrps + curSpaceGrpN ;


    if( spcgrpPtr->origin[0] == '\0' ) {
      if( (value = Py_BuildValue("sii", curSpaceGrp, spcgrpPtr->spcnum,
				 spcgrpPtr->setting[0])) == NULL ) return NULL;
    } else {
      if( (value = Py_BuildValue("siis", curSpaceGrp, spcgrpPtr->spcnum,
				 spcgrpPtr->setting[0], spcgrpPtr->origin))
	  == NULL ) return NULL;
    }

    if( (slist = PyList_New(2+spcgrpPtr->nsub)) == NULL ) return NULL ;
    PyList_SET_ITEM(slist, 0, value) ;

    if( (tlist = PyList_New(spcgrpPtr->ntrans)) == NULL ) {
      Py_DECREF(slist) ; return NULL ;
    }

    for( i=0 ; i<spcgrpPtr->ntrans ; i++ ) {
      ptPtr = spcgrpPtr->trans[i].ptstr ;
      if( ptPtr == NULL || (ptlen = strlen(ptPtr)) < 3 ) continue ;
      strcpy(coords, ptPtr+1) ;
      coords[ptlen-2] = '\0' ;
      if( (value = Py_BuildValue("s", coords)) == NULL ) return NULL ;
      PyList_SET_ITEM(tlist, i, value) ;
    }
    PyList_SET_ITEM(slist, 1, tlist) ;

    for( i=0 ; i<spcgrpPtr->nsub ; i++ ) {
      subgrpPtr = spcgrpPtr->sub + i ;
      if( (sublist = PyList_New(2)) == NULL ) {
	Py_DECREF(slist) ; return NULL ;
      }
      if( strlen(subgrpPtr->sitesym) < 12 ) sitesym = subgrpPtr->sitesym ;
      else sitesym = one ;
      if( (value = Py_BuildValue("ics",
				 subgrpPtr->npts, subgrpPtr->Wyck, sitesym))
	  == NULL ) return NULL ;
      PyList_SET_ITEM(sublist, 0, value) ;

      if( (ptlist = PyList_New(subgrpPtr->npts)) == NULL ) {
	Py_DECREF(slist) ; return NULL ;
      }
      for( j=0 ; j<subgrpPtr->npts ; j++ ) {
	ptPtr = subgrpPtr->pts[j].ptstr ;
	if( ptPtr == NULL || (ptlen = strlen(ptPtr)) < 3 ) continue ;
	strcpy(coords, ptPtr+1) ;
	coords[ptlen-2] = '\0' ;
	if( (value = Py_BuildValue("s", coords)) == NULL ) {
	  Py_DECREF(slist) ; return NULL ;
	}
	PyList_SET_ITEM(ptlist, j, value) ;
      }
      PyList_SET_ITEM(sublist, 1, ptlist) ;
      PyList_SET_ITEM(slist, 2+i, sublist) ;
    }
    return Py_BuildValue("Ns",slist,
	  " (spcgrpName,num,set),[trans],[[(npts,Wyck,sym),ptslist],..]") ;
    /* passing slist to Py_BuildValue use N to avoid incr its ref count */
  }

  if ( (arg = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
  if( PyInt_Check(arg) ) {
    ispc = PyInt_AsLong(arg) ;
    if( ispc < 1 || ispc > NspaceGrps ) {
      PyErr_SetString(PyExc_ValueError, "spacegroup number out of range") ;
      return NULL ;
    }
    iset = 1 ;
    if( narg > 1 ) {
      if ( (arg = PyTuple_GetItem(args, 1)) == NULL ) return NULL ;
      if( ! PyInt_Check(arg) ) {
	PyErr_SetString(PyExc_ValueError,
			"spacegroup arg after spc number is int setting") ;
	return NULL ;
      }
      iset = PyInt_AsLong(arg) ;
    }
    curSpaceGrpN = getSpaceGroupIndex(ispc, iset) ;
    spcgrpPtr = spcgrps + curSpaceGrpN ;
    strcpy(curSpaceGrp, spcgrps[curSpaceGrpN].symbol) ;
  } else if( PyString_Check(arg) ) {
    if( (str = PyString_AsString(arg)) == NULL ) return NULL ;
    strcpy(buf,str) ;
    str = stripstring(buf) ;
    for( i=1 ; i<strlen(str) ; i++ ) str[i] = tolower(str[i]) ;
    if ((value=PyDict_GetItemString((PyObject *)spaceGRPdict, str)) == NULL) {
      /* check for an alias */
      if ((value=PyDict_GetItemString((PyObject *)spacegrpALIASdict, str))
	  == NULL) {
	PyErr_SetString(PyExc_ValueError,
			"cant find that spacegroup name") ;
	return NULL ;
      }
      /* found alias */
      str = PyString_AsString(value) ;
      if ((value=PyDict_GetItemString((PyObject *)spaceGRPdict, str))
	  == NULL) {
	PyErr_SetString(PyExc_ValueError,
			"cant find that spacegroup name or alias") ;
	return NULL ;
      }
    }
    spcgrpPtr = (Spcgrp *)PyCObject_AsVoidPtr(value) ;
    strcpy(curSpaceGrp, spcgrpPtr->symbol) ;
    curSpaceGrpN = spcgrpPtr - spcgrps ;
  }

  if( atomslist.n > 0 && *spacereset ) {
    resetAtomsList(&atomslist) ;
  }
  ok = checkLattice(1) ;
  return Py_BuildValue("iiss", spcgrpPtr->spcnum, spcgrpPtr->setting[0],
		       curSpaceGrp, lattcheck[ok]) ;
}

static PyObject *
alias(PyObject *self, PyObject *args)
{
  static char location[] = "alias" ;
  /*
    alias(spacegrpname)
  */

  PyObject *key, *value, *tpl, *arg, *slist ;
  char buf[128] ;
  char *keystr, *valuestr ;
  int narg, pos, i ;

  narg = PyTuple_Size(args) ;

  if( narg < 1 ) {
    /* return list of all alias pairs */
    pos = 0 ;
    i = 0 ;
    if((slist = PyList_New(PyDict_Size((PyObject*)spacegrpALIASdict))) == NULL)
     return NULL ;
    while (PyDict_Next((PyObject*)spacegrpALIASdict, &pos, &key, &value)) {
      keystr = PyString_AsString(key) ;
      valuestr = PyString_AsString(value) ;
      if( (tpl = Py_BuildValue("ss", keystr, valuestr)) == NULL ) {
	Py_DECREF(slist) ; return NULL ;
      }
      PyList_SET_ITEM(slist, i, tpl) ;
      i++ ;
    }
    return slist ;
  }

  /* lookup the arg */
  if ( (arg = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
  if( ! PyString_Check(arg) ) {
    PyErr_SetString(PyExc_ValueError, "alias arg must be string") ;
    return NULL ;
  }
  keystr = PyString_AsString(arg) ;
  strcpy(buf,keystr) ;
  keystr = stripstring(buf) ;
  //for( i=1 ; i<strlen(str) ; i++ ) str[i] = tolower(str[i]) ;
  if ((value=PyDict_GetItemString((PyObject *)spacegrpALIASdict, keystr))
      == NULL) {
    PyErr_SetString(PyExc_ValueError,
		    "cant find that spacegroup name") ;
    return NULL ;
  }
  /* found alias */
  valuestr = PyString_AsString(value) ;
  return Py_BuildValue("s", valuestr) ;
}


/* ATOM DEFINITION commands */

static int atomLookupP(PyObject *key, PyObject **value, ATOMdef **aptr) ;
static int atomLookupS(char *key, PyObject **value, ATOMdef **aptr) ;
static ATOMdef *atomNewP(PyObject *key) ;
static ATOMdef *atomNewS(char *key) ;
static PyObject *atomdefTuple(char *name, ATOMdef *adef) ;
static PyObject *listFromDbls(int n, double *d) ;


static PyObject *
atomDef(PyObject *self, PyObject *args)
{
  /*
    atomDef( name/atomicnumFORreturnvalues br bi m brNonRes bresCoefs..",},
    with no args returns 2-tuple: 1st item is list of atomDef tuples
    2nd item is description of atomDef data
    so to call atomDef with its return value
    aa = atomDef()
    atomDef(aa[0])
  */

  PyObject *alist, *atom ;
  PyObject *arg, *key, *value, *newt, *src, *Clist ;

  int atomicnum ;
  int narg, i, istat ;
  int ff, dw, atomnum, natom, copy ;
  int new ;
  int pos ;
  double br, bi, m ;

  char *name ;
  ATOMdef *atomDefPtr ;
  static ATOMdef *atomDefPtrCopy = NULL ;

  narg = PyTuple_Size(args) ;

  if( narg < 1 ) {
    /* return list of all defined atoms */
    pos = 0 ;
    i = 0 ;
    if ( (alist = PyList_New(PyDict_Size((PyObject*)atomDEFdict))) == NULL )
	 return NULL ;
    while (PyDict_Next((PyObject*)atomDEFdict, &pos, &key, &value)) {
      atomDefPtr = (ATOMdef *) PyCObject_AsVoidPtr(value) ;
      if ( atomDefPtr == NULL ) {
	Py_DECREF(alist) ;
	if ( PyErr_Occurred() ) return NULL ;
	PyErr_SetString(PyExc_MemoryError,
			"failed to convert to ATOMdef pointer") ;
	return NULL ;
      }
      name = PyString_AsString(key) ;

      if( (atom = atomdefTuple(name, atomDefPtr)) == NULL ) return NULL ;
      PyList_SET_ITEM(alist, i, atom) ;
      i++ ;
    }
    return Py_BuildValue("Ns", alist,
			 " defined atoms: name, atomicnum, br, bi, mom") ;
  }

  src = args ;
  if ( (arg = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;

  if ( PyList_Check(arg) ) {
    /*
      called with list of tuples
      extract each tuple and call atomDef with each tuple
    */
    narg = PyList_Size(arg) ;
    for( i=0 ; i<narg ; i++ ) {
      if ( (newt = PyList_GetItem(arg, i)) == NULL ) return NULL ;
      if ( ! PyTuple_Check(newt) ) {
	PyErr_SetString(PyExc_ValueError,
			"atomDef list items must be tuples") ;
	return NULL ;
      }
      if ( PyTuple_Size(newt) < 2 ) {
	PyErr_SetString(PyExc_ValueError,
     "atomDef list items must be tuples with at least name bR") ;
	return NULL ;
      }
      if( ! atomDef(self, newt) ) return NULL ;
    }
    return Py_BuildValue("") ;
  }

  if( narg == 1 ) {
    /* return that atomdef */
    if( PyString_Check(arg) ) {
      if( (value = PyDict_GetItem((PyObject *)atomDEFdict, arg)) == NULL ) {
	PyErr_SetString(PyExc_ValueError,
			"that ATOM is not in ATOMdef dictionary") ;
	return NULL ;
      }
      atomDefPtr = (ATOMdef *) PyCObject_AsVoidPtr(value) ;
      if ( atomDefPtr == NULL ) {
	if ( PyErr_Occurred() ) return NULL ;
	PyErr_SetString(PyExc_ValueError,
			"failed to convert to ATOMdef pointer") ;
	return NULL ;
      }
      name = PyString_AsString(arg) ;
      if( (atom = atomdefTuple(name, atomDefPtr)) == NULL ) return NULL ;
      return atom ;
    } else if( PyInt_Check(arg) ) {
      /* return all definitions for that atomic num */
      atomicnum = PyInt_AsLong(arg) ;
      pos = 0 ;
      i = 0 ;
      if ( (alist = PyList_New(PyDict_Size((PyObject*)atomDEFdict))) == NULL )
	   return NULL ;
      while (PyDict_Next((PyObject*)atomDEFdict, &pos, &key, &value)) {
	atomDefPtr = (ATOMdef *) PyCObject_AsVoidPtr(value) ;
	if ( atomDefPtr == NULL ) {
	  Py_DECREF(alist) ;
	  if ( PyErr_Occurred() ) return NULL ;
	  PyErr_SetString(PyExc_ValueError,
			  "failed to convert to ATOMdef pointer") ;
	  return NULL ;
	}
	if( atomDefPtr->atomic != atomicnum ) continue ;
	name = PyString_AsString(key) ;
	if( (atom = atomdefTuple(name, atomDefPtr)) == NULL ) return NULL ;
	PyList_SET_ITEM(alist, i, atom) ;
	i++ ;
      }
      return Py_BuildValue("Ns", alist,
	  " defined atoms for that atomicnum: name, atomicnum, br, bi, mom") ;
    } else if( PyTuple_Check(arg) ) {
      src = arg ;
      narg = PyTuple_Size(arg) ;
      if( narg < 2 ) {
	PyErr_SetString(PyExc_ValueError,
			"atomDef called with tuple requires at least 2 args") ;
	return NULL ;
      }
      if ( (arg = PyTuple_GetItem(arg, 0)) == NULL ) return NULL ;
    } else {
      PyErr_SetString(PyExc_ValueError,
		      "invalid arg type to atomDef") ;
      return NULL ;
    }
  }

  /*
   * Handle bR bI moment [atomicnum] args
   */

  if ( ! PyString_Check(arg) ) {
    PyErr_SetString(PyExc_ValueError,
		    "first arg to atomDef must be atom name string") ;
    return NULL ;
  }
  name = PyString_AsString(arg) ;

  if ( (istat = atomLookupS(name, &value, &atomDefPtr)) < 0 ) return NULL ;
  if ( istat == 0 && (atomDefPtr = atomNewS(name)) == NULL ) return NULL ;

  for( i=1 ; i<narg ; i++ ) {
    if ( (value = PyTuple_GetItem(src, i)) == NULL ) return NULL ;
    if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError,
		      "args to atomDef after name must be numeric") ;
      return NULL ;
    }
    if( i==1 ) atomDefPtr->b.r = PyFloat_AsDouble(value) ;
    else if( i==2 ) atomDefPtr->b.i = PyFloat_AsDouble(value) ;
    else if( i==3 ) atomDefPtr->m = PyFloat_AsDouble(value) ;
    else if( i==4 && PyInt_Check(value) )
      atomDefPtr->atomic = PyInt_AsLong(value) ;
    else if( i>4 && i<14 ) atomDefPtr->bcoef[i-5] = PyFloat_AsDouble(value) ; 
  }
  atomDefPtr->nbcoef = 0 ;
  /* must be at least 5 bcoefs to define wavelength dep */
  if( narg >= 10 ) atomDefPtr->nbcoef = narg - 5 ;
  if( (atom = atomdefTuple(name, atomDefPtr)) == NULL ) return NULL ;
  return atom ;
}

static PyObject *MakeDblList(int n, double *dlst) ;
static int CopyAtomDef(ATOMdef *dest, ATOMdef *src) ;

static PyObject *
atomDefFF(PyObject *self, PyObject *args)
{
  static char location[] = "atomDefFF" ;
  /*
    atomDefFF( name/atomicnumber, j2frac, magffcoefs )
    with no args returns 2-tuple: 1st item is list of FF tuples
    2nd item is description of FF data
    so to call atomDefFF with its return value
    ff = atomDefFF()
    atomDefFF(ff[0])
  */

  PyObject *alist, *atom, *fflist ;
  PyObject *arg, *key, *value, *src, *newt ;

  int atomicnum ;
  int narg, i, istat, nret ;
  int ff, dw, atomnum, natom, copy ;
  int offset, isList ;
  int pos ;
  double br, bi, m ;

  char *name ;
  ATOMdef *atomDefPtr ;
  MagFormFactor *mf ;
  static ATOMdef *atomDefPtrCopy = NULL ;

  narg = PyTuple_Size(args) ;

  if( narg < 1 ) {
    /* return list of all defined atoms magnetic magformfactors */
    pos = 0 ;
    i = 0 ;
    nret = 0 ;
    while (PyDict_Next((PyObject*)atomDEFdict, &pos, &key, &value)) {
      atomDefPtr = (ATOMdef *) PyCObject_AsVoidPtr(value) ;
      if ( atomDefPtr == NULL ) {
	if ( PyErr_Occurred() ) return NULL ;
	PyErr_SetString(PyExc_MemoryError,
			"failed to convert to ATOMdef pointer") ;
	return NULL ;
      }
      if( (mf = atomDefPtr->ff) == NULL ) continue ;
      nret++ ;
    }

    if ( (alist = PyList_New(nret)) == NULL ) return NULL ;
    pos = 0 ;
    i = 0 ;
    while (PyDict_Next((PyObject*)atomDEFdict, &pos, &key, &value)) {
      atomDefPtr = (ATOMdef *) PyCObject_AsVoidPtr(value) ;
      if( (mf = atomDefPtr->ff) == NULL ) continue ;
      name = PyString_AsString(key) ;
      fflist = MakeDblList(mf->n, mf->coefs) ;
      if( (atom = Py_BuildValue("sdN", name, mf->j2frac, fflist)) == NULL )
	return NULL ;
      PyList_SET_ITEM(alist, i, atom) ;
      i++ ;
    }
    return Py_BuildValue("Ns", alist,
			 " defined atoms: name, j2frac, magFFcoefs") ;
  }

  src = args ;
  if ( (arg = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;

  if ( PyList_Check(arg) ) {
    /*
      called with list of tuples
      extract each tuple and call atomDefFF with each tuple
    */
    narg = PyList_Size(arg) ;
    for( i=0 ; i<narg ; i++ ) {
      if ( (newt = PyList_GetItem(arg, i)) == NULL ) return NULL ;
      if ( ! PyTuple_Check(newt) ) {
	PyErr_SetString(PyExc_ValueError,
			"atomDefFF list items must be tuples") ;
	return NULL ;
      }
      if ( PyTuple_Size(newt) < 2 ) {
	PyErr_SetString(PyExc_ValueError,
     "atomDefFF list items must be tuples with at least name j2frac") ;
	return NULL ;
      }
      if( ! atomDefFF(self, newt) ) return NULL ;
    }
    return Py_BuildValue("") ;
  }

  if( narg == 1 ) {
    /* return that atomdef FF */
    if( PyString_Check(arg) ) {
      if( (value = PyDict_GetItem((PyObject *)atomDEFdict, arg)) == NULL ) {
	PyErr_SetString(PyExc_ValueError,
			"that ATOM is not in ATOMdef dictionary") ;
	return NULL ;
      }
      atomDefPtr = (ATOMdef *) PyCObject_AsVoidPtr(value) ;
      if ( atomDefPtr == NULL ) {
	if ( PyErr_Occurred() ) return NULL ;
	PyErr_SetString(PyExc_ValueError,
			"failed to convert to ATOMdef pointer") ;
	return NULL ;
      }
      name = PyString_AsString(arg) ;
      if( (mf = atomDefPtr->ff) != NULL ) {
	fflist = MakeDblList(mf->n, mf->coefs) ;
	if( (atom = Py_BuildValue("sdN", name, mf->j2frac, fflist)) == NULL )
	  { Py_DECREF(fflist) ; return NULL ; }
      } else {
	if( (fflist = PyList_New(0)) == NULL ) return NULL ;
	if( (atom = Py_BuildValue("sdN", name, 0., fflist)) == NULL )
	  return NULL ;
      }
      return atom ;
    } else if( PyInt_Check(arg) ) {
      /* return all definition FFs for that atomic num */
      atomicnum = PyInt_AsLong(arg) ;
      pos = 0 ;
      i = 0 ;
      while (PyDict_Next((PyObject*)atomDEFdict, &pos, &key, &value)) {
	atomDefPtr = (ATOMdef *) PyCObject_AsVoidPtr(value) ;
	if ( atomDefPtr == NULL ) {
	  if ( PyErr_Occurred() ) return NULL ;
	  PyErr_SetString(PyExc_ValueError,
			  "failed to convert to ATOMdef pointer") ;
	  return NULL ;
	}
	if( atomDefPtr->atomic != atomicnum ) continue ;
	i++ ;
      }
      nret = i ;

      pos = 0 ;
      i = 0 ;
      if ( (alist = PyList_New(nret)) == NULL ) return NULL ;
      while (PyDict_Next((PyObject*)atomDEFdict, &pos, &key, &value)) {
	atomDefPtr = (ATOMdef *) PyCObject_AsVoidPtr(value) ;
	if( atomDefPtr->atomic != atomicnum ) continue ;
	name = PyString_AsString(key) ;
	if( (mf = atomDefPtr->ff) != NULL ) {
	  fflist = MakeDblList(mf->n, mf->coefs) ;
	} else if( (fflist = PyList_New(0)) == NULL ) {
	  Py_DECREF(alist) ;
	  return NULL ;
	}
	if( (atom = Py_BuildValue("sdN", name, mf->j2frac, fflist)) == NULL ) {
	  Py_DECREF(alist) ; return NULL ;
	}
	PyList_SET_ITEM(alist, i, atom) ;
      }
      return Py_BuildValue("Ns", alist,
	       " defined atoms FF for that atomicnum: name, j2frac, FFcoefs") ;
    } else if( PyTuple_Check(arg) ) {
      src = arg ;
      narg = PyTuple_Size(arg) ;
      if( narg < 2 ) {
	PyErr_SetString(PyExc_ValueError,
		    "atomdefFF called with tuple requires at least 2 args") ;
	return NULL ;
      }
      if ( (arg = PyTuple_GetItem(arg, 0)) == NULL ) return NULL ;
    } else {
      PyErr_SetString(PyExc_ValueError,
		      "invalid arg type to atomDefFF") ;
      return NULL ;
    }
  }

  /*
   *  narg > 1  Handle j2frac, FFcoef args
   */

  if ( ! PyString_Check(arg) ) {
    PyErr_SetString(PyExc_ValueError,
		    "first arg to atomDefFF must be atom name string") ;
    return NULL ;
  }
  name = PyString_AsString(arg) ;

  if ( (istat = atomLookupS(name, &value, &atomDefPtr)) < 0 ) return NULL ;
  if ( istat == 0 && (atomDefPtr = atomNewS(name)) == NULL ) return NULL ;
  if ( istat == 0 && atomDefPtr->atomic > 0 ) {
    if ( (istat = atomLookupS(ElemSymbols[atomDefPtr->atomic-1],
			      &value, &atomDefPtrCopy)) < 0 ) return NULL ;
    if ( istat > 0 ) if(! CopyAtomDef(atomDefPtr, atomDefPtrCopy)) return NULL;
  }

  if ( atomDefPtr->ff == NULL ) {
    atomDefPtr->ff = (MagFormFactor *)calloc(1,sizeof(MagFormFactor)) ;
    if ( memExc(isNULL(atomDefPtr->ff), location) ) return NULL ;
  }
  /* we know there is at least one data arg */
  if ( (value = PyTuple_GetItem(src, 1)) == NULL ) return NULL ;
  if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError,
		"args to atomDefFF after name must be numeric j2frac,coeffs") ;
      return NULL ;
  }
  atomDefPtr->ff->j2frac = PyFloat_AsDouble(value) ;

  if ( narg < 3 ) {
    fflist = MakeDblList(atomDefPtr->ff->n, atomDefPtr->ff->coefs) ;
    if( (atom = Py_BuildValue("sdN", name, atomDefPtr->ff->j2frac, fflist))
      == NULL ) { Py_DECREF(fflist) ; return NULL ; }
    return atom ;
  }

  if ( (arg = PyTuple_GetItem(src, 2)) == NULL ) return NULL ;
  //src = args ;
  isList = 0 ;
  offset = 2 ;
  narg -= 2 ;
  if ( PyList_Check(arg) ) {
    src = arg ;
    isList = 1 ;
    offset = 0 ;
    narg = PyList_Size(arg) ;
  } else if( PyTuple_Check(arg) ) {
    src = arg ;
    offset = 0 ;
    narg = PyTuple_Size(arg) ;
  }

  for( i=0 ; i<narg ; i++ ) {
    if ( ! isList ) {
      if ( (value = PyTuple_GetItem(src, i+offset)) == NULL ) return NULL ;
    } else {
      if ( (value = PyList_GetItem(src, i+offset)) == NULL ) return NULL ;
    }
    if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError,
		      "args to atomDefFF after name must be numeric") ;
      return NULL ;
    }
    atomDefPtr->ff->coefs[i] = PyFloat_AsDouble(value) ;
  }
  atomDefPtr->ff->n = narg ;

  fflist = MakeDblList(atomDefPtr->ff->n, atomDefPtr->ff->coefs) ;
  if( (atom = Py_BuildValue("sdN", name, atomDefPtr->ff->j2frac, fflist))
    == NULL ) { Py_DECREF(fflist) ; return NULL ; }
  return atom ;
}

static PyObject *MakeDWList(DebyeWallerFactor *dw) ;

static PyObject *
atomDefDW(PyObject *self, PyObject *args)
{
  static char location[] = "atomDefDW" ;
  /*
    atomDefDW( name/atomicnumber, dwcoefsORlist )
  */

  PyObject *alist, *atom, *fflist ;
  PyObject *arg, *key, *value, *src, *oldsrc ;

  int atomicnum ;
  int narg, i, j, row, col, ni, nret ;
  int ff, atomnum, natom, copy ;
  int isList, pos, istat ;
  double br, bi, m ;

  char *name ;
  ATOMdef *atomDefPtr ;
  static ATOMdef *atomDefPtrCopy = NULL ;
  DebyeWallerFactor *dw ;

  narg = PyTuple_Size(args) ;

  if( narg < 1 ) {
    /* return list of all defined atoms dw coeffs */
    pos = 0 ;
    i = 0 ;
    while (PyDict_Next((PyObject*)atomDEFdict, &pos, &key, &value)) {
      atomDefPtr = (ATOMdef *) PyCObject_AsVoidPtr(value) ;
      if ( atomDefPtr == NULL ) {
	if ( PyErr_Occurred() ) return NULL ;
	PyErr_SetString(PyExc_ValueError,
			"failed to convert to ATOMdef pointer") ;
	return NULL ;
      }
      if( atomDefPtr->dw == NULL ) continue ;
      i++ ;
    }
    pos = 0 ;
    nret = i ;
    i = 0 ;
    if ( (alist = PyList_New(nret)) == NULL ) return NULL ;
    while (PyDict_Next((PyObject*)atomDEFdict, &pos, &key, &value)) {
      atomDefPtr = (ATOMdef *) PyCObject_AsVoidPtr(value) ;
      if( atomDefPtr->dw == NULL ) continue ;
      name = PyString_AsString(key) ;
      fflist = MakeDWList(atomDefPtr->dw) ;
      if( (atom = Py_BuildValue("sN", name, fflist)) == NULL ) return NULL ;
      PyList_SET_ITEM(alist, i, atom) ;
      i++ ;
    }
    return Py_BuildValue("Ns", alist,
			 " defined atoms: name, DWcoefs") ;
  }

  src = args ;
  if ( (arg = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;

  if( narg == 1 ) {
    /* return that atomdef DW */
    if( PyString_Check(arg) ) {
      if( (value = PyDict_GetItem((PyObject *)atomDEFdict, arg)) == NULL ) {
	PyErr_SetString(PyExc_ValueError,
			"that ATOM is not in ATOMdef dictionary") ;
	return NULL ;
      }
      atomDefPtr = (ATOMdef *) PyCObject_AsVoidPtr(value) ;
      if ( atomDefPtr == NULL ) {
	if ( PyErr_Occurred() ) return NULL ;
	PyErr_SetString(PyExc_ValueError,
			"failed to convert to ATOMdef pointer") ;
	return NULL ;
      }
      name = PyString_AsString(arg) ;
      if( atomDefPtr->dw != NULL ) {
	fflist = MakeDWList(atomDefPtr->dw) ;
      } else if( (fflist = PyList_New(0)) == NULL ) {
	return NULL ;
      }
      if( (atom = Py_BuildValue("sN", name, fflist)) == NULL ) {
	Py_DECREF(fflist) ; return NULL ;
      }
      return atom ;
    } else if( PyInt_Check(arg) ) {
      /* return all definition DWs for that atomic num */
      atomicnum = PyInt_AsLong(arg) ;
      pos = 0 ;
      i = 0 ;
      while (PyDict_Next((PyObject*)atomDEFdict, &pos, &key, &value)) {
	atomDefPtr = (ATOMdef *) PyCObject_AsVoidPtr(value) ;
	if ( atomDefPtr == NULL ) {
	  if ( PyErr_Occurred() ) return NULL ;
	  PyErr_SetString(PyExc_ValueError,
			  "failed to convert to ATOMdef pointer") ;
	  return NULL ;
	}
	if( atomDefPtr->atomic != atomicnum ) continue ;
	i++ ;
      }
      pos = 0 ;
      nret = i ;
      i = 0 ;
      if ( (alist = PyList_New(nret)) == NULL ) return NULL ;
      while (PyDict_Next((PyObject*)atomDEFdict, &pos, &key, &value)) {
	atomDefPtr = (ATOMdef *) PyCObject_AsVoidPtr(value) ;
	if( atomDefPtr->atomic != atomicnum ) continue ;
	name = PyString_AsString(key) ;
	if( atomDefPtr->dw != NULL )
	  fflist = MakeDWList(atomDefPtr->dw) ;
	else if( (fflist = PyList_New(0)) == NULL ) return NULL ;
	if( (atom = Py_BuildValue("sO", name, fflist)) == NULL ) return NULL ;
	PyList_SET_ITEM(alist, i, atom) ;
	i++ ;
      }
      return Py_BuildValue("Os", alist,
		    " defined atoms DW for that atomicnum: name, DWmatrix") ;
    } else if( PyTuple_Check(arg) ) {
      src = arg ;
      narg = PyTuple_Size(arg) ;
      if( narg < 2 ) {
	PyErr_SetString(PyExc_ValueError,
		    "atomdefDW called with tuple requires at least 2 args") ;
	return NULL ;
      }
      if ( (arg = PyTuple_GetItem(src, 0)) == NULL ) return NULL ;
    } else {
      PyErr_SetString(PyExc_ValueError,
		      "invalid arg type to atomDefDW") ;
      return NULL ;
    }
  }

  /*
   *  narg > 1  Handle DW args
   *  src is the tuple and arg is first extracted arg which should be name
   */

  if ( ! PyString_Check(arg) ) {
    PyErr_SetString(PyExc_ValueError,
		    "first arg to atomDefDW must be atom name string") ;
    return NULL ;
  }
  name = PyString_AsString(arg) ;

  if ( (istat = atomLookupS(name, &value, &atomDefPtr)) < 0 ) return NULL ;
  if ( istat == 0 && (atomDefPtr = atomNewS(name)) == NULL ) return NULL ;
  if ( istat == 0 && atomDefPtr->atomic > 0 ) {
    if ( (istat = atomLookupS(ElemSymbols[atomDefPtr->atomic-1],
			      &value, &atomDefPtrCopy)) < 0 ) return NULL ;
    if ( istat > 0 ) if(! CopyAtomDef(atomDefPtr, atomDefPtrCopy)) return NULL;
  }

  if ( atomDefPtr->dw == NULL ) {
    atomDefPtr->dw = (DebyeWallerFactor *)calloc(1,sizeof(DebyeWallerFactor)) ;
    if ( memExc(isNULL(atomDefPtr->dw), location) ) return NULL ;
  }
  dw = atomDefPtr->dw ;

  if ( (arg = PyTuple_GetItem(src, 1)) == NULL ) return NULL ;

  if ( PyList_Check(arg) ) {
    src = PyList_AsTuple(arg) ;
  } else {
    oldsrc = src ;
    if( (src = PyTuple_New(narg-1)) == NULL ) return NULL ;
    for( i=1 ; i<narg ; i++ ) {
      if( (value = PyTuple_GetItem(oldsrc, i)) == NULL ) {
	Py_DECREF(src) ; return NULL ;
      }
      if( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
			"args to atomDefDW after name must be numeric") ;
	Py_DECREF(src) ;
	return NULL ;
      }
      PyTuple_SetItem(src, i-1, Py_BuildValue("d", PyFloat_AsDouble(value))) ;
      //src = PyTuple_GetSlice(src, 1, narg) ;
    }
  }
  /* src is new ref */
  narg = PyTuple_Size(src) ;
  dw->iso = 0 ;
  if ( narg < 3 ) {
    /* iso ONLY */
    if ( narg > 0 ) {
      if ( (value = PyTuple_GetItem(src, 0)) == NULL ) {
	Py_DECREF(src) ; return NULL ;
      }
      if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
			"args to atomDefDW after name must be numeric") ;
	Py_DECREF(src) ;
	return NULL ;
      }
      dw->u[0][0] = PyFloat_AsDouble(value) ;
    } else {
      dw->u[0][0] = 0. ;
    }
    dw->iso = 1 ;
  } else if( narg < 6 ) {
    /* diagonal ONLY */
    for( i=0 ; i<3 ; i++ ) {
      if ( (value = PyTuple_GetItem(src, i)) == NULL ) {
	Py_DECREF(src) ; return NULL ;
      }
      if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
			"args to atomDefDW after name must be numeric") ;
	Py_DECREF(src) ;
	return NULL ;
      }
      for( j=0 ; j<3 ; j++ ) dw->u[i][j] = 0. ;
      dw->u[i][i] = PyFloat_AsDouble(value) ;
    }
  } else if( narg < 9 ) {
    /* symm u11 u12 u13 u22 u23 u33 */
    ni = 0 ;
    for( i=0 ; i<9 ; i++ ) {
      if ( i == 3 || i == 6 || i == 7 ) continue ;
      if ( (value = PyTuple_GetItem(src, ni)) == NULL ) {
	Py_DECREF(src) ; return NULL ;
      }
      if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
			"args to atomDefDW after name must be numeric") ;
	Py_DECREF(src) ;
	return NULL ;
      }
      row = i/3 ;
      col = i%3 ;
      dw->u[row][col] = PyFloat_AsDouble(value) ;
      dw->u[col][row] = dw->u[row][col] ;
      ni++ ;
    }
  } else {
    /* full anis u all 9 elems */
    for( i=0 ; i<9 ; i++ ) {
      if ( (value = PyTuple_GetItem(src, i)) == NULL ) {
	Py_DECREF(src) ; return NULL ;
      }
      if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
			"args to atomDefDW after name must be numeric") ;
	Py_DECREF(src) ;
	return NULL ;
      }
      row = i/3 ;
      col = i%3 ;
      dw->u[row][col] = PyFloat_AsDouble(value) ;
    }
  }
  Py_DECREF(src) ;
  fflist = MakeDWList(dw) ;
  if( (atom = Py_BuildValue("sN", name, fflist)) == NULL ) return NULL ;
  return atom ;
}

static PyObject *
atomDefRead(PyObject *self, PyObject *args)
{
  static char location[] = "atomDefRead" ;
  /*
    atomDefRead( file )
  */

  PyObject *arg ;

  int narg, natom ;
  char *name ;

  narg = PyTuple_Size(args) ;
  arg = NULL ;
  if( narg > 0 ) if ( (arg = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
  if( narg < 1 || ! PyString_Check(arg) ) {
    PyErr_SetString(PyExc_ValueError, "atomDefRead requires file name") ;
    return NULL ;
  }
  name = PyString_AsString(arg) ;
  if( (natom = readAtomDefFile(name)) < 0 ) return NULL ;
  return Py_BuildValue("iss", natom, " atoms read from ", name) ;
}

static PyObject *
atomDefFFRead(PyObject *self, PyObject *args)
{
  static char location[] = "atomDefFFRead" ;
  /*
    atomDefFFRead( file )
  */

  PyObject *arg ;

  int narg, natom ;
  char *name ;

  narg = PyTuple_Size(args) ;
  arg = NULL ;
  if( narg > 0 ) if ( (arg = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
  if( narg < 1 || ! PyString_Check(arg) ) {
    PyErr_SetString(PyExc_ValueError, "atomDefFFRead requires file name") ;
    return NULL ;
  }
  name = PyString_AsString(arg) ;
  if( (natom = readAtomDefFileFF(name)) < 0 ) return NULL ;
  return Py_BuildValue("iss", natom, " atoms read from ", name) ;
}


static int writeAtomDefFile(char *fil) ;

static PyObject *
atomDefWrite(PyObject *self, PyObject *args)
{
  static char location[] = "atomDefWrite" ;
  /*
    atomDefWrite( file )
  */

  PyObject *arg ;

  int narg, natom ;
  char *name ;

  narg = PyTuple_Size(args) ;
  arg = NULL ;
  if( narg > 0 ) if ( (arg = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
  if( narg < 1 || ! PyString_Check(arg) ) {
    PyErr_SetString(PyExc_ValueError, "atomDefWrite requires file name") ;
    return NULL ;
  }
  name = PyString_AsString(arg) ;
  if( (natom = writeAtomDefFile(name)) < 0 ) return NULL ;
  return Py_BuildValue("iss", natom, " atoms written to ", name) ;
}

static PyObject *
atomDefCopy(PyObject *self, PyObject *args)
{
  static char location[] = "atomDefCopy" ;
  /*
    atomDefCopy( destName, srcName )
  */

  PyObject *dest, *src, *value ;

  int narg, istat ;
  char *destname, *srcname ;

  ATOMdef *atomDefdest, *atomDefsrc ;

  narg = PyTuple_Size(args) ;
  if( narg > 1 ) {
    if ( (dest = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
    if ( (src = PyTuple_GetItem(args, 1)) == NULL ) return NULL ;
  }
  if( narg < 2 || ! PyString_Check(dest) || ! PyString_Check(src) ) {
    PyErr_SetString(PyExc_ValueError,
		    "atomDefCopy requires dest and src names") ;
    return NULL ;
  }
  destname = PyString_AsString(dest) ;
  srcname = PyString_AsString(src) ;

  /* lookup the source */

  if ( (istat = atomLookupS(srcname, &value, &atomDefsrc)) < 0 ) return NULL ;
  if ( istat == 0 ) {
    PyErr_SetString(PyExc_ValueError,
		    "src atom NOT in dictionary for copying");
    return NULL ;
  }

  /* check to see if dest is already in dict */

  if ( (istat = atomLookupS(destname, &value, &atomDefdest)) < 0 ) return NULL;
  if ( istat == 0 && (atomDefdest = atomNewS(destname)) == NULL ) return NULL ;

  if( ! CopyAtomDef(atomDefdest, atomDefsrc) ) return NULL ;
  return Py_BuildValue("ssss", "copied ", srcname, " to ", destname) ;
}


/* atomPut COMMAND */

static int putAtom(char *id, ATOMgroup *atomgroup,
		    ATOMdef *atomDefPtr, int isub,
		    double r[3], double s[3], double m, double c) ;
static ATOMgroup *nextAtomGroup(ATOMSlist *atomslst) ;
static int SymLookup(char S) ;
static int dblsFromTuple(PyObject *tpl, int n, double *d) ;
static int dblsFromList(PyObject *lst, int n, double *d) ;
static PyObject *MakeAtomList(char *name, int isub, myATOM *a) ;

static PyObject *
atomPut(PyObject *self, PyObject *args)
{
  /*
    atomPut( name/atomicnum, isub/symmSymbol,  x y z sx sy sz m c dwcoefs )
    NB x y z sx sy sz are in reduced coordinates
    args can also be lists for multiple atoms or list of such lists
    so that this command can read its output
  */

  PyObject *arg, *atom, *value, *item, *newt ;
  PyObject *atomsList, *groupList, *atomList ;

  int narg, i, j, isub, nmiss, atomicnum, istat, nc, nret ;
  int xyztyp, sxyztyp, nextarg ;
  double r[3], s[3], m, c ;
  char buf[128] ;
  char *name, *sym ;

  ATOMdef *atomDefPtr, *atomDefdest ;
  ATOMSlist *atomslst ;
  ATOMgroup *atomgroup, *atomgroupi ;
  myATOM *atoms ;

  narg = PyTuple_Size(args) ;

  atomicnum = 0 ;
  name = NULL ;
  sym = NULL ;
  isub = -1 ;

  if( narg > 0 ) {
    if ( (atom = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
    if ( PyList_Check(atom) ) {
      /* is this a list of lists or each arg is a list */
      if ( (item = PyList_GetItem(atom, 0)) == NULL ) return NULL ;
      if ( PyList_Check(item) ) {
	/*
	 * list of lists
	 * so for each list in the list call putAtom recursively
	 */
	narg = PyList_Size(atom) ;
	for( i=0 ; i<narg ; i++ ) {
	  if ( (item = PyList_GetItem(atom, i)) == NULL ) return NULL ;
	  if ( ! PyList_Check(item) ) {
	    PyErr_SetString(PyExc_ValueError,
			    "arg to atomPut should have beem list") ;
	    return NULL ;
	  }
	  newt = PyList_AsTuple(item) ;
	  if ( ! PyTuple_Check(newt) ) {
	    PyErr_SetString(PyExc_ValueError,
			    "failed to make tuple from list in atomPut") ;
	    Py_DECREF(newt) ;
	    return NULL ;
	  }
	  if ( atomPut(self, newt) == NULL ) { Py_DECREF(newt) ; return NULL ;}
	  Py_DECREF(newt) ;
	}
      } else {
	/*
	 * args are each lists
	 * elements of atom List should be standard args to atomPut
	 */
	for( i=0 ; i<narg ; i++ ) {
	  if ( (item = PyTuple_GetItem(args, i)) == NULL ) return NULL ;
	  if ( ! PyList_Check(item) ) {
	    PyErr_SetString(PyExc_ValueError,
			    "arg to atomPut should have beem list") ;
	    return NULL ;
	  }
	  newt = PyList_AsTuple(item) ;
	  if ( ! PyTuple_Check(newt) ) {
	    PyErr_SetString(PyExc_ValueError,
			    "failed to make tuple from list in atomPut") ;
	    Py_DECREF(newt) ;
	    return NULL ;
	  }
	  if ( atomPut(self, newt) == NULL ) { Py_DECREF(newt); return NULL ; }
	  Py_DECREF(newt) ;
	}
      }
      return Py_BuildValue("") ;
    }

    if ( PyString_Check(atom) ) name = PyString_AsString(atom) ;
    else if ( PyInt_Check(atom) ) atomicnum = PyInt_AsLong(atom) ;
    else {
      PyErr_SetString(PyExc_ValueError, "invalid first arg to atomPut") ;
      return NULL ;
    }
  }

  xyztyp = 0 ;
  if( narg > 1 ) {
    if ( (atom = PyTuple_GetItem(args, 1)) == NULL ) return NULL ;
    if ( PyString_Check(atom) ) {
      sym = PyString_AsString(atom) ;
    } else if ( PyInt_Check(atom) ) {
      isub = PyInt_AsLong(atom) ;
    } else if ( PyTuple_Check(atom) ) {
      if( dblsFromTuple(atom, 3, r) < 0 ) return NULL ;
      xyztyp = 1 ;
      nextarg = 2 ;
    } else if ( PyList_Check(atom) ) {
      if( dblsFromList(atom, 3, r) < 0 ) return NULL ;
      xyztyp = 2 ;
      nextarg = 2 ;
    } else {
      PyErr_SetString(PyExc_ValueError, "invalid 2nd arg to atomPut") ;
      return NULL ;
    }
  }
  if( sym != NULL ) isub = SymLookup(sym[0]) ;

  atomgroup = atomslist.atomslist ;
  if( (xyztyp == 0 && narg < 3) || (narg < 2) ) {
    /*
      return list of atoms from atomlist
      matching any supplied name/atomicnum or isub/sym
      only return list same as input, i.e. first atom of group
      use the atomList cmnd to get complete listing of all atoms
    */
    nret = 0 ;
    for( i=0 ; i<atomslist.n ; i++ ) {
      atomgroupi = atomgroup + i ;
      if( atomgroupi->atom == NULL ) continue ;
      if( atomicnum > 0 && atomgroupi->atom->atomic != atomicnum ) continue ;
      if( name != NULL && strcmp(name, atomgroupi->atom->name) ) continue ;
      if( sym != NULL && sym[0] != atomgroupi->Wsym ) continue ;
      if( isub >= 0 && isub != atomgroupi->isub ) continue ;
      nret++ ;
    }
    if ( (atomsList = PyList_New(nret)) == NULL ) return NULL ;
    nret = 0 ;
    for( i=0 ; i<atomslist.n ; i++ ) {
      atomgroupi = atomgroup + i ;
      if( atomgroupi->atom == NULL ) continue ;
      if( atomicnum > 0 && atomgroupi->atom->atomic != atomicnum ) continue ;
      if( name != NULL && strcmp(name, atomgroupi->atom->name) ) continue ;
      if( sym != NULL && sym[0] != atomgroupi->Wsym ) continue ;
      if( isub >= 0 && isub != atomgroupi->isub ) continue ;

      if ( (groupList = MakeAtomList(atomgroupi->atom->name, atomgroupi->isub,
				     atomgroupi->atoms)) == NULL )
	{ Py_DECREF(atomsList) ; return NULL ; }

      PyList_SET_ITEM(atomsList, nret, groupList) ;
      /* the groupList ref is stolen */
      nret++ ;
    }
    return atomsList ;
  }

  if( xyztyp == 0 && narg > 2 ) {
    if ( (atom = PyTuple_GetItem(args, 2)) == NULL ) return NULL ;
    if ( PyTuple_Check(atom) ) {
      if( dblsFromTuple(atom, 3, r) < 0 ) return NULL ;
      xyztyp = 1 ;
      nextarg = 3 ;
    } else if ( PyList_Check(atom) ) {
      if( dblsFromList(atom, 3, r) < 0 ) return NULL ;
      xyztyp = 2 ;
      nextarg = 3 ;
    } else {
      if( narg < 5 ) {
	PyErr_SetString(PyExc_ValueError,
			"atomPut requires at least x, y, z after name, isub") ;
	return NULL ;
      }
      for( i=0 ; i<3 ; i++ ) {
	if ( (value = PyTuple_GetItem(args, i+2)) == NULL ) return NULL ;
	if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	  PyErr_SetString(PyExc_ValueError,
			  "x y z args to atomPut must be numeric") ;
	  return NULL ;
	}
	r[i] = PyFloat_AsDouble(value) ;
      }
      nextarg = 5 ;
    }    
  }

  if( name != NULL ) {
    if ( (istat = atomLookupS(name, &value, &atomDefPtr)) < 0 ) return NULL ;
    if ( istat == 0 ) {
      /* make sure 2nd char is lc and try to find a shorter version */
      strcpy(buf, name) ;
      buf[1] = tolower(buf[1]) ;
      nc = strlen(buf) ;
      while( nc > 0 && istat == 0 ) {
	buf[nc] = '\0' ;
	if ( (istat = atomLookupS(buf, &value, &atomDefPtr)) < 0 )
	  return NULL ;
	if ( istat > 0 ) {
	  if ( (atomDefdest = atomNewS(name)) == NULL ) return NULL ;
	  if( ! CopyAtomDef(atomDefdest, atomDefPtr) ) return NULL ;
	  atomDefPtr = atomDefdest ;
	  break ;
	}
	nc-- ;
      }
      if ( istat == 0 ) {
	PyErr_SetString(PyExc_ValueError,
			"can't find atomname for put command") ;
	return NULL ;
      }
    }
  } else if( atomicnum > 0 && atomicnum <= Nelem ) {
    if ( (istat = atomLookupS(ElemSymbols[atomicnum-1],
			      &value, &atomDefPtr)) < 0 ) return NULL ;
    if ( istat == 0 ) {

      PyErr_SetString(PyExc_ValueError,
		      "can't find atomname for atomicnum in put command") ;
      return NULL ;
    }
  } else {
    PyErr_SetString(PyExc_ValueError, "No valid atom spec in put command") ;
    return NULL ;
  }
  
  for( i=0 ; i<3 ; i++ ) s[i] = 0. ;

  if( nextarg < narg ) {
    if ( (value = PyTuple_GetItem(args, nextarg)) == NULL ) return NULL ;

    if ( PyTuple_Check(value) ) {
      if( dblsFromTuple(value, 3, s) < 0 ) return NULL ;
      nextarg++ ;
    } else if ( PyList_Check(value) ) {
      if( dblsFromList(value, 3, s) < 0 ) return NULL ;
      nextarg++ ;
    } else {
      if( nextarg + 2 >= narg ) {
	PyErr_SetString(PyExc_ValueError,
	 "atomPut requires at least sx, sy, sz after name, isub, xyz") ;
	return NULL ;
      }
      for( i=0 ; i<3 ; i++ ) {
	if ( (value = PyTuple_GetItem(args, nextarg+i)) == NULL ) return NULL ;
	if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	  PyErr_SetString(PyExc_ValueError,
			  "sx sy sz args to atomPut must be numeric") ;
	  return NULL ;
	}
	s[i] = PyFloat_AsDouble(value) ;
      }
      nextarg += 3 ;
    }
  }
    
  m = atomDefPtr->m ;
  if( nextarg < narg ) {
    if ( (value = PyTuple_GetItem(args, nextarg)) == NULL ) return NULL ;
    if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError,
		      "moment arg to atomPut must be numeric") ;
      return NULL ;
    }
    m = PyFloat_AsDouble(value) ;
    nextarg++ ;
  }
  c = 1. ;
  if( nextarg < narg ) {
    if ( (value = PyTuple_GetItem(args, nextarg)) == NULL ) return NULL ;
    if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError,
		      "occupation arg to atomPut must be numeric") ;
      return NULL ;
    }
    c = PyFloat_AsDouble(value) ;
  }

  if( (atomgroup = nextAtomGroup(&atomslist)) == NULL ) return NULL ;
  atomgroup->i = atomslist.n ;
  if( ! sym ) isub = -1 ;
  if( putAtom(name, atomgroup, atomDefPtr, isub, r, s, m, c) < 0 ) return NULL;
  //calcAtomAllQ(1, atomgroup, qlist) ;
  return Py_BuildValue("") ;
}

static int etyp = 0 ;
// editAtom will now be declared before methods
//static PyObject *editAtom(PyObject *self, PyObject *args) ;

static PyObject *
edit(PyObject *self, PyObject *args)
{
  PyObject *elst ;
  etyp = 0 ;
  elst = editAtom(self, args) ;
  return elst ;
}

static PyObject *
editPos(PyObject *self, PyObject *args)
{
  PyObject *elst ;
  etyp = 1 ;
  elst = editAtom(self, args) ;
  return elst ;
}
static PyObject *
editDir(PyObject *self, PyObject *args)
{
  PyObject *elst ;
  etyp = 2 ;
  elst = editAtom(self, args) ;
  return elst ;
}
static PyObject *
editMom(PyObject *self, PyObject *args)
{
  PyObject *elst ;
  etyp = 3 ;
  elst = editAtom(self, args) ;
  return elst ;
}
static PyObject *
editOcc(PyObject *self, PyObject *args)
{
  PyObject *elst ;
  etyp = 4 ;
  elst = editAtom(self, args) ;
  return elst ;
}
static PyObject *
editMag(PyObject *self, PyObject *args)
{
  PyObject *elst ;
  etyp = 5 ;
  elst = editAtom(self, args) ;
  return elst ;
}
static PyObject *
editFF(PyObject *self, PyObject *args)
{
  PyObject *elst ;
  etyp = 6 ;
  elst = editAtom(self, args) ;
  return elst ;
}
static PyObject *
editDW(PyObject *self, PyObject *args)
{
  PyObject *elst ;
  etyp = 7 ;
  elst = editAtom(self, args) ;
  return elst ;
}



static int sameAtomPos(myATOM *a, myATOM *b, double tol) ;

static PyObject *editByPos(PyObject *self, PyObject *args)
{
  /*
    editByPos(x y z [sx sy sz] [m] [c])
    NB x y z sx sy sz are in reduced coordinates
    NB all atoms in a group have the same pt sym and so should
    have same DW and FF so we leave DW and FF as part of atomDef
    If you have more than one site symm with a given atom type
    you can copy the def and make a new atomDef to edit
    a new FF or DW etc
  */

  PyObject *value ;
  PyObject *atomList ;

  int narg, i, j, k ;
  double toler, smag ;

  ATOMgroup *atomgroup, *atomgroupi ;
  myATOM *atom ;
  myATOM test ;

  narg = PyTuple_Size(args) ;

  if( narg < 3 ) {
    PyErr_SetString(PyExc_ValueError,
		    "must specify atom pos for editByPos") ;
    return NULL ;
  }
  for( i=0 ; i<3 ; i++ ) {
    if ( (value = PyTuple_GetItem(args, i)) == NULL ) return NULL ;
    if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError,
		      "editByPos first 3 args must be numeric coords") ;
      return NULL ;
    }
    test.r[i] = PyFloat_AsDouble(value) ;
  }

  /* find the atom with this position */

  toler = pow(10., (double)(*atomtolerance)) ;

  atomgroup = atomslist.atomslist ;
  for( i=0 ; i<atomslist.n ; i++ ) {
    atomgroupi = atomgroup + i ;
    if( atomgroupi->atom == NULL || atomgroupi->atoms == NULL ) continue ;
    for( j=0 ; j<atomgroupi->natoms ; j++ ) {
      atom = atomgroupi->atoms + j ;
      if( ! sameAtomPos( &test, atom, toler ) ) continue ;
      /* found same position so edit and return */
      if( narg >= 6 ) {
	for( j=0 ; j<3 ; j++ ) {
	  if ( (value = PyTuple_GetItem(args, j+3)) == NULL ) return NULL ;
	  if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	    PyErr_SetString(PyExc_ValueError,
			    "editByPos args 3-6 must be red numeric spin dir");
	    return NULL ;
	  }
	  atom->s[j] = PyFloat_AsDouble(value) ;
	}
	smag = 0. ;
	for( k=0 ; k<3 ; k++ ) {
	  atom->ss[k] =
	    atom->s[0]*auni[k] + atom->s[1]*buni[k] + atom->s[2]*cuni[k] ;
	  smag += atom->ss[k]*atom->ss[k] ;
	}
	if( smag > 0. ) {
	  smag = sqrt(smag) ;
	  for( k=0 ; k<3 ; k++ ) atom->ss[k] /= smag ;
	}
      }
      if( narg >= 7 ) {
	if ( (value = PyTuple_GetItem(args, 6)) == NULL ) return NULL ;
	if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	  PyErr_SetString(PyExc_ValueError,
			  "editByPos arg 7 must be numeric moment");
	  return NULL ;
	}
	atom->m = PyFloat_AsDouble(value) ;
      }
      if( narg >= 8 ) {
	if ( (value = PyTuple_GetItem(args, 7)) == NULL ) return NULL ;
	if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	  PyErr_SetString(PyExc_ValueError,
			  "editByPos arg 8 must be numeric occupancy");
	  return NULL ;
	}
	atom->c = PyFloat_AsDouble(value) ;
      }

      if ( (atomList = MakeAtomList(atomgroupi->atom->name, atomgroupi->isub,
				     atom)) == NULL )
	return NULL ;
      return atomList ;
    }
  }
  return Py_BuildValue("s", "NO atom at that position") ;
}

static int omitAtoms() ;

static PyObject *
atomDel(PyObject *self, PyObject *args)
{
  /*
    atomDel( name OR symmetrySymbol OR index or [indexStart, indexEnd])
    deletions are marked first so that any combination of indices and other
    args can be used
  */

  PyObject *arg, *atom, *value, *item ;
  PyObject *atomsList, *groupList, *atomList ;

  int indx1, indx2 ;
  int narg, i, j, nc, isub, nlist, ndel ;
  double r[3], s[3], m, c ;
  char *name, *sym ;

  ATOMdef *atomDefPtr ;
  ATOMSlist *atomslst ;
  ATOMgroup *atomgroup, *atomgroupi ;
  myATOM *atoms ;

  for( i=0 ; i<atomslist.n ; i++ ) atomslist.atomslist[i].omit = 0 ;

  narg = PyTuple_Size(args) ;
  ndel = 0 ;
  for( i=0 ; i<narg ; i++ ) {
    sym = NULL ;
    name = NULL ;
    indx1 = -1 ;
    if ( (arg = PyTuple_GetItem(args, i)) == NULL ) return NULL ;
    if ( PyString_Check(arg) ) {
      name = PyString_AsString(arg) ;
      if( name == NULL || (nc = strlen(name)) < 1 ) continue ;
      if( nc == 1 && islower(name[0]) ) sym = name ;
    } else if ( PyInt_Check(arg) ) {
      indx1 = PyInt_AsLong(arg) ;
      indx2 = indx1 ;
    } else if ( PyList_Check(arg) ) {
      nlist = PyList_Size(arg) ;
      if( nlist < 2 ) {
	PyErr_SetString(PyExc_ValueError,
			"del atom list requires start and end indices") ;
	return NULL ;
      }
      if ( (value = PyList_GetItem(arg, 0)) == NULL ) return NULL ;
      if ( ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
			"del atom list requires start and end int indices") ;
	return NULL ;
      }
      indx1 = PyInt_AsLong(value) ;
      if ( (value = PyList_GetItem(arg, 1)) == NULL ) return NULL ;
      if ( ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
			"del atom list requires start and end int indices") ;
	return NULL ;
      }
      indx2 = PyInt_AsLong(value) ;
    } else if ( PyTuple_Check(arg) ) {
      nlist = PyTuple_Size(arg) ;
      if( nlist < 2 ) {
	PyErr_SetString(PyExc_ValueError,
			"del atom list requires start and end indices") ;
	return NULL ;
      }
      if ( (value = PyTuple_GetItem(arg, 0)) == NULL ) return NULL ;
      if ( ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
			"del atom tuple requires start and end int indices") ;
	return NULL ;
      }
      indx1 = PyInt_AsLong(value) ;
      if ( (value = PyTuple_GetItem(arg, 1)) == NULL ) return NULL ;
      if ( ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
			"del atom tuple requires start and end int indices") ;
	return NULL ;
      }
      indx2 = PyInt_AsLong(value) ;
    }
    /* now go thru the atomlist and mark to omit */

    if( sym || name ) {
      for( i=0 ; i<atomslist.n ; i++ ) {
	if( sym && atomslist.atomslist[i].Wsym == sym[0] ) {
	  if( atomslist.atomslist[i].omit != 1 ) {
	    atomslist.atomslist[i].omit = 1 ;
	    ndel++ ;
	  }
	} else if( name && atomslist.atomslist[i].atom &&
		   strcmp(name,atomslist.atomslist[i].atom->name)==0 ) {
	  if( atomslist.atomslist[i].omit != 1 ) {
	    atomslist.atomslist[i].omit = 1 ;
	    ndel++ ;
	  }
	}
      }
    } else if( indx1 > 0 && indx1 <= atomslist.n ) {
      if( indx2 > atomslist.n ) indx2 = atomslist.n ;
      for( j=indx1-1 ; j<indx2 ; j++ ) {
	if( atomslist.atomslist[j].omit != 1 ) {
	  atomslist.atomslist[j].omit = 1 ;
	  ndel++ ;
	}
      }
    }
  }

  if( ndel < 1 ) return Py_BuildValue("s", "NO atoms deleted") ;

  omitAtoms() ;
  return Py_BuildValue("is", ndel, " atoms deleted") ;
}

static int pruneAtoms() ;

static PyObject *
atomPrune(PyObject *self, PyObject *args)
{
  int ndel ;
  ndel = pruneAtoms() ;
  return Py_BuildValue("is", ndel, " atomgroups deleted") ;
}

static PyObject *
atomList(PyObject *self, PyObject *args)
{
  /*
    atomList( name/atomicnum, isub/symmSymbol )
  */

  PyObject *arg, *atom, *value, *item ;
  PyObject *atomsList, *groupList, *atomLIST ;

  int narg, i, j, nc, isub, nmiss, atomicnum, istat, indx, nlst, nlist ;
  int indx1, indx2 ;
  double r[3], s[3], m, c ;
  char *name, *sym ;

  ATOMdef *atomDefPtr ;
  ATOMSlist *atomslst ;
  ATOMgroup *atomgroup, *atomgroupi ;
  myATOM *atoms ;

  narg = PyTuple_Size(args) ;

  atomicnum = 0 ;
  name = NULL ;
  sym = NULL ;
  isub = -1 ;


  if ( narg > 0 ) {
    for( i=0 ; i<atomslist.n ; i++ ) atomslist.atomslist[i].omit = 1 ;
    nlst = 0 ;
  } else {
    for( i=0 ; i<atomslist.n ; i++ ) atomslist.atomslist[i].omit = 0 ;
    nlst = atomslist.n ;
  }

  for( i=0 ; i<narg ; i++ ) {
    sym = NULL ;
    name = NULL ;
    indx1 = -1 ;
    if ( (arg = PyTuple_GetItem(args, i)) == NULL ) return NULL ;
    if ( PyString_Check(arg) ) {
      name = PyString_AsString(arg) ;
      if( name == NULL || (nc = strlen(name)) < 1 ) continue ;
      if( nc == 1 && islower(name[0]) ) sym = name ;
    } else if ( PyInt_Check(arg) ) {
      indx1 = PyInt_AsLong(arg) ;
      indx2 = indx1 ;
    } else if ( PyList_Check(arg) ) {
      nlist = PyList_Size(arg) ;
      if( nlist < 2 ) {
	PyErr_SetString(PyExc_ValueError,
			"atom list requires start and end indices") ;
	return NULL ;
      }
      if ( (value = PyList_GetItem(arg, 0)) == NULL ) return NULL ;
      if ( ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
			"atom list requires start and end int indices") ;
	return NULL ;
      }
      indx1 = PyInt_AsLong(value) ;
      if ( (value = PyList_GetItem(arg, 1)) == NULL ) return NULL ;
      if ( ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
			"atom list requires start and end int indices") ;
	return NULL ;
      }
      indx2 = PyInt_AsLong(value) ;
    }
    /* now go thru the atomlist and add to list those chosen */

    if( sym || name ) {
      for( i=0 ; i<atomslist.n ; i++ ) {
	if( sym && atomslist.atomslist[i].Wsym == sym[0] ) {
	  if( atomslist.atomslist[i].omit != 0 ) {
	    atomslist.atomslist[i].omit = 0 ;
	    nlst++ ;
	  }
	} else if( name && atomslist.atomslist[i].atom &&
		   strcmp(name,atomslist.atomslist[i].atom->name)==0 ) {
	  if( atomslist.atomslist[i].omit != 0 ) {
	    atomslist.atomslist[i].omit = 0 ;
	    nlst++ ;
	  }
	}
      }
    } else if( indx1 > 0 && indx1 <= atomslist.n ) {
      for( j=indx1-1 ; j<indx2 ; j++ ) {
	if( atomslist.atomslist[i].omit != 0 ) {
	  atomslist.atomslist[j].omit = 0 ;
	  nlst++ ;
	}
      }
    }
  }

  //if( nlst < 1 ) return Py_BuildValue("s", "NO atoms to list") ;

  /*
    return list of atoms from atomlist
    matching any supplied name/atomicnum or isub/sym
    make each atom list contents readable by put command
  */
  if ( (atomsList = PyList_New(nlst)) == NULL ) return NULL ;
  atomgroup = atomslist.atomslist ;
  nlst = 0 ;
  for( i=0 ; i<atomslist.n ; i++ ) {
    if( atomslist.atomslist[i].omit == 1 ) continue ;
    atomgroupi = atomgroup + i ;
    //if( atomgroupi->atom == NULL ) continue ;
    //if( atomicnum > 0 && atomgroupi->atom->atomic != atomicnum ) continue ;
    //if( name != NULL && strcmp(name, atomgroupi->atom->name) ) continue ;
    //if( sym != NULL && sym[0] != atomgroupi->Wsym ) continue ;
    //if( isub >= 0 && isub != atomgroupi->isub ) continue ;

    if ( (groupList = PyList_New(atomgroupi->natoms)) == NULL ) {
      Py_DECREF(atomsList) ; return NULL ;
    }
    atoms = atomgroupi->atoms ;
    for( j=0 ; j<atomgroupi->natoms ; j++ ) {

      if ( (atomLIST = MakeAtomList(atomgroupi->atom->name, atomgroupi->isub,
				     atoms+j)) == NULL )
	{ Py_DECREF(groupList) ; Py_DECREF(atomsList) ; return NULL ; }
      PyList_SET_ITEM(groupList, j, atomLIST) ;
    }
    PyList_SET_ITEM(atomsList, nlst, groupList) ;
    nlst++ ;
  }
  return atomsList ;
}

static PyObject *
atomPrint(PyObject *self, PyObject *args)
{
  /*
    atomPrint( name or symmSymbol or index or [istart,iend] )
    printf version of atomList
  */

  PyObject *arg, *atom, *value, *item ;
  PyObject *atomsList, *groupList, *atomLIST ;

  int narg, i, j, isub, nmiss, atomicnum, istat, indx, nlst, nlist, nc ;
  int indx1, indx2 ;
  double r[3], s[3], m, c ;
  char *name, *sym ;

  ATOMdef *atomDefPtr ;
  ATOMSlist *atomslst ;
  ATOMgroup *atomgroup, *atomgroupi ;
  myATOM *atoms ;

  narg = PyTuple_Size(args) ;

  atomicnum = 0 ;
  name = NULL ;
  sym = NULL ;
  isub = -1 ;


  if ( narg > 0 ) {
    for( i=0 ; i<atomslist.n ; i++ ) atomslist.atomslist[i].omit = 1 ;
    nlst = 0 ;
  } else {
    for( i=0 ; i<atomslist.n ; i++ ) atomslist.atomslist[i].omit = 0 ;
    nlst = atomslist.n ;
  }

  for( i=0 ; i<narg ; i++ ) {
    sym = NULL ;
    name = NULL ;
    indx1 = -1 ;
    if ( (arg = PyTuple_GetItem(args, i)) == NULL ) return NULL ;
    if ( PyString_Check(arg) ) {
      name = PyString_AsString(arg) ;
      if( name == NULL || (nc = strlen(name)) < 1 ) continue ;
      if( nc == 1 && islower(name[0]) ) sym = name ;
    } else if ( PyInt_Check(arg) ) {
      indx1 = PyInt_AsLong(arg) ;
      indx2 = indx1 ;
    } else if ( PyList_Check(arg) ) {
      nlist = PyList_Size(arg) ;
      if( nlist < 2 ) {
	PyErr_SetString(PyExc_ValueError,
			"atom list requires start and end indices") ;
	return NULL ;
      }
      if ( (value = PyList_GetItem(arg, 0)) == NULL ) return NULL ;
      if ( ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
			"atom list requires start and end int indices") ;
	return NULL ;
      }
      indx1 = PyInt_AsLong(value) ;
      if ( (value = PyList_GetItem(arg, 1)) == NULL ) return NULL ;
      if ( ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
			"atom list requires start and end int indices") ;
	return NULL ;
      }
      indx2 = PyInt_AsLong(value) ;
    }
    /* now go thru the atomlist and add to list those chosen */

    if( sym || name ) {
      for( i=0 ; i<atomslist.n ; i++ ) {
	if( sym && atomslist.atomslist[i].Wsym == sym[0] ) {
	  if( atomslist.atomslist[i].omit != 0 ) {
	    atomslist.atomslist[i].omit = 0 ;
	    nlst++ ;
	  }
	} else if( name && atomslist.atomslist[i].atom &&
		   strcmp(name,atomslist.atomslist[i].atom->name)==0 ) {
	  if( atomslist.atomslist[i].omit != 0 ) {
	    atomslist.atomslist[i].omit = 0 ;
	    nlst++ ;
	  }
	}
      }
    } else if( indx1 > 0 && indx1 <= atomslist.n ) {
      for( j=indx1-1 ; j<indx2 ; j++ ) {
	if( atomslist.atomslist[i].omit != 0 ) {
	  atomslist.atomslist[j].omit = 0 ;
	  nlst++ ;
	}
      }
    }
  }

  if( nlst < 1 ) return Py_BuildValue("s", "NO atoms to list") ;

  /*
    return list of atoms from atomlist
  */

  atomgroup = atomslist.atomslist ;
  for( i=0 ; i<atomslist.n ; i++ ) {
    if( atomslist.atomslist[i].omit == 1 ) continue ;
    atomgroupi = atomgroup + i ;
    printf("atomGroup# %i  atom=%s isub=%i natoms=%i\n",
	   i+1, atomgroupi->atom->name,atomgroupi->isub,atomgroupi->natoms) ;

    atoms = atomgroupi->atoms ;
    for( j=0 ; j<atomgroupi->natoms ; j++ ) {
      printf("%9g %9g %9g %9g %9g %9g %9g %9g\n",
	     atoms[j].r[0],atoms[j].r[1],atoms[j].r[2],
	     atoms[j].s[0],atoms[j].s[1],atoms[j].s[2],
	     atoms[j].m,atoms[j].c) ;
    }
  }
  return Py_BuildValue("") ;
}


static PyObject *MakeFlagList(int n, Flag *flst) ;
static int ParseFlagsFromArgs(int nflags, Flag *flags,
			      int narg, PyObject *args) ;

static PyObject *
atomPutFlags(PyObject *self, PyObject *args)
{
  int narg ;

  narg = PyTuple_Size(args) ;
 
  if ( narg < 1 ) {
    /* just return flag list, name/intvalue pairs */
    return MakeFlagList(nputflags, putflags) ;
  }
  if ( ! ParseFlagsFromArgs(nputflags, putflags, narg, args) )
    return NULL ;
  return MakeFlagList(nputflags, putflags) ;
}


/* calc COMMANDS */

static PyObject *
calcFlags(PyObject *self, PyObject *args)
{
  int narg ;

  narg = PyTuple_Size(args) ;
 
  if ( narg < 1 ) {
    /* just return flag list, name/intvalue pairs */
    return MakeFlagList(ncalcflags, calcflags) ;
  }
  if ( ! ParseFlagsFromArgs(ncalcflags, calcflags, narg, args) )
    return NULL ;
  return MakeFlagList(ncalcflags, calcflags) ;
}

static PyObject *
lorentzFlags(PyObject *self, PyObject *args)
{
  int narg ;

  narg = PyTuple_Size(args) ;
 
  if ( narg < 1 ) {
    /* just return flag list, name/intvalue pairs */
    return MakeFlagList(nlorentzflags, lorentzflags) ;
  }
  if ( ! ParseFlagsFromArgs(nlorentzflags, lorentzflags, narg, args) )
    return NULL ;
  return MakeFlagList(nlorentzflags, lorentzflags) ;
}

static PyObject *
calcColumns(PyObject *self, PyObject *args)
{
  int narg ;

  narg = PyTuple_Size(args) ;
 
  if ( narg < 1 ) {
    /* just return flag list, name/intvalue pairs */
    return MakeFlagList(ncalccolumns, calccolumns) ;
  }
  if ( ! ParseFlagsFromArgs(ncalccolumns, calccolumns, narg, args) )
    return NULL ;
  return MakeFlagList(ncalccolumns, calccolumns) ;
}

static PyObject *
unknownColumns(PyObject *self, PyObject *args)
{
  int narg ;

  narg = PyTuple_Size(args) ;
 
  if ( narg < 1 ) {
    /* just return flag list, name/intvalue pairs */
    return MakeFlagList(nunknowncolumns, unknowncolumns) ;
  }
  if ( ! ParseFlagsFromArgs(nunknowncolumns, unknowncolumns, narg, args) )
    return NULL ;
  return MakeFlagList(nunknowncolumns, unknowncolumns) ;
}

static int writeCalc(FILE *fpt) ;
static int outColumnsList( int *ilst ) ;
static void addOut( Qs *qs ) ;
static void calcQ( ATOMSlist *atomslst, Qs *qs ) ;
static void sortQmag() ;
static void docalc() ;
static FILE *openFileFromArgs(PyObject *args, char *fstr) ;


static PyObject *
calcFile(PyObject *self, PyObject *args)
{
  int narg ;
  PyObject *newargs ;
  FILE *fpt ;

  narg = PyTuple_Size(args) ;
  newargs = PyTuple_GetSlice(args, 0, narg) ;
  if( (fpt = openFileFromArgs(newargs, "calcFile")) == NULL ) {
    Py_DECREF(newargs) ;return NULL ;
  }
  Py_DECREF(newargs) ;
  fprintf(fpt, "%s\n", progid) ;

  docalc() ;

  if( ! writeCalc(fpt) ) { fclose(fpt) ; return NULL ; }
  fclose(fpt) ;
  return Py_BuildValue("") ;
}

static int writeSetup(FILE *fpt) ;
static PyObject *
setupFile(PyObject *self, PyObject *args)
{
  FILE *fpt ;

  if( (fpt = openFileFromArgs(args, "setupFile")) == NULL ) return NULL ;
  fprintf(fpt, "%s\n", progid) ;
  if( ! writeSetup(fpt) ) { fclose(fpt) ; return NULL ; }
  fclose(fpt) ;
  return Py_BuildValue("") ;
}

static int writeAtoms(FILE *fpt) ;
static PyObject *
atomFile(PyObject *self, PyObject *args)
{
  FILE *fpt ;

  if( (fpt = openFileFromArgs(args, "atomFile")) == NULL ) return NULL ;
  fprintf(fpt, "%s\n", progid) ;
  if( ! writeAtoms(fpt) ) { fclose(fpt) ; return NULL ; }
  fclose(fpt) ;
  return Py_BuildValue("") ;
}

static int writeQ(FILE *fpt) ;
static PyObject *
QFile(PyObject *self, PyObject *args)
{
  FILE *fpt ;

  if( (fpt = openFileFromArgs(args, "QFile")) == NULL ) return NULL ;
  fprintf(fpt, "%s\n", progid) ;
  if( ! writeQ(fpt) ) { fclose(fpt) ; return NULL ; }
  fclose(fpt) ;
  return Py_BuildValue("") ;
}

static int writeFlags(FILE *fpt) ;
static PyObject *
saveFlags(PyObject *self, PyObject *args)
{
  FILE *fpt ;

  if( (fpt = fopen(flagfile, "w")) == NULL ) {
    PyErr_SetString(PyExc_ValueError, "saveFlags failed to open file");
    return NULL ;
  }
  fprintf(fpt, "%s\n", progid) ;
  if( ! writeFlags(fpt) ) { fclose(fpt) ; return NULL ; }
  fclose(fpt) ;
  return Py_BuildValue("") ;
}

static PyObject *
allFlags(PyObject *self, PyObject *args)
{
  if( ! writeFlags(NULL) ) return NULL ;
  return Py_BuildValue("") ;
}



static PyObject *
calc(PyObject *self, PyObject *args)
{
  /*
    output select columns using calcColumns command
    columns: h k l Q Qoutofplane CSp CSm SQ++ SQ-- SQ+- SQ-+ CStotal 2thet omeg
  */

  PyObject *value, *ptrList, *lblList, *QList ;

  int narg, i, j ;
  int ncols ;
  int cols[16] ;
  char dptr[16] ;


  /*
   * pass any args through to calcFlags
   */

  narg = PyTuple_Size(args) ;
  if( narg > 0 ) if( calcFlags(self, args) == NULL ) return NULL ;

  docalc() ;

  /*
   * the rest of the work is figuring out what to return and print
   * main choice is to look at print flag
   * if print flag TRUE then printf the output else return Python Lists
   * secondary choice is whether to return pointers to the selected columns
   * or the actual data
   */

  if( *calcprint ) {
    if( ! writeCalc(NULL) ) return NULL ;
    return Py_BuildValue("") ;
  }

  ncols = outColumnsList(cols) ;

  if( *calcpointers ) {
    if ( (ptrList = PyList_New(ncols)) == NULL ) return NULL ;
    for( i=0 ; i<ncols ; i++ ) {
      sprintf(dptr, "%p", *(outcols[cols[i]].d)) ;
      PyList_SET_ITEM(ptrList, i, Py_BuildValue("s", dptr)) ;
    }
  } else if( *perQ ) {
    /* must return list of lists, one list per column or one list per Q */
    if ( (ptrList = PyList_New(NQ)) == NULL ) return NULL ;
    for( i=0 ; i<NQ ; i++ ) {
      if ( (QList = PyList_New(ncols)) == NULL ) {
	Py_DECREF(ptrList) ; return NULL ;
      }
      for( j=0 ; j<ncols ; j++ )
	PyList_SET_ITEM(QList, j,
			Py_BuildValue("d", (*(outcols[cols[j]].d))[i])) ;

      PyList_SET_ITEM(ptrList, i, QList) ;
    }
  } else {
    if ( (ptrList = PyList_New(ncols)) == NULL ) return NULL ;
    for( j=0 ; j<ncols ; j++ ) {
      if ( (QList = PyList_New(NQ)) == NULL ) {
	Py_DECREF(ptrList) ; return NULL ;
      }
      for( i=0 ; i<NQ ; i++ )
	PyList_SET_ITEM(QList, i,
			Py_BuildValue("d", (*(outcols[cols[j]].d))[i])) ;
      PyList_SET_ITEM(ptrList, j, QList) ;
    }
  }

  if ( (lblList = PyList_New(ncols)) == NULL ) {
    Py_DECREF(ptrList) ; return NULL ;
  }
  for( j=0 ; j<ncols ; j++ )
    PyList_SET_ITEM(lblList, j, Py_BuildValue("s", outcols[cols[j]].l)) ;

  return Py_BuildValue("iNN", NQ, ptrList, lblList) ;
}

/* Q commands */

static PyObject *QaddIndex(int indx, PyObject *self, PyObject *args) ;

static PyObject *
Qadd(PyObject *self, PyObject *args)
{
  PyObject *result ;
  /*
    args: h k l hs ks ls N
    or lists [h, k, l, {hs, ks, ls, N}], ...
    or       [[h,k,l,{}],...]
  */

  result = QaddIndex(0, self, args) ;
  return result ;
}
static PyObject *
Qi(PyObject *self, PyObject *args)
{
  /*
    first arg is an integer index in qlist
    other args: h k l hs ks ls N
    or lists [h, k, l, {hs, ks, ls, N}], ...
    or       [[h,k,l,{}],...]
  */

  int narg, indx ;
  PyObject *result, *value, *newargs ;

  narg = PyTuple_Size(args) ;
  if( narg < 1 ) {
    result = QaddIndex(0, self, args) ;
    return result ;
  }

  if ( (value = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
  if ( ! PyInt_Check(value) ) {
    PyErr_SetString(PyExc_ValueError, "first arg to Qi must be int index") ;
    return NULL ;
  }
  indx = PyInt_AsLong(value) ;
  /* now put the remaining args into a tuple and call QaddIndex */
  newargs = PyTuple_GetSlice(args, 1, narg) ;
  result = QaddIndex(indx, self, newargs) ;
  Py_DECREF(newargs) ;
  return result ;
}
static PyObject *
Qindices(PyObject *self, PyObject *args)
{
  /*
    make multiplex Qlists from hkl ranges
    args: [hi hf {hs}] {[ki kf {ks}] {[li lf {ls}]}}
  */

  int i, j, k, narg, nr, im, i1, i2 ;
  int N[3], Nt[3] ;
  double HKL[3][3], Q[3], QS[3], temp ;
  PyObject *arg, *result, *value, *newt ;

  narg = PyTuple_Size(args) ;
  if( narg < 1 ) {
    result = QaddIndex(0, self, args) ;
    return result ;
  }

  if ( (arg = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
  if ( PyList_Check(arg) ) {
    if ( (nr = PyList_Size(arg)) < 2 ) {
      PyErr_SetString(PyExc_ValueError,
  "Qindices list arg requires at least start and end for hrange [hi hf {hs}]");
      return NULL ;
    }
    HKL[0][2] = 1. ;
    for( i=0 ; i<nr ; i++ ) {
      if ( (value = PyList_GetItem(arg, i)) == NULL ) return NULL ;
      if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
	    "Qindices hrange list must be numeric values [hi hf {hs}]");
	return NULL ;
      }
      HKL[0][i] = PyFloat_AsDouble(value) ;
    }
  } else if( PyFloat_Check(arg) || PyInt_Check(arg) ) {
    HKL[0][0] = HKL[0][1] = PyFloat_AsDouble(arg) ;
    HKL[0][2] = 0. ;
  } else {
    PyErr_SetString(PyExc_ValueError,
		    "Qindices single values must be numeric");
    return NULL ;
  }

  HKL[1][0] = HKL[0][0] ;
  HKL[1][1] = HKL[0][1] ;
  HKL[1][2] = HKL[0][2] ;
  if( narg > 1 ) {
    if ( (arg = PyTuple_GetItem(args, 1)) == NULL ) return NULL ;
    if ( PyList_Check(arg) ) {
      if ( (nr = PyList_Size(arg)) < 2 ) {
	PyErr_SetString(PyExc_ValueError,
  "Qindices range list requires start and end for krange [ki kf {ks}]");
	return NULL ;
      }
      HKL[1][2] = 1. ;
      for( i=0 ; i<nr ; i++ ) {
	if ( (value = PyList_GetItem(arg, i)) == NULL ) return NULL ;
	if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	  PyErr_SetString(PyExc_ValueError,
	      "Qindices krange list must be numeric values [ki kf {ks}]");
	  return NULL ;
	}
	HKL[1][i] = PyFloat_AsDouble(value) ;
      }
    } else if( PyFloat_Check(arg) || PyInt_Check(arg) ) {
      HKL[1][0] = HKL[1][1] = PyFloat_AsDouble(arg) ;
      HKL[1][2] = 0. ;
    } else {
      PyErr_SetString(PyExc_ValueError,
		      "Qindices single values must be numeric");
      return NULL ;
    }
  }

  HKL[2][0] = HKL[1][0] ;
  HKL[2][1] = HKL[1][1] ;
  HKL[2][2] = HKL[1][2] ;
  if( narg > 2 ) {
    if ( (arg = PyTuple_GetItem(args, 2)) == NULL ) return NULL ;
    if ( PyList_Check(arg) ) {
      if ( (nr = PyList_Size(arg)) < 2 ) {
	PyErr_SetString(PyExc_ValueError,
     "Qindices range list requires start and end for lrange [li lf {ls}]");
	return NULL ;
      }
      HKL[2][2] = 1. ;
      for( i=0 ; i<nr ; i++ ) {
	if ( (value = PyList_GetItem(arg, i)) == NULL ) return NULL ;
	if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	  PyErr_SetString(PyExc_ValueError,
	    "Qindices lrange list must be numeric values [li lf {ls}]");
	  return NULL ;
	}
	HKL[2][i] = PyFloat_AsDouble(value) ;
      }
    } else if( PyFloat_Check(arg) || PyInt_Check(arg) ) {
      HKL[2][0] = HKL[2][1] = PyFloat_AsDouble(arg) ;
      HKL[2][2] = 0. ;
    } else {
      PyErr_SetString(PyExc_ValueError,
		      "Qindices single values must be numeric");
      return NULL ;
    }
  }

  N[0] = N[1] = N[2] = 1 ;
  for( i=0 ; i<3 ; i++ ) {
    if( HKL[i][1] < HKL[i][0] ) {
      temp = HKL[i][1] ;
      HKL[i][1] = HKL[i][0] ;
      HKL[i][0] = temp ;
    }
    if( HKL[i][2] < 0 ) HKL[i][2] = -HKL[i][2] ;
    if( HKL[i][2] > 0. ) N[i] = (int)((HKL[i][1] - HKL[i][0])/HKL[i][2]) + 1 ;
  }

  im = 0 ; i1 = 1 ; i2 = 2 ;
  if( N[1] > N[0] && N[1] >= N[2] ) { im = 1 ; i1 = 0 ; i2 = 2 ; }
  if( N[2] > N[0] && N[2] > N[1] ) { im = 2 ; i1 = 0 ; i2 = 1 ; }

  Q[im] = HKL[im][0] ;
  QS[im] = HKL[im][2] ;
  QS[i1] = QS[i2] = 0. ;
  Nt[im] = N[im] ;
  Nt[i1] = Nt[i2] = 0 ;

  /* make the im qlist for each i1 i2 index */
  for( i=0 ; i<N[i1] ; i++ ) {
    Q[i1] = HKL[i1][0] + i*HKL[i1][2] ;
    for( j=0 ; j<N[i2] ; j++ ) {
      Q[i2] = HKL[i2][0] + j*HKL[i2][2] ;
      /* setup the tuple for Qadd */
      if( (newt = PyTuple_New(9)) == NULL ) return NULL ;
      for( k=0 ; k<3 ; k++ ) {
	PyTuple_SetItem(newt, k, PyFloat_FromDouble(Q[k])) ;
	PyTuple_SetItem(newt, k+3, PyFloat_FromDouble(QS[k])) ;
	PyTuple_SetItem(newt, k+6, PyInt_FromLong(Nt[k])) ;
      }
      if( (result = QaddIndex(0, self, newt)) == NULL ) {
	Py_DECREF(newt) ; return NULL ;
      }
      Py_DECREF(newt) ;
      Py_DECREF(result) ;
    }
  }
  return Py_BuildValue("") ;
}


static PyObject *
Qdelete(PyObject *self, PyObject *args)
{
  /*
    args are integer index in qlist or start, end
  */

  int i, j, narg, istart, iend, nnew ;
  PyObject *value ;

  narg = PyTuple_Size(args) ;
  if( narg < 1 ) {
    PyErr_SetString(PyExc_ValueError,
		    "args to Qdel must be int index or istart and iend");
    return NULL ;
  }
  if( qlist.n < 1 ) {
    PyErr_SetString(PyExc_ValueError, "empty qlist");
    return NULL ;
  }

  if ( (value = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
  if ( ! PyInt_Check(value) ) {
    PyErr_SetString(PyExc_ValueError,
		    "args to Qdel must be int index or istart and iend");
    return NULL ;
  }
  istart = PyInt_AsLong(value) ;
  iend = istart ;
  if ( narg > 1 ) {
    if ( (value = PyTuple_GetItem(args, 1)) == NULL ) return NULL ;
    if ( ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError,
		      "args to Qdel must be int index or istart and iend");
      return NULL ;
    }
    iend = PyInt_AsLong(value) ;
  }

  if( istart < 1 || istart > qlist.n ) {
    PyErr_SetString(PyExc_ValueError,
		    "Qdel index start must be 1 to Nqlist");
    return NULL ;
  }
  if( iend > 0 && iend < istart ) iend = istart ;
  if( iend <= 0 || iend > qlist.n ) iend = qlist.n ;

  nnew = qlist.n - (iend - istart + 1) ;

  /* free the Qs storage for the ones del */
  for( j=istart ; j<=iend ; j++ ) {
    free(qlist.qlist[j-1].q) ;
    qlist.qlist[j-1].q = NULL ;
  }
  /* copy remaining up to freed index positions */
  for( i=iend, j=istart-1 ; i<qlist.n ; i++, j++ )
    *(qlist.qlist+j) = *(qlist.qlist+i) ;
  /* set qs to null for zombies at end of list */
  for( i=nnew ; i<qlist.n ; i++ ) {
    qlist.qlist[i].n = 0 ;
    qlist.qlist[i].q = NULL ;
  }
  qlist.n = nnew ;
  return Py_BuildValue("") ;
}

static void sortQlist(Qlist *qlst) ;

static PyObject *
Qsort(PyObject *self, PyObject *args)
{
  sortQlist(&qlist) ;
  return Py_BuildValue("") ;
}

static PyObject *printQlist(int indx, PyObject *Qlst) ;

static PyObject *
Qlisting(PyObject *self, PyObject *args)
{
  int narg, indx1, indx2 ;
  int i, nq ;
  static char Qheading[] =
 "Qindex#       H          K          L        Hstep      Kstep      Lstep     N"
    ;
  PyObject *value, *values, *newargs, *result ;

  narg = PyTuple_Size(args) ;
  indx1 = 0 ;

  if( narg < 1 ) {
    printf("%s\n", Qheading) ;
    if ( (result = QaddIndex(0, self, args)) == NULL ) return NULL ;
    nq = PyList_Size(result) ;
    for( i=0 ; i<nq ; i++ ) {
      if ( (value = PyList_GetItem(result, i)) == NULL ) {
	Py_DECREF(result) ; return NULL ;
      }
      if ( ! printQlist(i+1, value) ) { Py_DECREF(result) ; return NULL ; }
    }
    Py_DECREF(result) ;
    return Py_BuildValue("") ;
  }

  if ( (value = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
  if ( ! PyInt_Check(value) ) {
    PyErr_SetString(PyExc_ValueError,
		    "first arg to Qlisting must be int index") ;
  }
  indx1 = PyInt_AsLong(value) ;
  indx2 = indx1 ;

  if ( narg > 1 ) {
    if ( (value = PyTuple_GetItem(args, 1)) == NULL ) return NULL ;
    if ( ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError,
		      "second arg to Qlisting must be final int index") ;
    }
    indx2 = PyInt_AsLong(value) ;
  }

  if( indx1 < 1 || indx1 > qlist.n ) {
    PyErr_SetString(PyExc_ValueError,
		    "Qlisting index start must be 1 to Nqlist");
    return NULL ;
  }
  if( indx2 > 0 && indx2 < indx1 ) indx2 = indx1 ;
  if( indx2 <= 0 || indx2 > qlist.n ) indx2 = qlist.n ;

  printf("%s\n", Qheading) ;
  for( i=indx1 ; i<=indx2 ; i++ ) {
    if( (newargs = PyTuple_New(0)) == NULL ) return NULL ;
    if ( (result = QaddIndex(i, self, newargs)) == NULL ) {
      Py_DECREF(newargs) ; return NULL ;
    }
    Py_DECREF(newargs) ;
    /* result should be list with a single sublist */
    if ( (value = PyList_GetItem(result, 0)) == NULL ) {
      Py_DECREF(result) ;return NULL ;
    }
    if( ! printQlist(i, value) ) {
      Py_DECREF(result) ; return NULL ;
    }
    Py_DECREF(result) ;
  }
  return Py_BuildValue("") ;
}

static int readStructureFile(char *fil, char *typ, int ph) ;

static PyObject *
readFile(PyObject *self, PyObject *args)
{
  int narg, ph ;
  PyObject *value ;
  char *fil, *typ ;

  typ = NULL ;
  narg = PyTuple_Size(args) ;
  if( narg < 1 ) {
    PyErr_SetString(PyExc_ValueError, "readFile requires filename arg") ;
    return NULL ;
  }
  if ( (value = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
  if ( ! PyString_Check(value) ) {
    PyErr_SetString(PyExc_ValueError, "readFile requires filename arg") ;
    return NULL ;
  }
  fil = PyString_AsString(value) ;
  ph = -1 ;
  if( narg > 1 ) {
    if ( (value = PyTuple_GetItem(args, 1)) == NULL ) return NULL ;
    if ( PyString_Check(value) ) {
      typ = PyString_AsString(value) ;
    } else if( PyInt_Check(value) ) {
      ph = PyInt_AsLong(value) ;
    } else {
      PyErr_SetString(PyExc_ValueError,
	 "readFile second arg must be string file type or int phase index") ;
      return NULL ;
    }
  }

  if( narg > 2 ) {
    if ( (value = PyTuple_GetItem(args, 2)) == NULL ) return NULL ;
    if ( ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError,
		      "readFile third arg must be int phase index") ;
      return NULL ;
    }
    ph = PyInt_AsLong(value) ;
  }
  if( readStructureFile(fil, typ, ph) < 0 ) return NULL ;
  return Py_BuildValue("") ;
}

/*****************************************************************/
/* auxilliary functions */
/*****************************************************************/


static void initQ(Q *qpt) ;
static int nextQ(Qlist *qlst) ;

static PyObject *QaddIndex(int indx, PyObject *self, PyObject *args)
{
  /*
    args: h k l hs ks ls N
    or lists [h, k, l, {hs, ks, ls, N}], ...
    or       [[h,k,l,{}],...]
  */

  PyObject *qsList, *qList ;
  PyObject *arg, *item, *value ;
  PyObject *newt ;

  int narg, i, j, nq, nql ;
  int istart, iend ;
  Qs *qs ;

  narg = PyTuple_Size(args) ;
  if( narg < 1 ) {
    /* return a list of the Q lists */

    qs = qlist.qlist ;
    if( indx < 1 || indx > qlist.n ) {
      istart = 0 ;
      iend = qlist.n ;
    } else {
      istart = indx - 1 ;
      iend = indx ;
    }
    nq = iend - istart ;
    if ( ( qsList = PyList_New(nq) ) == NULL ) return NULL ;
    nq = 0 ;
    for( i=istart ; i<iend ; i++ ) {
      nql = 3 ;
      if ( qs[i].n > 0 ) nql = 7 ;
      if ( ( qList = PyList_New(nql) ) == NULL ) {
	Py_DECREF(qsList) ; return NULL ;
      }
      PyList_SET_ITEM(qList, 0, Py_BuildValue("d", qs[i].hkl[0])) ;
      PyList_SET_ITEM(qList, 1, Py_BuildValue("d", qs[i].hkl[1])) ;
      PyList_SET_ITEM(qList, 2, Py_BuildValue("d", qs[i].hkl[2])) ;
      if ( qs[i].n > 0 ) {
	PyList_SET_ITEM(qList, 3, Py_BuildValue("d", qs[i].step[0])) ;
	PyList_SET_ITEM(qList, 4, Py_BuildValue("d", qs[i].step[1])) ;
	PyList_SET_ITEM(qList, 5, Py_BuildValue("d", qs[i].step[2])) ;
	PyList_SET_ITEM(qList, 6, Py_BuildValue("i", qs[i].n)) ;
      }
      PyList_SET_ITEM(qsList, nq, qList) ;
      nq++ ;
    }
    return qsList ;
  }

  if ( (arg = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
  if ( PyList_Check(arg) ) {
    /* is this a list of lists, or each arg is a list */
    if ( (item = PyList_GetItem(arg, 0)) == NULL ) return NULL ;
    if ( PyList_Check(item) ) {
      /*
       * list of lists
       * so for each list in the list call QaddIndex recursively
       */
      narg = PyList_Size(arg) ;
      for( i=0 ; i<narg ; i++ ) {
	if ( (item = PyList_GetItem(arg, i)) == NULL ) return NULL ;
	if ( ! PyList_Check(item) ) {
	  PyErr_SetString(PyExc_ValueError,
			  "arg to Qadd should have beem list") ;
	  return NULL ;
	}
	newt = PyList_AsTuple(item) ;
	if ( ! PyTuple_Check(newt) ) {
	  PyErr_SetString(PyExc_ValueError,
			  "failed to make tuple from list in Qadd") ;
	  Py_DECREF(newt) ;
	  return NULL ;
	}
	if ( QaddIndex(indx, self, newt) == NULL ) {
	  Py_DECREF(newt) ; return NULL ;
	}
	Py_DECREF(newt) ;
      }
    } else {
      /*
       * args are each lists
       * elements of the List should be standard args to Q
       */
      for( i=0 ; i<narg ; i++ ) {
	if ( (item = PyTuple_GetItem(args, i)) == NULL ) return NULL ;
	if ( ! PyList_Check(item) ) {
	  PyErr_SetString(PyExc_ValueError,
			  "arg to Qadd should have beem list") ;
	  return NULL ;
	}
	newt = PyList_AsTuple(item) ;
	if ( ! PyTuple_Check(newt) ) {
	  PyErr_SetString(PyExc_ValueError,
			  "failed to make tuple from list in Qadd") ;
	  Py_DECREF(newt) ;
	  return NULL ;
	}
	if ( QaddIndex(indx, self, newt) == NULL ) {
	  Py_DECREF(newt) ;
	  return NULL ;
	}
	Py_DECREF(newt) ;
      }
    }
    return Py_BuildValue("") ;
  }

  /*
    get here if args are from tuple just doubles defining a q or qrange
    there must be at least 3 (h k l)
  */
  if( narg < 3 ) {
    PyErr_SetString(PyExc_ValueError, "Qadd requires at least 3 args h, k, l");
    return NULL ;
  }

  if( indx < 1 || indx > qlist.n ) {
    if( (indx = nextQ(&qlist)) < 0 ) return NULL ;
  } else {
    indx-- ;
  }
  qs = qlist.qlist + indx ;

  qs->n = 1 ;
  for( i=0 ; i<3 ; i++ ) {
    if ( (value = PyTuple_GetItem(args, i)) == NULL ) return NULL ;
    if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError, "h k l args to Qadd must be numeric") ;
      return NULL ;
    }
    qs->hkl[i] = PyFloat_AsDouble(value) ;
    qs->step[i] = 0. ;
  }
  for( i=3 ; i<(narg<=6?narg:6) ; i++ ) {
    if ( (value = PyTuple_GetItem(args, i)) == NULL ) return NULL ;
    if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError, "step args to Qadd must be numeric") ;
      return NULL ;
    }
    qs->step[i-3] = PyFloat_AsDouble(value) ;
  }
  if( narg >= 7 ) {
    if ( (value = PyTuple_GetItem(args, 6)) == NULL ) return NULL ;
    if ( ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError, "npts arg to Qadd must be int") ;
      return NULL ;
    }
    qs->n = PyInt_AsLong(value) ;
  }

  if( qs->n < 1 ) qs->n = 1 ;

  qs->q = (Q *)realloc(qs->q, qs->n*sizeof(Q)) ;

  for( i=0 ; i<qs->n ; i++ ) {
    for( j=0 ; j<3 ; j++ ) qs->q[i].hkl[j] =
			     qs->hkl[j] + (double)i*qs->step[j] ;
    initQ(qs->q + i) ;
  }
  return Py_BuildValue("") ;
}

static PyObject *printQlist(int indx, PyObject *Qlst)
{
  int j, Nq ;
  double hkl ;
  PyObject *values ;
  printf("%7d", indx) ;
  /* each of the sublists has the seven values: h k l hs ks ls N */
  for( j=0 ; j<6 ; j++ ) {
    if ( (values = PyList_GetItem(Qlst, j)) == NULL ) return NULL ;
    hkl = PyFloat_AsDouble(values) ;
    printf(" %10g", hkl) ;
  }
  if ( (values = PyList_GetItem(Qlst, 6)) == NULL ) return NULL ;
  Nq = PyInt_AsLong(values) ;
  printf(" %7d\n", Nq) ;
  return Qlst ;
}


static int SymLookup(char S)
{
  int i ;
  Spcgrp *spc ;
  spc = spcgrps + curSpaceGrpN ;
  for( i=0 ; i<spc->nsub ; i++ ) {
    if( spc->sub[i].Wyck == S ) return i ;
  }
  return 0 ;
}

static int putAtom(char *id, ATOMgroup *atomgroup,
		   ATOMdef *atomDefPtr, int isub,
		   double r[3], double s[3], double m, double c)
{
  static char location[] = "putAtom" ;
  int i, j, k, l, ntrans ;
  int natoms, newatom ;
  int Na, Nb, Nc, itemp ;
  int a1, a2, b1, b2, c1, c2, ii[3], or[3] ;
  double smag, diff, magdiff ;
  double atomDiffTolerance ;
  DCMPLX b ;

  Spcgrp *spcgrp ;
  Subgrp *subgrp ;
  Subpt *subpt ;
  myATOM *atom, *atomp ;

  atomgroup->done = 0 ;
  if( atomgroup->name != NULL ) free(atomgroup->name) ;
  if( id == NULL ) if( ! setString(&(atomgroup->name),atomDefPtr->name) )
    return -1 ;
  else if( ! setString(&(atomgroup->name), id) ) return -1 ;
  atomgroup->atom = atomDefPtr ;
  atomgroup->spcgrp = spcgrps[curSpaceGrpN].spcnum ;
  atomgroup->setting = spcgrps[curSpaceGrpN].setting[0] ;

  if( isub < 0 && *usespacegroup ) isub = 0 ;
  /* recall first subgroup is most general positions */

  spcgrp = NULL ;
  subgrp = NULL ;
  ntrans = 0 ;

  atomDiffTolerance = pow(10., (double)(*atomtolerance)) ;

  if( *usespacegroup ) {
    /*
      have to generate all atoms for this symmetry
      the site multiplicity will be ntrans * npts in the subgrp
    */
    if( curSpaceGrpN < 0 || curSpaceGrpN >= NspaceGrpsLoaded ) return ;
    spcgrp = spcgrps + curSpaceGrpN ;
    subgrp = spcgrp->sub + isub ;
    ntrans = 1 ;
    if( spcgrp->ntrans > 0 ) ntrans = spcgrp->ntrans ;
    atomgroup->mult = ntrans*subgrp->npts ;
    if( *multicell ) {
      if( *amax < *amin ) { itemp = *amin ; *amin = *amax ; *amax = itemp ; }
      if( *bmax < *bmin ) { itemp = *bmin ; *bmin = *bmax ; *bmax = itemp ; }
      if( *cmax < *cmin ) { itemp = *cmin ; *cmin = *cmax ; *cmax = itemp ; }
      Na = (*amax - *amin) + 1 ;
      Nb = (*bmax - *bmin) + 1 ;
      Nc = (*cmax - *cmin) + 1 ;
      atomgroup->mult *= Na*Nb*Nc ;
      a1 = *amin ;
      a2 = *amax ;
      b1 = *bmin ;
      b2 = *bmax ;
      c1 = *cmin ;
      c2 = *cmax ;
    } else {
      a1 = a2 = b1 = b2 = c1 = c2 = 0 ;
    }

    atomgroup->Wsym = subgrp->Wyck ;
    atomgroup->ptsym = subgrp->sitesym ;
    atomgroup->isub = isub ;
    //printf("putAtom spcgrp#=%d isub=%d symbol=%c mult=%d\n",
    //   spcgrps[curSpaceGrpN].spcnum, isub, subgrp->Wyck, atomgroup->mult) ;
  } else {
    atomgroup->mult = 1 ;
    atomgroup->natoms = 1 ;
  }
  atomgroup->atoms =
    (myATOM*)realloc(atomgroup->atoms, atomgroup->mult * sizeof(myATOM)) ;
  if ( memExc(isNULL(atomgroup->atoms), location) ) return -1 ;

  /* compute b for current ki */
  //b = bresonance(ki, atomDefPtr->b, atomDefPtr->Eres, atomDefPtr->HW) ;
  /* try GSAS */
  b = bGSAS(atomDefPtr->b, Ei, atomDefPtr->nbcoef, atomDefPtr->bcoef) ;

  if( ! *usespacegroup ) {
    atom = atomgroup->atoms ;
    smag = 0. ;
    for( k=0 ; k<3 ; k++ ) {
      atom->r[k] = r[k] ;
      atom->s[k] = s[k] ;
      atom->ss[k] = s[0]*auni[k] + s[1]*buni[k] + s[2]*cuni[k] ;
      smag += atom->ss[k]*atom->ss[k] ;
    }
    if( smag > 0. ) {
      smag = sqrt(smag) ;
      for( k=0 ; k<3 ; k++ ) atom->ss[k] /= smag ;
    } else {
      m = 0. ;
    }
    atom->c = c ;
    atom->m = m ;
    atom->b = b ;
    if( *autoprune ) pruneAtoms() ;
    return 1 ;
  }

  atom = atomgroup->atoms ;
  natoms = 0 ;
  //printf("putAtom input coords xyz=%g %g %g\n", r[0],r[1],r[2]) ;

  /* calculate the originating atom cell indices */
  for( i=0 ; i<3 ; i++ ) or[i] = (int)r[i] ;

  for( ii[0] = a1 ; ii[0] <= a2 ; ii[0]++ ) {
    for( ii[1] = b1 ; ii[1] <= b2 ; ii[1]++ ) {
      for( ii[2] = c1 ; ii[2] <= c2 ; ii[2]++ ) {

	for( i=0 ; i<ntrans ; i++ ) {
	  //if( spcgrp->ntrans > 0 ) printf("putAtom for trans=%g %g %g\n",
	  //			    spcgrp->trans[i].t[0],
	  //			    spcgrp->trans[i].t[1],
	  //			    spcgrp->trans[i].t[2]) ;
	  for( j=0 ; j<subgrp->npts ; j++ ) {
	    subpt = subgrp->pts + j ;
	    atom->c = c ;
	    atom->m = m ;
	    atom->b = b ;
	    smag = 0. ;
	    //printf(" coefs and constant for x y z\n") ;
	    for( k=0 ; k<3 ; k++ ) {
	      if( spcgrp->ntrans > 0 ) atom->r[k] = spcgrp->trans[i].t[k] ;
	      else atom->r[k] = 0. ;
	      atom->s[k] = s[k] ;
	      atom->ss[k] = s[0]*auni[k] + s[1]*buni[k] + s[2]*cuni[k] ;
	      smag += atom->ss[k]*atom->ss[k] ;
	      
	      for( l=0 ; l<3 ; l++ ) atom->r[k] += subpt->m[k][l]*r[l] ;
	      atom->r[k] += subpt->t[k] ;
	      /* first make all coordinates in initial unit cell */
	      while( atom->r[k] < (double)or[k] ) atom->r[k] += 1. ;
	      while( atom->r[k] >= (double)or[k]+1.) atom->r[k] -= 1. ;
	      //printf(" %g %g %g %g\n",
	      //  subpt->m[k][0],subpt->m[k][1],subpt->m[k][2],subpt->t[k]) ;
	      /* Now do multicell indices if requested */
	      atom->r[k] += ii[k] ;
	    }
	    /* check that this atom position has not already been added */
	    newatom = 1 ;
	    for( k=0 ; k<natoms ; k++ ) {
	      atomp = atomgroup->atoms + k ;
	      magdiff = 0. ;
	      for( l=0 ; l<3 ; l++ ) {
		diff = atomp->r[l] - atom->r[l] ;
		magdiff += diff*diff ;
	      }
	      if( magdiff > 0. ) magdiff = sqrt(magdiff) ;
	      if( magdiff < atomDiffTolerance ) {
		newatom = 0 ;
		break ;
	      }
	    }
	    if( ! newatom ) continue ;
	    ++natoms ;
	    if( smag > 0. ) {
	      smag = sqrt(smag) ;
	      for( k=0 ; k<3 ; k++ ) atom->ss[k] /= smag ;
	    }
	    atom++ ;
	  }
	}
      }
    }
  }
  /* realloc storage in case some atoms were rejected */
  atomgroup->natoms = natoms ;
  if( atomgroup->mult != natoms ) {
    atomgroup->atoms =
      (myATOM*)realloc(atomgroup->atoms, atomgroup->natoms * sizeof(myATOM)) ;
    if ( memExc(isNULL(atomgroup->atoms), location) ) return -1 ;
  }
  if( *autoprune ) pruneAtoms() ;

  return 1 ;
}

static PyObject *editAtom(PyObject *self, PyObject *args)
{
  /*
    NB because this functions is called by all registered editAtom calls
    it must return a new reference!
    I had thought that because it wasnt a registered method I could change the
    standard arg list to have etyp the first arg.
    But this causes neg ref on cleanup after the function call (I think)
    So I'm going to make this a registered function and etype a global
    variable.
    atomEdit( name/symmSymbol/index,  [x y z] [sx sy sz] [m] [c] DW FF)
    NB x y z sx sy sz are in reduced coordinates
    etyp=0 enter all x y z sx sy sz m c  >= 3 args
    etyp=1           x y z               >= 3 args
    etyp=2           sx sy sz            >= 3 args
    etyp=3           m                   >= 1 arg
    etyp=4           c                   >= 1 arg
    etyp=5           [m0x qxh qxk qxl phx,r0x,r0y,r0z],[m0y...],[m0z,...]
    etyp=6           atomdefName OR S,L,[coefs]
    etyp=7           atomdefName OR DWcoefs
    NB all atoms in a group have the same pt sym and so should
    have same DW and FF so typically DW and FF set from atomDef
    If you have more than one site symm with a given atom type
    you can copy the def and make a new atomDef to edit
    a new FF or DW etc
    etyp=5 defines spin direction and moment from modulation
    mx(r) = m0x cos(qx.(r-r0) + phx)
    my(r) = m0y cos(qy.(r-r0) + phy)
    mz(r) = m0z cos(qz.(r-r0) + phz)
  */

  PyObject *arg, *atom, *value, *gname, *item, *newt, *lst ;
  PyObject *rsrc, *ssrc ;
  PyObject *atomsList, *groupList, *atomLIST ;
  PyObject *newargs, *j2frac, *Cval, *ffval, *result ;

  int narg, nl, i, j, k, l, isub, nmiss, istat, next, nnew, nret ;
  int roff, soff, moff, ooff ;
  int indx1, indx2 ;
  double r[3], s[3], m, c, cosarg ;
  double m0[3], ph[3], dr[3], mxyz[3] ;
  static double qmx[3], qmy[3], qmz[3] ;
  static double r0x[3], r0y[3], r0z[3], lastr0[3] ;
  static double *qmhkl[3] = { qmx, qmy, qmz } ;
  static double *r0[3] = { r0x, r0y, r0z } ;
  char *name, *sym ;

  ATOMdef *atomDefPtr ;
  ATOMSlist *atomslst ;
  ATOMgroup *atomgroup, *atomgroupi ;
  myATOM *atoms, *atomi ;

  narg = PyTuple_Size(args) ;

  name = NULL ;
  sym = NULL ;
  isub = -1 ;
  indx1 = indx2 = 0 ;

  if( narg < 1 ) {
    PyErr_SetString(PyExc_ValueError,
		    "must specify atom name, symmSymb or index for atomEdit") ;
    return NULL ;
  }
  if( etyp < 3 && narg < 2 ) {
    PyErr_SetString(PyExc_ValueError,
	       "atomEdit Pos or Dir requires atom-name and a list or 3 flts") ;
    return NULL ;
  }
  if( etyp >= 3 && narg < 2 ) {
    PyErr_SetString(PyExc_ValueError,
	       "atomEdit Mom Occ Mag FF DW requires atom-name and 1 arg min") ;
    return NULL ;
  }


  if ( (atom = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;
  if ( PyString_Check(atom) ) {
    name = PyString_AsString(atom) ;
    if ( strlen(name) == 1 && islower(name[0]) ) { sym = name ; name = NULL ; }
  } else if ( PyInt_Check(atom) ) {
    indx1 = PyInt_AsLong(atom) ;
    indx2 = indx1 ;
  } else if ( PyTuple_Check(atom) ) {
    if( PyTuple_Size(atom) < 2 ) {
      PyErr_SetString(PyExc_ValueError,
		  "editAtom first arg tuple values are index range min max") ;
      return NULL ;
    }
    if ( (value = PyTuple_GetItem(atom, 0)) == NULL ) return NULL ;
    indx1 = PyInt_AsLong(value) ;
    if ( (value = PyTuple_GetItem(atom, 1)) == NULL ) return NULL ;
    indx2 = PyInt_AsLong(value) ;
    if ( indx2 < indx1 ) { nmiss = indx1 ; indx1 = indx2 ; indx2 = nmiss ; }
  } else if ( PyList_Check(atom) ) {
    /*
      recursively call editAtom foreach member of the list
    */
    if( (newargs = PyTuple_New(narg)) == NULL ) return NULL ;
    /*
      first copy the args 1 to narg-1
      NB for some reason PyTuple_GetSlice doesnt copy tuple
      contents but appears to borrow refs
    */
    for( i=1 ; i<narg ; i++ ) {
      if ( (value = PyTuple_GetItem(args, i)) == NULL ) {
	Py_DECREF(newargs) ; return NULL ;
      }
      PyTuple_SetItem(newargs, i, Py_BuildValue("O",value)) ;
    }
    nl = PyList_Size(atom) ;
    for( i=0 ; i<nl ; i++ ) {
      if ( (value = Py_BuildValue("O", PyList_GetItem(atom, i))) == NULL ) {
	Py_DECREF(newargs) ; return NULL ;
      }
      PyTuple_SetItem(newargs, 0, value) ;
      if ( (result = editAtom(NULL, newargs)) == NULL ) {
	Py_DECREF(newargs) ; return NULL ;
      }
      Py_DECREF(result) ;
    }
    Py_DECREF(newargs) ;
    return Py_BuildValue("s", "editAtom list processed OK") ;
  } else {
    PyErr_SetString(PyExc_ValueError, "invalid first arg to atomEdit") ;
    return NULL ;
  }

  if( sym != NULL ) isub = SymLookup(sym[0]) ;

  /* return a list of the atoms edited */

  atomgroup = atomslist.atomslist ;
  nret = 0 ;
  for( i=0 ; i<atomslist.n ; i++ ) {
    atomgroupi = atomgroup + i ;
    if( atomgroupi->atom == NULL ) continue ;
    if( indx1 > 0 && (indx2 < i+1 || indx1 > i+1) ) continue ;
    if( name != NULL && strcmp(name, atomgroupi->atom->name) ) continue ;
    if( sym != NULL && sym[0] != atomgroupi->Wsym ) continue ;
    if( isub >= 0 && isub != atomgroupi->isub ) continue ;
    nret++ ;
  }

  if ( (atomsList = PyList_New(nret)) == NULL ) return NULL ;
  nret = 0 ;
  for( i=0 ; i<atomslist.n ; i++ ) {
    atomgroupi = atomgroup + i ;
    if( atomgroupi->atom == NULL ) continue ;
    if( indx1 > 0 && (indx2 < i+1 || indx1 > i+1) ) continue ;
    if( name != NULL && strcmp(name, atomgroupi->atom->name) ) continue ;
    if( sym != NULL && sym[0] != atomgroupi->Wsym ) continue ;
    if( isub >= 0 && isub != atomgroupi->isub ) continue ;

    atoms = atomgroupi->atoms ;
    /* default values from atoms[0] */
    for( j=0 ; j<3 ; j++ ) {
      r[j] = atoms[0].r[j] ;
      s[j] = atoms[0].s[j] ;
    }
    m = atoms[0].m ;
    c = atoms[0].c ;

    if( etyp < 3 ) {
      roff = soff = moff = ooff = -1 ;
      if ( (value = PyTuple_GetItem(args, 1)) == NULL ) {
	Py_DECREF(atomsList) ; return NULL ;
      }
      if ( PyList_Check(value) ) {
	/* arg1 list */
	rsrc = PyList_AsTuple(value) ;
	if( PyTuple_Size(rsrc) < 3 ) {
	  PyErr_SetString(PyExc_ValueError, "coord list must have 3 floats") ;
	  Py_DECREF(rsrc) ; Py_DECREF(atomsList) ;
	  return NULL ;
	}
	roff = 0 ;
	if( narg > 2 ) {
	  if ( (value = PyTuple_GetItem(args, 2)) == NULL ) {
	    Py_DECREF(rsrc) ; Py_DECREF(atomsList) ;
	    return NULL ;
	  }
	  if ( PyList_Check(value) ) {
	    /* arg2=momdir also list */
	    ssrc = PyList_AsTuple(value) ;
	    if( PyTuple_Size(ssrc) < 3 ) {
	      PyErr_SetString(PyExc_ValueError,
			      "coord list must have 3 floats") ;
	      Py_DECREF(rsrc) ; Py_DECREF(ssrc) ; Py_DECREF(atomsList) ;
	      return NULL ;
	    }
	    soff = 0 ;
	    if( narg > 3 ) moff = 3 ;
	    if( narg > 4 ) ooff = 4 ;
	  } else {
	    /* arg1 list arg2-4 floats dir */
	    ssrc = PyTuple_New(3) ;
	    for( j=0 ; j<3 ; j++ ) {
	      if ( j+2 < narg ) {
		if ( (value = PyTuple_GetItem(args, j+2)) == NULL ) {
		  Py_DECREF(rsrc) ; Py_DECREF(ssrc) ; Py_DECREF(atomsList) ;
		  return NULL ;
		}
		if( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
		  Py_DECREF(rsrc) ; Py_DECREF(ssrc) ; Py_DECREF(atomsList) ;
		  PyErr_SetString(PyExc_ValueError, "momdir must be numeric") ;
		  return NULL ;
		}
		PyTuple_SetItem(ssrc, j,
				Py_BuildValue("d", PyFloat_AsDouble(value))) ;
	      } else {
		PyTuple_SetItem(ssrc, j, Py_BuildValue("d", 0.)) ;
	      }
	    }
	    soff = 2 ;
	    if( narg > 5 ) moff = 5 ;
	    if( narg > 6 ) ooff = 6 ;
	  }
	}
      } else {
	/* args1-3 float NB no isub as for putAtom */
	rsrc = PyTuple_New(3) ;
	for( j=0 ; j<3 ; j++ ) {
	  if ( j+1 < narg ) {
	    if ( (value = PyTuple_GetItem(args, j+1)) == NULL ) {
	      Py_DECREF(rsrc) ; Py_DECREF(atomsList) ;
	      return NULL ;
	    }
	    if( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	      Py_DECREF(rsrc) ; Py_DECREF(atomsList) ;
	      PyErr_SetString(PyExc_ValueError, "pos must be numeric") ;
	      return NULL ;
	    }
	    PyTuple_SetItem(rsrc, j,
			    Py_BuildValue("d", PyFloat_AsDouble(value))) ;
	  } else {
	    PyTuple_SetItem(rsrc, j, Py_BuildValue("d", 0.)) ;
	  }
	}
	roff = 1 ;
	if( narg > 4 ) {
	  if ( (value = PyTuple_GetItem(args, 4)) == NULL ) return NULL ;
	  if ( PyList_Check(value) ) {
	    /* args1-3 float arg 4 list */
	    ssrc = PyList_AsTuple(value) ;
	    if( PyTuple_Size(ssrc) < 3 ) {
	      PyErr_SetString(PyExc_ValueError,
			      "coord list must have 3 floats") ;
	      Py_DECREF(atomsList) ;
	      Py_DECREF(rsrc) ; Py_DECREF(ssrc) ;
	      return NULL ;
	    }
	    soff = 0 ;
	    if( narg > 5 ) moff = 5 ;
	    if( narg > 6 ) ooff = 6 ;
	  } else {
	    /* args1-3 float args4-6 float */
	    ssrc = PyTuple_New(3) ;
	    for( j=0 ; j<3 ; j++ ) {
	      if ( j+4 < narg ) {
		if ( (value = PyTuple_GetItem(args, j+4)) == NULL ) {
		  Py_DECREF(rsrc) ; Py_DECREF(ssrc) ; Py_DECREF(atomsList) ;
		  return NULL ;
		}
		if( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
		  Py_DECREF(rsrc) ; Py_DECREF(ssrc) ; Py_DECREF(atomsList) ;
		  PyErr_SetString(PyExc_ValueError, "momdir must be numeric") ;
		  return NULL ;
		}
		PyTuple_SetItem(ssrc, j,
				Py_BuildValue("d", PyFloat_AsDouble(value))) ;
	      } else {
		PyTuple_SetItem(ssrc, j, Py_BuildValue("d", 0.)) ;
	      }
	    }
	    soff = 4 ;
	    if( narg > 7 ) moff = 7 ;
	    if( narg > 6 ) ooff = 8 ;
	  }
	}
      }

      for( j=0 ; j<3 ; j++ ) {
	if ( (value = PyTuple_GetItem(rsrc, j)) == NULL ) {
	  Py_DECREF(atomsList) ;
	  Py_DECREF(rsrc) ; if( soff >= 0 ) Py_DECREF(ssrc) ;
	  return NULL ;
	}
	if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	  PyErr_SetString(PyExc_ValueError,
		     "args to atomEdit must be numeric or list of numeric") ;
	  Py_DECREF(atomsList) ;
	  Py_DECREF(rsrc) ; if( soff >= 0 ) Py_DECREF(ssrc) ;
	  return NULL ;
	}
	if( etyp == 2 ) s[j] = PyFloat_AsDouble(value) ;
	else r[j] = PyFloat_AsDouble(value) ;
      }
      Py_DECREF(rsrc) ;

      if( etyp == 0 && soff >= 0 ) {
	for( j=0 ; j<3 ; j++ ) {
	  if ( (value = PyTuple_GetItem(ssrc, j)) == NULL ) {
	    Py_DECREF(atomsList) ;
	    Py_DECREF(ssrc) ;
	    return NULL ;
	  }
	  if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	    PyErr_SetString(PyExc_ValueError,
		     "args to atomEdit must be numeric or lists of numeric") ;
	    Py_DECREF(atomsList) ;
	    Py_DECREF(ssrc) ;
	    return NULL ;
	  }
	  s[j] = PyFloat_AsDouble(value) ;
	}
	Py_DECREF(ssrc) ;
      }
      if( etyp == 0 && moff >= 0 ) {
	if ( (value = PyTuple_GetItem(args, moff)) == NULL ) {
	  Py_DECREF(atomsList) ;
	  return NULL ;
	}
	if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	  PyErr_SetString(PyExc_ValueError,
			  "moment arg to atomEdit must be numeric") ;
	  Py_DECREF(atomsList) ;
	  return NULL ;
	}
	m = PyFloat_AsDouble(value) ;
      }
      if( etyp == 0 && ooff >= 0 ) {
	if ( (value = PyTuple_GetItem(args, ooff)) == NULL ) {
	  Py_DECREF(atomsList) ;
	  return NULL ;
	}
	if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	  PyErr_SetString(PyExc_ValueError,
			  "occupancy arg to atomEdit must be numeric") ;
	  Py_DECREF(atomsList) ;
	  return NULL ;
	}
	c = PyFloat_AsDouble(value) ;
      }

    } else if( etyp == 3 ) {
      if ( (value = PyTuple_GetItem(args, 1)) == NULL ) {
	Py_DECREF(atomsList) ;
	return NULL ;
      }
      if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
			"moment arg to atomEdit must be numeric") ;
	Py_DECREF(atomsList) ;
	return NULL ;
      }
      m = PyFloat_AsDouble(value) ;
      for( j=0 ; j<atomgroupi->natoms ; j++ ) {
	atomi = atomgroupi->atoms + j ;
	atomi->m = m ;
      }
    } else if( etyp == 4 ) {
      if ( (value = PyTuple_GetItem(args, 1)) == NULL ) {
	Py_DECREF(atomsList) ;
	return NULL ;
      }
      if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	PyErr_SetString(PyExc_ValueError,
			"occupancy arg to atomEdit must be numeric") ;
	Py_DECREF(atomsList) ;
	return NULL ;
      }
      c = PyFloat_AsDouble(value) ;
    } else if( etyp == 5 ) {
      for( j=0 ; j<3 ; j++ ) {
	m0[j] = 0. ;
	for( k=0 ; k<3 ; k++ ) { qmhkl[j][k] = 0. ; r0[j][k] = 0. ; }
	ph[j] = 0. ;
	if ( j+1 < narg ) {
	  if ( (lst = PyTuple_GetItem(args, j+1)) == NULL ) return NULL ;
	  if ( ! PyList_Check(lst) ) {
	    PyErr_SetString(PyExc_ValueError,
	"editMag args to atomEdit after first must be lists [m h k l phase x y z] OR [m (h k l) phase (x y z)]") ;
	    Py_DECREF(atomsList) ;
	    return NULL ;
	  }

	  nnew = PyList_Size(lst) ;
	  if( nnew > 0 ) {
	  /* m h k l phase OR m (h k l) phase */
	    if ( (value = PyList_GetItem(lst, 0)) == NULL ) {
	      Py_DECREF(atomsList) ;
	      return NULL ;
	    }
	    if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	      PyErr_SetString(PyExc_ValueError,
			      "m0 in editMag list must be numeric") ;
	      Py_DECREF(atomsList) ;
	      return NULL ;
	    }
	    m0[j] = PyFloat_AsDouble(value) ;
	  }
	  next = 0 ;
	  if( nnew > 1 ) {
	    if ( (value = PyList_GetItem(lst, 1)) == NULL ) {
	      Py_DECREF(atomsList) ;
	      return NULL ;
	    }
	    if ( PyTuple_Check(value) ) {
	      if( dblsFromTuple(value, 3, qmhkl[j]) < 0 ) {
		Py_DECREF(atomsList) ;
		return NULL ;
	      }
	      if( nnew > 2 ) next = 2 ;
	    } else if( nnew > 3 ) {
	      for( k=0 ; k<3 ; k++ ) {
		if ((value = PyList_GetItem(lst, k+1)) == NULL) {
		  Py_DECREF(atomsList) ;
		  return NULL ;
		}
		if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
		  PyErr_SetString(PyExc_ValueError,
			     "qmag values in editMag list must be numeric") ;
		  Py_DECREF(atomsList) ;
		  return NULL ;
		}
		qmhkl[j][k] = PyFloat_AsDouble(value) ;
	      }
	      if( nnew > 4 ) next = 4 ;
	    }
	  }
	  if( next > 0 ) {
	    if ( (value = PyList_GetItem(lst, next)) == NULL ) {
	      Py_DECREF(atomsList) ;
	      return NULL ;
	    }
	    if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	      PyErr_SetString(PyExc_ValueError,
			      "phase values in editMag list must be numeric") ;
	      Py_DECREF(atomsList) ;
	      return NULL ;
	    }
	    ph[j] = PyFloat_AsDouble(value) ;
	    next++ ;
	    if( next >= nnew ) next = 0 ;
	  }
	  if( next > 0 ) {
	    if ( (value = PyList_GetItem(lst, next)) == NULL ) {
	      Py_DECREF(atomsList) ;
	      return NULL ;
	    }
	    if ( PyTuple_Check(value) ) {
	      if( dblsFromTuple(value, 3, r0[j]) < 0 ) {
		Py_DECREF(atomsList) ;
		return NULL ;
	      }
	    } else if( next+2 < nnew ) {
	      for( k=0 ; k<3 ; k++ ) {
		if ((value=PyList_GetItem(lst, k+next)) == NULL) {
		  Py_DECREF(atomsList) ;
		  return NULL;
		}
		if ( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
		  PyErr_SetString(PyExc_ValueError,
			     "r0 values in editMag list must be numeric") ;
		  Py_DECREF(atomsList) ;
		  return NULL ;
		}
		r0[j][k] = PyFloat_AsDouble(value) ;
	      }
	    }
	  }
	}
      }

      /* now foreach atom in this atomgroup edit m sx sy sz */
      for( j=0 ; j<atomgroupi->natoms ; j++ ) {
	/*
	  use m0 qmhkl ph and r0 to compute the components along
	  the direct lattice vectors
	*/
	atomi = atomgroupi->atoms + j ;
	for( k=0 ; k<3 ; k++ ) {
	  for( l=0 ; l<3 ; l++ ) dr[l] = atomi->r[l] - r0[k][l] ;
	  cosarg = TWOPI * dotpro(qmhkl[k], dr) + ph[k] ;
	  s[k] = m0[k]*cos(cosarg) ;
	}
	/* compute the total moment and s in sample Cartesian ss */
	m = 0. ;
	for( k=0 ; k<3 ; k++ ) {
	  mxyz[k] = s[0]*auni[k] + s[1]*buni[k] + s[2]*cuni[k] ;
	  m += mxyz[k]*mxyz[k] ;
	}
	if( m > 0. ) {
	  m = sqrt(m) ;
	  for( k=0 ; k<3 ; k++ ) mxyz[k] /= m ;
	}
	atomi->m = m ;
	for( k=0 ; k<3 ; k++ ) {
	  atomi->s[k] = s[k] ;
	  atomi->ss[k] = mxyz[k] ;
	}
      }
    } else if( etyp == 6 ) {
      /* edit FF */
      if ( (value = PyTuple_GetItem(args, 1)) == NULL ) {
	Py_DECREF(atomsList) ;
	return NULL ;
      }
      if ( PyString_Check(value) ) {
	/* srcname */
	if( (newargs = PyTuple_New(1)) == NULL ) {
	  Py_DECREF(atomsList) ;
	  return NULL ;
	}
	/* call atomDefFF with this name to return that FF */
	PyTuple_SetItem(newargs, 0, Py_BuildValue("O",value)) ;
	if( (ffval = atomDefFF(NULL, newargs)) == NULL ) {
	  Py_DECREF(newargs) ; Py_DECREF(atomsList) ;
	  return NULL ;
	}
	Py_DECREF(newargs) ;
	if( ! PyTuple_Check(ffval) || PyTuple_Size(ffval) != 3 ) {
	  j2frac = Py_BuildValue("d", 0.) ;
	  Cval = PyList_New(0) ;
	} else {
	  j2frac = PyTuple_GetItem(ffval, 1) ;
	  if( ! PyFloat_Check(j2frac) && ! PyInt_Check(j2frac) )
	    j2frac = Py_BuildValue("d", 0.) ;
	  else
	    j2frac = Py_BuildValue("d", PyFloat_AsDouble(j2frac)) ;
	  Cval = PyTuple_GetItem(ffval, 2) ;
	  if( ! PyList_Check(Cval) )
	    Cval = PyList_New(0) ;
	  else
	    Cval = PyList_GetSlice(Cval, 0, PyList_Size(Cval)) ;
	}
	Py_DECREF(ffval) ;
      } else if( PyFloat_Check(value) || PyInt_Check(value) ) {
	/* j2frac [coefs] */
	j2frac = PyTuple_GetItem(args, 1) ;
	if( ! PyFloat_Check(j2frac) && ! PyInt_Check(j2frac) )
	  j2frac = Py_BuildValue("d", 0.) ;
	else
	  j2frac = Py_BuildValue("d", PyFloat_AsDouble(j2frac)) ;
	if( narg > 2 ) {
	  Cval = PyTuple_GetItem(args, 2) ;
	  if( ! PyList_Check(Cval) )
	    Cval = PyList_New(0) ;
	  else
	    Cval = PyList_GetSlice(Cval, 0, PyList_Size(Cval)) ;
	} else {
	  Cval = NULL ;	  
	}
      } else {
	/* error */
	PyErr_SetString(PyExc_ValueError, "editFF bad args") ;
	Py_DECREF(atomsList) ;
	return NULL ;
      }
      if( Cval != NULL ) {
	if( (newargs = PyTuple_New(3)) == NULL ) {
	  Py_DECREF(atomsList) ; return NULL ;
	}
      } else {
	if( (newargs = PyTuple_New(2)) == NULL ) {
	  Py_DECREF(atomsList) ; return NULL ;
	}
      }
      gname = Py_BuildValue("s", atomgroupi->atom->name) ;
      PyTuple_SetItem(newargs, 0, gname) ;
      PyTuple_SetItem(newargs, 1, j2frac) ;
      if( Cval != NULL ) PyTuple_SetItem(newargs, 2, Cval) ;

      if( (result = atomDefFF(NULL, newargs)) == NULL ) {
	Py_DECREF(atomsList) ;
	Py_DECREF(newargs) ;
	return NULL ;
      }
      Py_DECREF(result) ; Py_DECREF(newargs) ;
    } else if( etyp == 7 ) {
      /* edit DW */
      if ( (value = PyTuple_GetItem(args, 1)) == NULL ) {
	Py_DECREF(atomsList) ;	
	return NULL ;
      }
      if ( PyString_Check(value) ) {
	/* srcname */
	if( (newargs = PyTuple_New(1)) == NULL ) {
	  Py_DECREF(atomsList) ;	
	  return NULL ;
	}
	PyTuple_SetItem(newargs, 0, Py_BuildValue("O",value)) ;
	if( (ffval = atomDefDW(NULL, newargs)) == NULL ) {
	  Py_DECREF(newargs) ; Py_DECREF(atomsList) ;
	  return NULL ;
	}
	Py_DECREF(newargs) ;
	if( ! PyTuple_Check(ffval) || PyTuple_Size(ffval) != 2 ) {
	  Cval = PyList_New(0) ;
	} else {
	  Cval = PyTuple_GetItem(ffval, 1) ;
	  if( ! PyList_Check(Cval) )
	    Cval = PyList_New(0) ;
	  else
	    Cval = PyList_GetSlice(Cval, 0, PyList_Size(Cval)) ;
	}
	Py_DECREF(ffval) ;
      } else if( PyList_Check(value) ) {
	/* [coefs] */
	Cval = PyTuple_GetItem(args, 1) ;
	if( ! PyList_Check(Cval) )
	  Cval = PyList_New(0) ;
	else
	  Cval = PyList_GetSlice(Cval, 0, PyList_Size(Cval)) ;
      } else if( PyFloat_Check(value) || PyInt_Check(value) ) {
	Cval = PyList_New(narg-1) ;
	for( j=1 ; j<narg ; j++ ) {
	  if ( (value = PyTuple_GetItem(args, j)) == NULL ) {
	    Py_DECREF(atomsList) ; Py_DECREF(Cval) ;
	    return NULL ;
	  }
	  if( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
	    PyErr_SetString(PyExc_ValueError,
			    "editDW args after 1st must be numeric or list") ;
	    Py_DECREF(atomsList) ; Py_DECREF(Cval) ;
	    return NULL ;
	  }
	  PyList_SetItem(Cval, j-1,
			 Py_BuildValue("d", PyFloat_AsDouble(value))) ;
	}
      } else {
	/* error */
	PyErr_SetString(PyExc_ValueError, "editDW bad args") ;
	Py_DECREF(atomsList) ;
	return NULL ;
      }
      if( (newargs = PyTuple_New(2)) == NULL ) {
	Py_DECREF(atomsList) ;
	return NULL ;
      }
      gname = Py_BuildValue("s", atomgroupi->atom->name) ;
      PyTuple_SetItem(newargs, 0, gname) ;
      PyTuple_SetItem(newargs, 1, Cval) ;
      if( (result = atomDefDW(NULL, newargs)) == NULL ) {
	Py_DECREF(newargs) ; Py_DECREF(atomsList) ;
	return NULL ;
      }
      Py_DECREF(result) ; Py_DECREF(newargs) ;
    }

    if( etyp != 5 && etyp != 3 && etyp < 6 ) {
      if( ! sym ) isub = -1 ;
      if( putAtom(name, atomgroupi, atomgroupi->atom, isub, r, s, m, c) < 0 )
	return NULL;
    }

    if( etyp < 6 ) {
      if ( (groupList = MakeAtomList(atomgroupi->atom->name, atomgroupi->isub,
				     atoms)) == NULL ) {
	Py_DECREF(atomsList) ;
	return NULL ;
      }
    } else if( etyp == 6 ) {
      /* FF */
      if( (newargs = PyTuple_New(1)) == NULL ) {
	Py_DECREF(atomsList) ;	
	return NULL ;
      }
      PyTuple_SetItem(newargs, 0, Py_BuildValue("s",atomgroupi->atom->name)) ;
      if ( (groupList = atomDefFF(NULL, newargs)) == NULL ) {
	Py_DECREF(newargs) ; Py_DECREF(atomsList) ;
	return NULL ;
      }
      Py_DECREF(newargs) ;
    } else if( etyp == 7 ) {
      /* DW */
      if( (newargs = PyTuple_New(1)) == NULL ) {
	Py_DECREF(atomsList) ;	
	return NULL ;
      }
      PyTuple_SetItem(newargs, 0, Py_BuildValue("s",atomgroupi->atom->name)) ;
      if ( (groupList = atomDefDW(NULL, newargs)) == NULL ) {
	Py_DECREF(newargs) ; Py_DECREF(atomsList) ;
	return NULL ;
      }
      Py_DECREF(newargs) ;
    }

    PyList_SET_ITEM(atomsList, nret, groupList) ;
    nret++ ;
  }
  return atomsList ;
}



static void resetAtomsList(ATOMSlist *atomslst)
{
  ATOMgroup *atomgroup ;
  myATOM *atom ;
  double r[3], s[3] ;
  int i, j ;
  for( i=0 ; i<atomslst->n ; i++ ) {
    atomgroup = atomslst->atomslist + i ;
    /* use the first atom which is the original coordinate generator */
    atom = atomgroup->atoms ;
    for( j=0 ; j<3 ; j++ ) { r[j] = atom->r[j] ; s[j] = atom->s[j] ; }
    putAtom(atomgroup->atom->name, atomgroup, atomgroup->atom, 0,
	    r, s, atom->m, atom->c) ;
  }
}

static int getSpaceGroupIndex(int ispc, int iset)
{
  int i, j ;
  for( i=0 ; i<NspaceGrpsLoaded ; i++ ) {
    if( spcgrps[i].spcnum == ispc ) {
      for( j=0 ; j<6 ; j++ ) {
	if( spcgrps[i].setting[j] == 0 ) break ;
	if( spcgrps[i].setting[j] == iset ) return i ;
      }
    }
  }
  /* failed to find that setting so return first instance of that spcnum */
  for( i=0 ; i<NspaceGrpsLoaded ; i++ ) {
    if( spcgrps[i].spcnum == ispc ) return i ;
  }
  return -1 ;
}



static void addOut( Qs *qs )
{
  int i, j, N ;
  Q *qpt ;
  N = NQ + qs->n ;
  if( N > NQalloc ) {
    SQh = (double*)realloc(SQh, N*sizeof(double)) ;
    SQk = (double*)realloc(SQk, N*sizeof(double)) ;
    SQl = (double*)realloc(SQl, N*sizeof(double)) ;
    SQq = (double*)realloc(SQq, N*sizeof(double)) ;
    SQp = (double*)realloc(SQp, N*sizeof(double)) ;
    SQcs = (double*)realloc(SQcs, N*sizeof(double)) ;
    SQcsp = (double*)realloc(SQcsp, N*sizeof(double)) ;
    SQcsm = (double*)realloc(SQcsm, N*sizeof(double)) ;
    SQpp = (double*)realloc(SQpp, N*sizeof(double)) ;
    SQmm = (double*)realloc(SQmm, N*sizeof(double)) ;
    SQpm = (double*)realloc(SQpm, N*sizeof(double)) ;
    SQmp = (double*)realloc(SQmp, N*sizeof(double)) ;
    SQtt = (double*)realloc(SQtt, N*sizeof(double)) ;
    SQom = (double*)realloc(SQom, N*sizeof(double)) ;
    SQlor = (double*)realloc(SQlor, N*sizeof(double)) ;
    NQalloc = N ;
  }
  for( i=NQ, j=0 ; i<N ; i++, j++ ) {
    qpt = qs->q + j ;
    SQh[i] = qpt->hkl[0] ;
    SQk[i] = qpt->hkl[1] ;
    SQl[i] = qpt->hkl[2] ;
    SQq[i] = qpt->q ;
    SQp[i] = qpt->outofplane ;
    SQcs[i] = qpt->csp + qpt->csm ;
    SQcsp[i] = qpt->csp ;
    SQcsm[i] = qpt->csm ;
    SQpp[i] = qpt->pp ;
    SQmm[i] = qpt->mm ;
    SQpm[i] = qpt->pm ;
    SQmp[i] = qpt->mp ;
    SQtt[i] = qpt->tt ;
    SQom[i] = qpt->om ;
    SQlor[i] = qpt->lor ;
  }
  NQ = N ;
}

static void initQsums(Q *qpt)
{
  static DCMPLX czero = {0., 0.} ;
  qpt->S0 = qpt->Sx = qpt->Sy = qpt->Sz = czero ;
  qpt->Sp0 = qpt->Spx = qpt->Spy = qpt->Spz = czero ;
  qpt->csp = qpt->csm = qpt->pp = qpt->mm = qpt->pm = qpt->mp = 0. ;
}
static void initQ(Q *qpt)
{
  int i ;
  double arg1, arg2, alpha, omega, delta, twotheta, theta ;
  double sint, sin2t ;
  static double omvec[3] = { 0., 0., 1. } ;
  //static double RADTODEG = 57.295779513 ;

  qpt->q = recvec( qpt->hkl, qpt->qun ) ;
  if( *powderavg ) {
    /*
      recompute orientation rhsetup and pbcalc so that q is in scatt plane
     */
    rhsetup( qpt->hkl, ip3.hkl, ip2.hkl ) ;
    pbcalc() ;
  }

  if( fabs(qpt->q) < smallestQ ) {
    /* define qpt-qun perp to ip2 and ip3 so still can have non-zero mag vec */
    vecpro(ip3.qun, ip2.qun, qpt->qun) ;
  }

  qpt->outofplane = dotpro( qpt->qun, ip3.qun ) ;
  if( fabs(qpt->outofplane) < smallestQ ) qpt->outofplane = 0. ;

  arg1 = qpt->q/ki/2. ;

  if( arg1 <= 1. ) {
    delta = acos(arg1) ;
    theta = asin(arg1) ;
    twotheta = 2.*theta ;
    qpt->tt = RADTODEG*twotheta ;
    if( fabs(qpt->q) < smallestQ ) {
      sint = sin2t = 1. ;
    } else {
      sint = sin(theta) ;
      sin2t = sin(twotheta) ;
    }
    if( *TTlorentz ) {
      if( *Qlorentz ) qpt->lor = 1./sint ; /* equal Q-steps */
      else qpt->lor = 1./sin2t ; /* angle steps */
      if( *Plorentz ) qpt->lor *= 1./sin2t ;
    } else {
      /* omega scans irrelevant for powder */
      if( *Qlorentz ) qpt->lor = 1./cos(theta) ; /* equal Q-steps */
      else qpt->lor = 1./sin2t ; /* angle steps */
    }
  } else {
    delta = 0. ;
    qpt->lor = 0. ;
    qpt->tt = 0. ;
  }

  arg1 = dotpro( qpt->qun, ip2.qun ) ;
  arg2 = dotpro( qpt->qun, ip1.qun ) ;
  
  if( fabs(arg1) < 1.e-6 && fabs(arg2) < 1.e-6 )
    alpha = 0. ;
  else
    alpha = atan2(arg1,arg2) ;
  
  /* alpha is signed angle of q from hkl1 towards hkl2 */
  /* delta is the angle between ki and q */
  
  omega = fmod(PIOVER2 - delta + alpha, TWOPI) ;
  qpt->om = RADTODEG*omega ;

  if( pbq == 0 ) {
    omvec[0] = cos(omega) ;
    omvec[1] = sin(omega) ;
    
    /* calc pol beam components in xtal Cartesian coordinates */
    matvec( trx, omvec, qpt->pbx ) ;
    matvec( try, omvec, qpt->pby ) ;
    matvec( trz, omvec, qpt->pbz ) ;

  } else {
    for( i=0 ; i<3 ; i++ )
      {
	qpt->pbx[i] = pbx[i] ;
	if( pbq > 0 ) { qpt->pbz[i] = qpt->qun[i] ; }
	else { qpt->pbz[i] = -qpt->qun[i] ; }
      }
    vecpro( qpt->pbz, qpt->pbx, qpt->pby ) ;
  }
}

static int outColumnsList( int *ilst )
{
  /*
   * go thru the calccolumns and make the ordered list of columns to write
   * to do this successively find the lowest remaining index
   */
  int i, ncols, inext, more ;

  ncols = 0 ;
  inext = 1 ;
  more = 1 ;
  while( more ) {
    more = 0 ;
    for( i=0 ; i<ncalccolumns ; i++ ) {
      if( calccolumns[i].i == inext ) {
	ilst[ncols] = i ;
	ncols++ ;
	more = 1 ;
	break ;
      } else if( calccolumns[i].i > inext ) {
	more = 1 ;
      }
    }
    inext++ ;
  }
  return ncols ;
}

static int writeSetup(FILE *fpt)
{
  PyObject *self, *args, *result, *result2, *result3, *result4 ;
  FILE *f ;
  
  int cols[16] ;
  int i, j, ncols ;

  f = fpt ;
  if( f == NULL ) f = stdout ;

  /*
    header:
    write lattice orient wavelength polarize
    write the flag values for magonly nuconly powderavg
    write spacegroup or NONE depending on ALLPOSITIONS flag
    write atoms in atomlist one per commented line
   */

  self = NULL ;
  if ( (args = PyTuple_New(0)) == NULL ) return 0 ;

  if ( (result = lattice(self, args)) == NULL ) {
    Py_DECREF(args) ; return 0 ;
  }
  fprintf(f, "lattice(") ; PyObject_Print(result, f, 0) ; fprintf(f, ")\n") ;
  Py_DECREF(result) ;

  if ( (result = orient(self, args)) == NULL ) {
    Py_DECREF(args) ; return 0 ;
  }
  fprintf(f, "orient(") ; PyObject_Print(result, f, 0) ; fprintf(f, ")\n") ;
  Py_DECREF(result) ;

  if ( (result = polarize(self, args)) == NULL ) {
    Py_DECREF(args) ; return 0 ;
  }
  fprintf(f, "polarize(") ; PyObject_Print(result, f, 0) ; fprintf(f, ")\n");
  Py_DECREF(result) ;

  if ( (result = wavelength(self, args)) == NULL ) {
    Py_DECREF(args) ; return 0 ;
  }
  fprintf(f, "wavelength(") ; PyObject_Print(result,f, 0) ; fprintf(f,")\n");
  Py_DECREF(result) ;

  if ( (result = spacegroup(self, args)) == NULL ) {
    Py_DECREF(args) ; return 0 ;
  }
  result2 = PyTuple_GetItem(result, 0) ;
  /* PyList_GetItem I dont own the ref */
  result3 = PyList_GetItem(result2, 0) ;
  result4 = PyTuple_GetItem(result3, 0) ;
  fprintf(f, "spacegroup(") ; PyObject_Print(result4,f, 0) ; fprintf(f,")\n");
  Py_DECREF(result) ;
  Py_DECREF(args) ;
  return 1 ;
}

static int writeAtoms(FILE *fpt)
{
  /*
    rewrite this to segment line output, one line per atomgroup in the atomlist
    if spacegroup is inactive get return from atomList(igroup) and use this as
    arg to put command. This is what to do for magnetic structures where the
    symmetry of spacegroup is destroyed and put() return will no longer
    regenerate the magnetic structure. Recall put only returns first atom of
    each group assumes rest are generated by spacegroup.
  */
  int i, j, nl ;
  PyObject *self, *args, *result, *alist, *atomLIST ;
  FILE *f ;

  f = fpt ;
  if( f == NULL ) f = stdout ;
  self = NULL ;
  if ( *usespacegroup ) {
    if ( (args = PyTuple_New(0)) == NULL ) return 0 ;
    if ( (result = atomPut(self, args)) == NULL ) {
      Py_DECREF(args) ; return 0 ;
    }
    Py_DECREF(args) ;
    nl = PyList_Size(result) ;
    for( i=0 ; i<nl ; i++ ) {
      if( (atomLIST = PyList_GetItem(result, i)) == NULL ) {
	Py_DECREF(result) ; return 0 ;
      }
      fprintf(f, "put(") ; PyObject_Print(atomLIST, f, 0) ; fprintf(f,")\n");
      /* I dont own the atomLIST ref returned from PyList_GetItem */
    }
    Py_DECREF(result) ;
  } else {
    for( i=0 ; i<atomslist.n ; i++ ) {
      if ( (args = PyTuple_New(1)) == NULL ) return 0 ;
      PyTuple_SET_ITEM(args, 0, Py_BuildValue("i", i+1)) ;
      if ( (result = atomList(self, args)) == NULL ) {
	Py_DECREF(args) ; return 0 ;
      }
      Py_DECREF(args) ;
      if ( PyList_Size(result) < 1 ) {
	PyErr_SetString(PyExc_ValueError, "atomList returned empty list") ;
	Py_DECREF(result) ; return 0 ;
      }
      alist = PyList_GetItem(result, 0) ;
      nl = PyList_Size(alist) ;
      for( j=0 ; j<nl ; j++ ) {
	if( (atomLIST = PyList_GetItem(alist, j)) == NULL ) return 0 ;
	fprintf(f, "put(") ; PyObject_Print(atomLIST, f, 0) ; fprintf(f,")\n");
      }
      Py_DECREF(result) ;
    }
  }
  return 1 ;
}

static int writeQ(FILE *fpt)
{
  PyObject *self, *args, *result ;
  FILE *f ;
  f = fpt ;
  if( f == NULL ) f = stdout ;
  self = NULL ;
  if ( (args = PyTuple_New(0)) == NULL ) return 0 ;
  if ( (result = Qadd(self, args)) == NULL ) return 0 ;
  fprintf(f, "Q(") ; PyObject_Print(result,f, 0) ; fprintf(f,")\n");
  return 1 ;
}

static int writeCalc(FILE *fpt)
{
  PyObject *self, *args, *result ;
  FILE *f ;
  double Qmax ;
  
  int cols[16] ;
  int i, j, ncols ;

  f = fpt ;
  if( f == NULL ) f = stdout ;

  /*
    header:
    write lattice orient wavelength polarize
    write the flag values for magonly nuconly powderavg
    write spacegroup or NONE depending on ALLPOSITIONS flag
    write atoms in atomlist one per commented line
   */


  if( *printSetup ) if( ! writeSetup(f) ) return 0 ;
  if( *printAtoms ) if( ! writeAtoms(f) ) return 0 ;

  self = NULL ;
  if ( (args = PyTuple_New(0)) == NULL ) return 0 ;
  if ( (result = calcFlags(self, args)) == NULL ) {
    Py_DECREF(args) ; return 0 ;
  }
  Py_DECREF(args) ;
  fprintf(f, "calcflags(") ; PyObject_Print(result,f, 0) ; fprintf(f,")\n");
  Py_DECREF(result) ;

  if ( (args = PyTuple_New(0)) == NULL ) return 0 ;
  if ( (result = calcColumns(self, args)) == NULL ) {
    Py_DECREF(args) ; return 0 ;
  }
  Py_DECREF(args) ;
  fprintf(f, "calccolumns(") ; PyObject_Print(result,f, 0) ; fprintf(f,")\n");
  Py_DECREF(result) ;

  if ( (args = PyTuple_New(0)) == NULL ) return 0 ;
  if ( (result = lorentzFlags(self, args)) == NULL ) {
    Py_DECREF(args) ; return 0 ;
  }
  Py_DECREF(args) ;
  fprintf(f, "lorentz(") ; PyObject_Print(result,f, 0) ; fprintf(f,")\n");
  Py_DECREF(result) ;


  if ( (ncols = outColumnsList(cols)) < 1 ) return 1 ;

  /* construct the column labels line */
  for( j=0 ; j<ncols ; j++ ) fprintf(f, "%s", outcols[cols[j]].l) ;
  fprintf(f, "\n") ;

  Qmax = 10. ;
  if( lambda > 0. ) Qmax = 2.*TWOPI/lambda ;
  for( i=0 ; i<NQ ; i++ ) {
    /* dont print any resuls for Q > Qmax = 4Pi/lambda */
    if( (*(outcols[3].d))[i] > Qmax ) continue ;
    for( j=0 ; j<ncols ; j++ ) fprintf(f, "%12g", (*(outcols[cols[j]].d))[i]) ;
    fprintf(f, "\n") ;
  }
  return 1 ;
}

static int writeFlags(FILE *fpt)
{
  PyObject *self, *args, *result ;
  FILE *f ;
  
  int cols[16] ;
  int i, j, ncols ;

  f = fpt ;
  if( f == NULL ) f = stdout ;

  /*
    header:
    write lattice orient wavelength polarize
    write the flag values for magonly nuconly powderavg
    write spacegroup or NONE depending on ALLPOSITIONS flag
    write atoms in atomlist one per commented line
   */

  self = NULL ;
  if ( (args = PyTuple_New(0)) == NULL ) return 0 ;

  if ( (result = calcFlags(self, args)) == NULL ) {
    Py_DECREF(args) ; return 0 ;
  }

  fprintf(f, "calcflags(") ; PyObject_Print(result, f, 0) ; fprintf(f,")\n");
  Py_DECREF(result) ;

  if ( (result = lorentzFlags(self, args)) == NULL ) {
    Py_DECREF(args) ; return 0 ;
  }
  fprintf(f, "lorentz(") ; PyObject_Print(result,f, 0) ; fprintf(f,")\n");
  Py_DECREF(result) ;

  if ( (result = calcColumns(self, args)) == NULL ) {
    Py_DECREF(args) ; return 0 ;
  }
  fprintf(f, "calccolumns(") ; PyObject_Print(result,f, 0) ; fprintf(f,")\n");
  Py_DECREF(result) ;

  if ( (result = unknownColumns(self, args)) == NULL ) {
    Py_DECREF(args) ; return 0 ;
  }
  fprintf(f, "unknowncolumns("); PyObject_Print(result,f, 0); fprintf(f,")\n");
  Py_DECREF(result) ;

  if ( (result = atomPutFlags(self, args)) == NULL ) {
    Py_DECREF(args) ; return 0 ;
  }
  fprintf(f, "putflags(") ; PyObject_Print(result,f, 0) ; fprintf(f,")\n");
  Py_DECREF(result) ;

  Py_DECREF(args) ;
  return 1 ;
}

static FILE *openFileFromArgs(PyObject *args, char *fstr)
{
  int narg ;
  PyObject *arg ;
  char *str, *dot ;
  char errbuf[64] ;
  char filebuf[512] ;
  FILE *fpt ;

  narg = PyTuple_Size(args) ;
 
  if ( narg < 1 ) {
    strcpy(errbuf, "filename arg required for ") ;
    strcat(errbuf, fstr) ;
    PyErr_SetString(PyExc_ValueError, errbuf) ;
    return NULL ;
  }
  if ( (arg = PyTuple_GetItem(args, 0)) == NULL ) return NULL ;  
  if ( ! PyString_Check(arg) ) {
    strcpy(errbuf, "filename string arg required for ") ;
    strcat(errbuf, fstr) ;
    PyErr_SetString(PyExc_ValueError, errbuf);
    return NULL ;
  }
  if( (str = PyString_AsString(arg)) == NULL || strlen(str) < 1 ) {
    strcpy(errbuf, "failed to convert arg to string in ") ;
    strcat(errbuf, fstr) ;
    PyErr_SetString(PyExc_ValueError, errbuf);
    return NULL ;
  }

  if( (dot = strrchr(str, '.')) == NULL || strcmp(dot+1,"py") ) {
    /* fix up so extension is .py */
    if( strlen(str) > 507 ) {
      PyErr_SetString(PyExc_ValueError, "filename is too long") ;
      return NULL ;
    }
    strcpy(filebuf, str) ;
    strcat(filebuf, ".py") ;
    str = filebuf ;
  }
  if( (fpt = fopen(str, "w")) == NULL ) {
    strcpy(errbuf, "failed to open file in ") ;
    strcat(errbuf, fstr) ;
    PyErr_SetString(PyExc_ValueError, errbuf);
    return NULL ;
  }
  return fpt ;
}

static double dwCalc(Q *qpt, DebyeWallerFactor *dw)
{
  int i ;
  double shkl[3], tv[3] ;
  if( qpt == NULL || dw == NULL ) return 1. ;
  if( dw->iso ) return exp(-0.5*dw->u[0][0]*qpt->q*qpt->q) ;
  for( i=0 ; i<3 ; i++ ) shkl[i] = qpt->hkl[i]*rlat[i] ;
  matvec( dw->u, shkl, tv ) ;
  return exp(-0.5*(dotpro(tv, shkl))) ;
}

static double ffCalc( Q *qpt, MagFormFactor *ff )
{
  int i ;
  double ssq ;
  double ffval, ffval2 ;
  double l2s ;
  static double fourpi = 12.5663706144 ;

  /*
   * normalized mag form factor as per J. Brown
   * dipole approx f(Q) = j0(Q) + (L/(L+2S)) j2(Q) 
   * OR            f(Q) = j0(Q) + (2/g - 1) j2(Q)
   * where g is the Lande g-factor = 2 for spin only and 1 for orbital only
   * GSAS manual has (1 - g/2)??? for this factor and leaves out the s^2 ???
   * A. Zheludev's web site (2/g - 1) for rare-earths with moment gJ
   * for transition metals however use
   *  (1 - 2/g) where g is effective gyro-mag ratio so moment is gS
   * j0 = Aexp(-as^2) + Bexp(-bs^2) + Cexp(-cs^2) + D
   * j2 = same * s^2
   * s = sinthet/lambda  q =4pi sinthet/lambda
   * OR s = q/4pi
   */

  if( qpt == NULL || ff == NULL ) return 1. ;
  ffval = 0 ;
  ssq = qpt->q / fourpi ;
  ssq *= ssq ;

  for( i=0 ; i<7 ; i+=2 )
    {
      if( i >= ff->n ) return ffval ;
      if( i == 6 ) ffval += ff->coefs[i] ;
      else ffval += ff->coefs[i]*exp(-(ff->coefs[i+1]*ssq)) ;
    }

  ffval2 = 0. ;
  for( i=7 ; i<14 ; i+=2 )
    {
      if( i >= ff->n ) return (ffval + ff->j2frac*ssq*ffval2) ;
      if( i == 13 ) ffval2 += ff->coefs[i] ;
      else ffval2 += ff->coefs[i]*exp(-(ff->coefs[i+1]*ssq)) ;
    }
  return (ffval + ff->j2frac*ssq*ffval2) ;
}

static void docalc()
{
  int i ;
  if ( *resetcalc ) {
    if( *selectQ > 0 && selectQlist.n > 0 ) {
      for( i=0 ; i<selectQlist.n ; i++ )
	{
	  if( selectQlist.i[i] < 1 || selectQlist.i[i] > qlist.n )
	    continue ;
	  calcQ(&atomslist, qlist.qlist + selectQlist.i[i] - 1) ;
	}
    } else {
      for( i=0 ; i<qlist.n ; i++ ) calcQ(&atomslist, qlist.qlist + i) ;
    }
  }

  NQ = 0 ;
  if( *selectQ > 0 && selectQlist.n > 0 ) {
    for( i=0 ; i<selectQlist.n ; i++ )
      {
	if( selectQlist.i[i] < 1 || selectQlist.i[i] > qlist.n )
	  continue ;
	addOut( qlist.qlist + selectQlist.i[i] - 1 ) ;
      }
  } else {
    for( i=0 ; i<qlist.n ; i++ ) addOut( qlist.qlist + i ) ;
  }

  if( *sortQ ) sortQmag() ;
}

static void cscalc( Q *qpt )
{
  DCMPLX yi ;
  static DCMPLX icmplx = { 0., 1. } ;
  yi = Cmul( icmplx, qpt->Sy ) ;
  qpt->pm = Cabsq( Cadd(qpt->Sx, yi) ) ;
  qpt->mp = Cabsq( Csub(qpt->Sx, yi) ) ;

  if( *powderavg ) {
    yi = Cmul( icmplx, qpt->Spy ) ;
    qpt->pm = (qpt->pm + Cabsq( Cadd(qpt->Spx, yi) ))/2. ;
    qpt->mp = (qpt->mp + Cabsq( Csub(qpt->Spx, yi) ))/2. ;
    qpt->pp = Cabsq(qpt->S0) + (Cabsq(qpt->Sz) + Cabsq(qpt->Spz))/2. ;
    qpt->mm = qpt->pp ;
  } else {
    qpt->pp = Cabsq( Cadd(qpt->S0, qpt->Sz) ) ;
    qpt->mm = Cabsq( Csub(qpt->S0, qpt->Sz) ) ;
  }
  if( *lorentzcalc ) {
    qpt->pm *= qpt->lor ;
    qpt->mp *= qpt->lor ;
    qpt->pp *= qpt->lor ;
    qpt->mm *= qpt->lor ;
  }
  /*
    divide all polarize beam CS by 2 so that adding them will give
    standard unpolarized beam result
  */
  qpt->pm /= 2. ;
  qpt->mp /= 2. ;
  qpt->pp /= 2. ;
  qpt->mm /= 2. ;

  if(qpt->pm < smallestSQ) qpt->pm = 0. ;
  if(qpt->mp < smallestSQ) qpt->mp = 0. ;
  if(qpt->pp < smallestSQ) qpt->pp = 0. ;
  if(qpt->mm < smallestSQ) qpt->mm = 0. ;
  qpt->csp = qpt->pp + qpt->pm ;
  qpt->csm = qpt->mp + qpt->mm ;
  if(qpt->csp < smallestSQ) qpt->csp = 0. ;
  if(qpt->csm < smallestSQ) qpt->csm = 0. ;
  return ;
}

static void calcAtomOneQ(int add, ATOMgroup *atomgroup, Qs *qs)
{
  /*
    add or subtract from phase sums for all Q in current Qlist
    For powder avg where average the cross-section for all rotations
    of xtal about the q-vector, can show that spin-flip cross sections
    are just avg of normal and case where magvec is rotated 90deg to q x magvec
    Also for non-flip CSs the cross term between S0 and Sz averages to zero and
    one averages only the magnetic part |Sz|^2 using the perp magvec.
    This is because the CS is averaged over the rotation of the magnetic vector
    where CSflip has terms  SUMatomsi,j faci facj pu.mi(phi) pu.mj(phi)
    where pu is polarization component x,y,z
    mi(phi) = cos(phi) m0 + sin(phi) (qhat x m0) 
    describes the rotation of magvec about q.
    Thus the cos^2 and sin^2 terms survive the average and each give 1/2
   */
  int i, k, ff, dw ;
  double dwfac, fffac, rfac ;
  double mvx, mvy, mvz ;
  double pmvx, pmvy, pmvz ;
  double rvec[3], magvec[3], pmagvec[3] ;
  DCMPLX fac ;

  Q *qpt ;
  myATOM *atom ;
  ATOMdef *atomdef ;

  atomdef = atomgroup->atom ;
  //b = atomdef->b ;
  /*
    b imaginary values are tabulated at 1.8 A (2200 m/s)
    see M.&L. appendix A for k dep of b
  */

  ff = 0 ;
  if( atomdef->ff != NULL && atomdef->ff->n > 0 ) ff = 1 ;
  dw = 0 ;
  if( atomdef->dw != NULL ) dw = 1 ;

  for( i=0 ; i<atomgroup->natoms ; i++ ) {
    atom = atomgroup->atoms + i ;

    for( k=0 ; k<qs->n ; k++ ) {
      qpt = qs->q + k ;
      dwfac = 1. ;
      if( dw ) dwfac = dwCalc(qpt, atomdef->dw) ;
      fac = RCmul(dwfac ,RCexp((double)TWOPI*dotpro(qpt->hkl, atom->r))) ;
      
      if( ! (*magonly) ) {
	if( add ) qpt->S0 =
		    Cadd( qpt->S0, Cmul(fac, RCmul(atom->c, atom->b)) ) ;
	else      qpt->S0 = 
		    Csub( qpt->S0, Cmul(fac, RCmul(atom->c, atom->b)) ) ;
      }
      
      if( ! (*nuconly) )
	{
	  fffac = 1. ;
	  if( ff ) fffac = ffCalc(qpt, atomdef->ff) ;
	  rfac = Bohrmagtob * atom->c * atom->m * fffac ;
	  if( rfac > 0. )
	    {
	      fac = RCmul( rfac, fac ) ;
	      vecpro( atom->ss, qpt->qun, rvec ) ;
	      vecpro( qpt->qun, rvec, magvec ) ;
	      
	      mvx = dotpro( qpt->pbx, magvec ) ;
	      mvy = dotpro( qpt->pby, magvec ) ;
	      mvz = dotpro( qpt->pbz, magvec ) ;

	      if( *powderavg ) {
		vecpro( qpt->qun, magvec, pmagvec ) ;
		pmvx = dotpro( qpt->pbx, pmagvec ) ;
		pmvy = dotpro( qpt->pby, pmagvec ) ;
		pmvz = dotpro( qpt->pbz, pmagvec ) ;
		if( add ) {
		  qpt->Spx = Cadd( qpt->Spx, RCmul( pmvx, fac ) ) ;
		  qpt->Spy = Cadd( qpt->Spy, RCmul( pmvy, fac ) ) ;
		  qpt->Spz = Cadd( qpt->Spz, RCmul( pmvz, fac ) ) ;
		} else {
		  qpt->Spx = Csub( qpt->Spx, RCmul( pmvx, fac ) ) ;
		  qpt->Spy = Csub( qpt->Spy, RCmul( pmvy, fac ) ) ;
		  qpt->Spz = Csub( qpt->Spz, RCmul( pmvz, fac ) ) ;
		}
	      }

	      if( add ) {
		qpt->Sx = Cadd( qpt->Sx, RCmul( mvx, fac ) ) ;
		qpt->Sy = Cadd( qpt->Sy, RCmul( mvy, fac ) ) ;
		qpt->Sz = Cadd( qpt->Sz, RCmul( mvz, fac ) ) ;
	      } else {
		qpt->Sx = Csub( qpt->Sx, RCmul( mvx, fac ) ) ;
		qpt->Sy = Csub( qpt->Sy, RCmul( mvy, fac ) ) ;
		qpt->Sz = Csub( qpt->Sz, RCmul( mvz, fac ) ) ;
	      }
	    }
	}
      cscalc( qpt ) ;
    }
  }
}

static void calcQ( ATOMSlist *atomslst, Qs *qs )
{
  /*
    calc new phase sums for this Qseq using all atoms
   */
  int i ;
  if( atomslst->n <= 0 ) return ;
  for( i=0 ; i<qs->n ; i++ ) { initQsums(qs->q + i) ; initQ(qs->q + i) ; }
  for( i=0 ; i<atomslst->n ; i++ )
    calcAtomOneQ(1, atomslst->atomslist+i, qs) ;
}


static int readQFile(char *fil)
{
  /* need to fix this */
  char buf[256], obuf[256] ;
  static char sqcmnd[] = { "sq q " } ;
  FILE *fpt ;
  if( (fpt = fopen(fil, "r")) == NULL ) return 0 ;

  while( fgets(buf, 512, fpt) ) {
    if( buf[0] == '#' ) continue ;
    if( strlen(buf) < 1 ) continue ;
    /*
      parse any sq commands or atom commands
      like sq lattice OR sq spacegroup
    */
    if( strncmp(buf, "sq", 2) == 0 ) {
      //if( Tcl_GlobalEval(interp, buf) != TCL_OK ) return 0 ;
      continue ;
    }
    if( strncmp(buf, "atom", 4) == 0 ) {
      //if( Tcl_GlobalEval(interp, buf) != TCL_OK ) return 0 ;
      continue ;
    }
    /* else try to parse h k l or seq spec */
    strcpy(obuf, sqcmnd) ;
    if( buf[0] == 'Q' ) strcat(obuf, buf+1) ;
    else strcat(obuf, buf) ;
    //if( Tcl_GlobalEval(interp, obuf) ) return 0 ;
  }
  fclose(fpt) ;
  return 1 ;
}


static int writeQFile(char *fil, Qlist *qlist)
{
  int i ;
  FILE *fpt ;

  Qs *qs ;

  if( (fpt = fopen(fil, "w")) == NULL ) return 0 ;

  for( i=0 ; i<qlist->n ; i++ ) {
    qs = qlist->qlist + i ;
    fprintf(fpt, "Qadd(%g, %g, %g, %g, %g, %g, %d)\n",
	    qs->hkl[0],qs->hkl[1],qs->hkl[2],
	    qs->step[0],qs->step[1],qs->step[2],qs->n) ;
  }
  fclose(fpt) ;
  return 1 ;
}

static int nextQ(Qlist *qlst)
{
  static char location[] = "nextQ" ;
  int n ;
  static ATOMgroup *zerogroupPtr = NULL ;

  if( zerogroupPtr == NULL ) {
    zerogroupPtr = (ATOMgroup *)calloc(1, sizeof(ATOMgroup)) ;
    if ( memExc(isNULL(zerogroupPtr), location) ) return -1 ;
  }
  qlst->n++ ;
  n = qlst->n ;
  if( n > qlst->nalloc ) {
    qlst->qlist = (Qs *)realloc(qlst->qlist, n*sizeof(Qs)) ;
    if ( memExc(isNULL(qlst->qlist), location) ) return -1 ;
    qlst->nalloc = n ;
    qlst->qlist[n-1].n = 0 ;
    qlst->qlist[n-1].q = NULL ;
  }
  return (n - 1) ;
}

static void swapDbls( double *a, double *b )
{
  double temp ;
  temp = *a ;
  *a = *b ;
  *b = temp ;
}
static void sortQmag()
{
  /* sort the SQ arrays on SQq. There are NQ elements */
  int i, j ;
  for( i=0 ; i<NQ-1 ; i++ ) {
    for( j=i+1 ; j<NQ ; j++ ) {
      if( SQq[i] > SQq[j] ) {
	swapDbls(SQh+i, SQh+j) ;
	swapDbls(SQk+i, SQk+j) ;
	swapDbls(SQl+i, SQl+j) ;
	swapDbls(SQq+i, SQq+j) ;
	swapDbls(SQp+i, SQp+j) ;
	swapDbls(SQcsp+i, SQcsp+j) ;
	swapDbls(SQcsm+i, SQcsm+j) ;
	swapDbls(SQpp+i, SQpp+j) ;
	swapDbls(SQpm+i, SQpm+j) ;
	swapDbls(SQmm+i, SQmm+j) ;
	swapDbls(SQmp+i, SQmp+j) ;
	swapDbls(SQcs+i, SQcs+j) ;
	swapDbls(SQtt+i, SQtt+j) ;
	swapDbls(SQom+i, SQom+j) ;
	swapDbls(SQlor+i, SQlor+j) ;
      }
    }
  }
}

static void sortQlist(Qlist *qlst)
{
  /*
    qlist is a list of qlists
    sort on the first qmag in each qlist
   */

  int i, j, n ;
  Qs temp ;

  Qs *qlisti, *qlistj ;
  Q  *qli, *qlj ;

  Qs *qlist0 = qlst->qlist ;
  n = qlst->n ;
  for( i=0 ; i<n-1 ; i++ ) {
    qlisti = qlist0 + i ; 
    qli = qlisti->q ;
    for( j=i+1 ; j<n ; j++ ) {
      qlistj = qlist0 + j ;
      qlj = qlistj->q ;
      if( qlj->q <= qli->q ) {
	temp = *qlisti ;
	*qlisti = *qlistj ;
	*qlistj = temp ;
	qli = qlisti->q ;
      }
    }
  }
}


static int symbolInSpcGrp(int spcGrpN, char *symbol)
{
  int i ;
  Subgrp *subgrp ;
  //printf("symbolInSpcGrp symbol=%c\n", symbol[0]) ;
  for( i=0 ; i<spcgrps[spcGrpN].nsub ; i++ ) {
    subgrp = &(spcgrps[spcGrpN].sub[i]) ;
    //printf("check against Wyck=%c for subgrp %d\n",subgrp->Wyck,i) ;
    if( subgrp->Wyck == symbol[0] ) return i ;
  }
  return -1 ;
}


static ATOMgroup *nextAtomGroup(ATOMSlist *atomslst)
{
  static char location[] = "nextAtomGroup" ;
  int n ;

  static ATOMgroup *zerogroupPtr = NULL ;
  if( zerogroupPtr == NULL )
    zerogroupPtr = (ATOMgroup *)calloc(1, sizeof(ATOMgroup)) ;

  atomslst->n++ ;
  n = atomslst->n ;
  if( n > atomslst->nalloc ) {
    atomslst->atomslist =
      (ATOMgroup *)realloc(atomslst->atomslist, n*sizeof(ATOMgroup)) ;
    if ( memExc(isNULL(atomslst->atomslist), location) ) return NULL ;
    atomslst->nalloc = n ;
    atomslst->atomslist[n-1] = *zerogroupPtr ;
  }
  return atomslst->atomslist + n - 1 ;
}

static int matchString(char *chk, char **types)
{
  int i, nc, nchk, it ;
  char *type ;
  char chklc[64], typelc[64] ;

  if( chk == NULL || types == NULL ) return -1 ;
  nc = strlen(chk) ;
  for( i=0 ; i<nc ; i++ ) chklc[i] = tolower(chk[i]) ;
  chklc[nc] = '\0' ;
  type = *types ;
  it = 0 ;
  while( type != NULL ) {
    nchk = strlen(type) ;
    for( i=0 ; i<nchk ; i++ ) typelc[i] = tolower(type[i]) ;
    typelc[nchk] = '\0' ;
    if( ! strncmp(chklc, typelc, nc) ) return it + 1 ;
    it++ ;
    type = types[it] ;
  }
  return -1 ;
}
static int readCheckType(char *fil)
{
  static char *test[] = {
    "SGR", "SPCGRP", "S GRUP", "RGNR", "!ATOM", "_ATOM", "CRS", "SYMM", "SPGP",
    "ATOMPUT"
  } ;
  static char *alt1[] = {
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, "SPGR", "PUT"
  } ;
  static char *alt2[] = {
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, "SGRP", "SPACEGROUP"
  } ;
  static int ntypes = 10 ;

  FILE *fpt ;
  char buf[512] ;
  char *str ;
  int i, nc ;

  if( (fpt = fopen(fil, "r")) == NULL ) return -1 ;
  /*
    "icsd",       SGR
    "lazy",       spcgrp
    "Xtal-3D",    S GRUP
    "powdercell", RGNR
    "fullprof",   !Atom
    "cif",        _atom
    "gsas",       CRS
    "shel-X",     SYMM
    "DRAWxtl",    spgp
    "pysq"        atom  OR put
  */

  while( fgets(buf, 511, fpt) ) {
    if( buf[0] == '#' ) continue ;
    str = buf ;
    while( *str != '\0' ) { *str = toupper(*str) ; str++ ; }
    str = buf ;
    while( *str == ' ' || *str == '\t' ) str++ ;
    for( i=0 ; i<ntypes ; i++ ) {
      nc = strlen(test[i]) ;
      if( ! strncmp(test[i], str, nc) ) { fclose(fpt) ; return i + 1 ; }
      if( alt1[i] != NULL ) {
	nc = strlen(alt1[i]) ;
	if( ! strncmp(alt1[i], str, nc) ) { fclose(fpt) ; return i + 1 ; }
      }
      if( alt2[i] != NULL ) {
	nc = strlen(alt2[i]) ;
	if( ! strncmp(alt2[i], str, nc) ) { fclose(fpt) ; return i + 1 ; }
      }
    }
  }
  fclose(fpt) ;
  return 0 ;
}

static int readUNKNOWN(char *fil, int ph) ;
static int readICSD(char *fil, int ph) ;
static int readLAZY(char *fil, int ph) ;
static int readXTAL3D(char *fil, int ph) ;
static int readPOWDERCELL(char *fil, int ph) ;
static int readFULLPROF(char *fil, int ph) ;
static int readCIF(char *fil, int ph) ;
static int readGSAS(char *fil, int ph) ;
static int readSHELX(char *fil, int ph) ;
static int readDRAWXTL(char *fil, int ph) ;


static int badorder = 0 ;

static int readStructureFile(char *fil, char *typ, int ph)
{
  int ityp, ir ;
  static char *filetypes[] = {
    "icsd",
    "lazy",
    "Xtal-3D",
    "powdercell",
    "fullprof",
    "cif",
    "gsas",
    "shel-X",
    "DRAWxtl",
    "pysq",
    NULL
  } ;
  static int (*readfun[])(char *file, int ph) = {
    readUNKNOWN,
    readICSD,
    readLAZY,
    readXTAL3D,
    readPOWDERCELL,
    readFULLPROF,
    readCIF,
    readGSAS,
    readSHELX,
    readDRAWXTL,
    readPYSQ
  } ;


  ityp = -1 ;
  if( typ != NULL ) ityp = matchString(typ, filetypes) ;
  if( ityp < 0 ) ityp = readCheckType(fil) ;
  if( ityp < 0 ) {
    PyErr_SetString(PyExc_ValueError, "open failure from readCheckType") ;
    return -1 ;  /* open error */
  }
  if( ityp == 0 ) {
    /*
      for unknown file type first try readPYSQ
      if import fails, try readUNKNOWN i.e. procede with ityp == 0
    */
    if( readPYSQ(fil, ph) >= 0 ) return 1 ;
    if( PyErr_Occurred() ) PyErr_Clear() ;
  }

  if( (ir = readfun[ityp](fil, ph)) < 0 ) {
    badorder = 0 ;
    return ir ;
  }
  /*
    have to reread the file if cmnd order is bad
    N.B. atoms are reset when a new spacegroup is read
    but we go ahead and reread the whole file anyway
    This is probably not necessary because spacegroup
    resets any existing atom list.
  */
  if( badorder ) {
    if( (ir = readfun[ityp](fil, ph)) < 0 ) {
      badorder = 0 ;
      return ir ;
    }
    badorder = 0 ;
  }
  return 1 ;
}

#define READERR { fclose(fpt) ; return -1 ; }
#define READRET { fclose(fpt) ; if( (cel && !acel) || (grp && !agrp) ) badorder = 1 ; return 1 ; }


static int lineparse(char *line, WordList *wordlist, int parenflag) ;
static PyObject *tupleFromDbls(int n, double *d) ;
static void clearAtomList() ;
static int atomNumLookup(char *buf) ;


static int readUNKNOWN(char *fil, int ph)
{
  /*
    only hope to read atoms from unknown format file using
    the unknowncolumns designations
    i.e. we arent looking for lattice or spacegroup.
    Hopefully the user can look at the file and figure out
    which columns correspond to which data item.
    The atomname has to parse as a valid element.

    look for lines with
    atomname as word number unknowncolumns[0].i
    X        as word number unknowncolumns[1].i
    Y        as word number unknowncolumns[2].i
    Z        as word number unknowncolumns[3].i
    Biso     as word number unknowncolumns[4].i
    occ      as word number unknowncolumns[5].i
    mom      as word number unknowncolumns[6].i
    sx       as word number unknowncolumns[7].i
    sy       as word number unknowncolumns[8].i
    sz       as word number unknowncolumns[9].i

    recall that sx,sy,sz for this program should be along the
    direct space axes which may not be orthogonal
  */

  PyObject *tpl, *stpl, *mval, *oval, *name, *name2 ;
  PyObject *args, *dwargs ;
  FILE *fpt ;
  int i, nb, nd, nw, ready ;
  int cel, grp, acel, agrp ;
  int *isn ;
  char buf[512] ;
  char *str, *end ;
  double dbls[6], s[3], Uiso ;
  Flag *uk ;

  cel = grp = 0 ;
  acel = agrp = 1 ;
  s[0] = s[1] = s[2] = 0. ;

  if( (fpt = fopen(fil, "r")) == NULL ) {
    PyErr_SetString(PyExc_ValueError, "readUNKNOWN failed to open file") ;
    return -1 ;
  }

  clearAtomList() ;
  uk = unknowncolumns ;
  while(1) {
    if( ! fgets(buf, 511, fpt) ) break ;
    if ( lineparse(str, &wordlist, 0) < 0 ) READERR ;
    nw = wordlist.nwords ;
    isn = wordlist.isnum ;
    /* has to have atomname X Y Z at minimum */
    if ( uk[0].i < 1 || nw < uk[0].i ) continue ;
    if ( uk[1].i < 1 || nw < uk[1].i ) continue ;
    if ( uk[2].i < 1 || nw < uk[2].i ) continue ;
    if ( uk[3].i < 1 || nw < uk[3].i ) continue ;

    /* check that name looks like an element */
    if( atomNumLookup(wordlist.words[uk[0].i].word) <= 0 ) continue ;

    /* check that the X Y Z columns are numeric */
    if ( ! isn[uk[1].i-1] ) continue ;
    if ( ! isn[uk[2].i-1] ) continue ;
    if ( ! isn[uk[3].i-1] ) continue ;

    if( (args = PyTuple_New(5)) == NULL ) return -1 ;
    if ( (tpl = PyTuple_New(3)) == NULL ) { Py_DECREF(args) ; READERR ; }

    for( i=0 ; i<3 ; i++ )
      PyTuple_SetItem(tpl, i, Py_BuildValue("d", wordlist.nums[uk[i+1].i-1])) ;

    name = Py_BuildValue("s", wordlist.words[uk[0].i-1].word) ;
    PyTuple_SET_ITEM(args, 0, name) ;
    PyTuple_SET_ITEM(args, 1, tpl) ;

    if ( uk[7].i > 0 && uk[7].i <= nw && isn[uk[7].i-1] &&
	 uk[8].i > 0 && uk[8].i <= nw && isn[uk[8].i-1] &&
	 uk[9].i > 0 && uk[9].i <= nw && isn[uk[9].i-1] ) {
      if ( (stpl = PyTuple_New(3)) == NULL ) { Py_DECREF(args) ; READERR ; }
      for( i=0 ; i<3 ; i++ )
	PyTuple_SetItem(stpl, i,
			Py_BuildValue("d", wordlist.nums[uk[i+7].i-1]));
    } else {
      if( (stpl = tupleFromDbls(3, s)) == NULL ) { Py_DECREF(args); READERR ; }
    }
    if ( uk[6].i > 0 && uk[6].i <= nw && isn[uk[6].i-1] )
      mval = Py_BuildValue("d", wordlist.nums[uk[6].i-1]) ;
    else
      mval = Py_BuildValue("d", 0.) ;
    PyTuple_SET_ITEM(args, 2, stpl) ;
    PyTuple_SET_ITEM(args, 3, mval) ;

    if ( uk[5].i > 0 && uk[5].i <= nw && isn[uk[5].i-1] )
      oval = Py_BuildValue("d", wordlist.nums[uk[5].i-1]) ;
    else
      oval = Py_BuildValue("d", 1.) ;
    PyTuple_SET_ITEM(args, 4, oval) ;

    if(atomPut(NULL, args) == NULL) { Py_DECREF(args) ; READERR ; }
    Py_DECREF(args) ;

    if( (dwargs = PyTuple_New(2)) == NULL ) return -1 ;
    name2 = Py_BuildValue("s", wordlist.words[uk[0].i-1].word) ;
    PyTuple_SET_ITEM(dwargs, 0, name2) ;

    if ( uk[4].i > 0 && uk[4].i <= nw && isn[uk[4].i-1] )
      /* convert Biso to Uiso */
      Uiso = wordlist.nums[uk[4].i-1]/8./3.14159/3.14159 ;
    else
      Uiso = 0. ;
    PyTuple_SET_ITEM(dwargs, 1, Py_BuildValue("d", Uiso)) ;

    /*
      NB this will set DW for all atoms of a given element
      ie not allowing for different site symmetries
      to do that would require Copying atomdefs to create diff sites
    */
    if(atomDefDW(NULL, dwargs) == NULL) { Py_DECREF(dwargs) ; READERR ; }
    Py_DECREF(dwargs) ;

    if( ! cel ) acel = 0 ;
    if( ! grp ) agrp = 0 ;
  }

  fclose(fpt) ;
  if( (cel && !acel) || (grp && !agrp) ) badorder = 1 ;
  return 1 ;
}

/*
    readUNKNOWN,
    readICSD,
    readLAZY,
    readXTAL3D,
    readPOWDERCELL,
    readFULLPROF,
    readCIF,
    readGSAS,
    readSHELX,
    readDRAWXTL,
    readPYSQ
*/



static int readICSD(char *fil, int ph)
{
  /* I assume keywords for ICSD are in caps */
  /* need to add parsing of atomname num B= Biso lines */
  PyObject *tpl, *stpl, *mval, *oval, *name, *name2 ;
  PyObject *args, *dwargs ;
  PyObject *result ;
  FILE *fpt ;
  int nb, nd, ready ;
  int cel, grp, acel, agrp ;
  char buf[512] ;
  char *str, *end ;
  double dbls[6], s[3], Uiso ;

  cel = grp = 0 ;
  acel = agrp = 1 ;
  s[0] = s[1] = s[2] = 0. ;

  if( (fpt = fopen(fil, "r")) == NULL ) {
    PyErr_SetString(PyExc_ValueError, "readICSD failed to open file") ;
    return -1 ;
  }

  ready = 0 ;
  while(1) {
    if( ! ready ) if( ! fgets(buf, 511, fpt) ) break ;
    ready = 0 ;
    str = buf ;
    while( *str == ' ' ) str++ ;
    if( ! strncmp(str, "TITL", 4) ) {
      nd = strlen(str+4) ;
      if( strlen(str) > 511 ) str[511] = '\0' ;
      strcpy(titlebuf, str+4) ;
      while( fgets(buf, 511, fpt) ) {
	str = buf ;
	nb = 0 ;
	while( *str == ' ' ) { str++ ; nb++ ; }
	if ( nb < 2 ) { ready = 1 ; break ; }
	str-- ;
	nd += strlen(str) ;
	if( nd > 511 ) break ;
	strcat(titlebuf, str) ;
      }
      tpl = Py_BuildValue("(s)", titlebuf) ;
      if( (result = title(NULL, tpl)) == NULL ) { Py_DECREF(tpl) ; READERR ; }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
    } else if( ! strncmp(str, "CELL", 4) ) {
      str += 3 ;
      nd = 0 ;
      while( nd < 6 && (str = strchr(str+1,'=')) != NULL ) {
	if( sscanf(str+1, "%lf", dbls+nd) < 1 ) break ;
	nd++ ;
      }
      //if ( lineparse(buf, &wordlist, 0) < 0 ) return -1 ;
      tpl = tupleFromDbls(nd, dbls) ;
      if( (result = lattice(NULL, tpl)) == NULL ) {
	Py_DECREF(tpl) ; READERR ;
      }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
      cel = 1 ;
    } else if( ! strncmp(str, "SGR", 3) ) {
      str += 3 ;
      if( (end = strchr(str, '(')) != NULL ) *end = '\0' ;
      tpl = Py_BuildValue("(s)", str) ;
      if( (result = spacegroup(NULL, tpl)) == NULL) {
	Py_DECREF(tpl) ; READERR ;
      }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
      grp = 1 ;
    } else if( ! strncmp(str, "PARM", 4) ) {
      /* clear the atomlist */
      clearAtomList() ;


      while( fgets(buf, 511, fpt) ) {
	str = buf ;
	nb = 0 ;
	while( *str == ' ' ) { str++ ; nb++ ; }
	if ( nb < 2 ) { ready = 1 ; break ; }
	if ( lineparse(str, &wordlist, 0) < 0 ) READERR ;
	/* word 1 is atomName, word 5-7 coords, word 8 occ */
	if ( wordlist.nwords < 7 ) continue ;
	name = Py_BuildValue("s", wordlist.words[0].word) ;
	if( atomNumLookup(wordlist.words[0].word) <= 0 ) {
	  Py_DECREF(name) ; continue ;
	}
	if( (args = PyTuple_New(5)) == NULL ) { Py_DECREF(name) ; READERR ; }
	if( (tpl = tupleFromDbls(3, wordlist.nums+4)) == NULL ) {
	  Py_DECREF(name) ; Py_DECREF(args) ;
	  READERR ;
	}
	if( (stpl = tupleFromDbls(3, s)) == NULL ) {
	  Py_DECREF(name) ; Py_DECREF(tpl) ; Py_DECREF(args) ; READERR ;
	}
	mval = Py_BuildValue("d", 0.) ;
	if ( wordlist.nwords > 7 ) oval = Py_BuildValue("d", wordlist.nums[7]);
	else oval = Py_BuildValue("d", 1.) ;
	PyTuple_SetItem(args, 0, name) ;
	PyTuple_SetItem(args, 1, tpl) ;
	PyTuple_SetItem(args, 2, stpl) ;
	PyTuple_SetItem(args, 3, mval) ;
	PyTuple_SetItem(args, 4, oval) ;
	if( (result = atomPut(NULL, args)) == NULL) {
	  Py_DECREF(args) ; READERR ;
	}
	Py_DECREF(args) ;
	Py_DECREF(result) ;
	if( ! cel ) acel = 0 ;
	if( ! grp ) agrp = 0 ;
      }

    } else if( ! strncmp(str, "RVAL", 4) ) {

      while( fgets(buf, 511, fpt) ) {
	str = buf ;
	nb = 0 ;
	while( *str == ' ' ) { str++ ; nb++ ; }
	if ( nb < 2 ) { ready = 1 ; break ; }
	if ( lineparse(str, &wordlist, 0) < 0 ) { Py_DECREF(dwargs); READERR ;}
	/* word 1 is atomName, word 4 is Biso */
	if ( wordlist.nwords < 4 ) continue ;
	name2 = Py_BuildValue("s", wordlist.words[0].word) ;
	if( atomNumLookup(wordlist.words[0].word) <= 0 ) {
	  Py_DECREF(name2) ; continue ;
	}
	if( (dwargs = PyTuple_New(2)) == NULL ) {
	  Py_DECREF(name2) ; READERR ;
	}
	PyTuple_SetItem(dwargs, 0, name2) ;
	if ( wordlist.isnum[3] )
	  /* convert Biso to Uiso */
	  Uiso = wordlist.nums[3]/8./3.14159/3.14159 ;
	else
	  Uiso = 0. ;
	/*
	  although Py_Build creates a new ref it is immediately stolen
	  by SetItem
	*/
	PyTuple_SetItem(dwargs, 1, Py_BuildValue("d", Uiso)) ;
	/*
	  NB this will set DW for all atoms of a given element
	  ie not allowing for different site symmetries
	  to do that would require Copying atomdefs to create diff sites
	*/
	if( (result = atomDefDW(NULL, dwargs)) == NULL) {
	  Py_DECREF(dwargs) ; READERR ;
	}
	Py_DECREF(dwargs) ;
	Py_DECREF(result) ;
      }
    }
  }
  fclose(fpt) ;
  if( (cel && !acel) || (grp && !agrp) ) badorder = 1 ;
  return 1 ;
}

static int readLAZY(char *fil, int ph)
{
  PyObject *tpl, *stpl, *mval, *oval ;
  PyObject *args, *dwargs, *name, *name2, *result ;
  FILE *fpt ;
  int nb, nd, ready ;
  int cel, grp, acel, agrp ;
  char buf[512] ;
  char *str, *end ;
  double dbls[6], s[3], Uiso ;

  cel = grp = 0 ;
  acel = agrp = 1 ;
  s[0] = s[1] = s[2] = 0. ;


  if( (fpt = fopen(fil, "r")) == NULL ) {
    PyErr_SetString(PyExc_ValueError, "readLAZY failed to open file") ;
    return -1 ;
  }
  clearAtomList() ;
  ready = 0 ;
  while(1) {
    if( ! ready ) if( ! fgets(buf, 511, fpt) ) break ;
    ready = 0 ;
    str = buf ;
    while( *str == ' ' ) str++ ;

    if( ! strncmp(str, "TITLE", 5) ) {
      nd = strlen(str+5) ;
      if( strlen(str) > 511 ) str[511] = '\0' ;
      strcpy(titlebuf, str+5) ;
      tpl = Py_BuildValue("(s)", titlebuf) ;
      if( (result = title(NULL, tpl)) == NULL ) { Py_DECREF(tpl) ; READERR ; }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
    } else if( ! strncmp(str, "CELL", 4) ) {
      if ( lineparse(str+4, &wordlist, 0) < 0 ) { fclose(fpt) ; return -1 ; }
      if ( wordlist.nwords > 0 ) {
	tpl = tupleFromDbls(wordlist.nwords, wordlist.nums) ;
	if( (result = lattice(NULL, tpl)) == NULL ) {
	  Py_DECREF(tpl) ; READERR ;
	}
	Py_DECREF(tpl) ;
	Py_DECREF(result) ;
      }
      cel = 1 ;
    } else if( ! strncmp(str, "SPCGRP", 6) ) {
      tpl = Py_BuildValue("(s)", str+6) ;
      if( (result = spacegroup(NULL, tpl)) == NULL) {
	Py_DECREF(tpl) ; READERR ;
      }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
      grp = 1 ;
    } else if( ! strncmp(str, "ATOM", 4) ) {
      if ( lineparse(str+4, &wordlist, 0) < 0 ) READERR ;
      /* word 1 is atomName, word 3-5 coords, 6 Biso 7 occ */
      if ( wordlist.nwords < 5 ) continue ;
      if( (args = PyTuple_New(5)) == NULL ) READERR ;
      name = Py_BuildValue("s", wordlist.words[0].word) ;
      if( atomNumLookup(wordlist.words[0].word) <= 0 ) {
	Py_DECREF(name) ; continue ;
      }
      if( (tpl = tupleFromDbls(3, wordlist.nums+2)) == NULL ) {
	Py_DECREF(name) ; READERR ;
      }
      if ( wordlist.nwords > 6 ) oval = Py_BuildValue("d", wordlist.nums[6]);
      else oval = Py_BuildValue("d", 1.) ;
      if( (stpl = tupleFromDbls(3, s)) == NULL ) {
	Py_DECREF(name) ; Py_DECREF(tpl) ; Py_DECREF(oval) ; READERR ;
      }
      mval = Py_BuildValue("d", 0.) ;

      PyTuple_SetItem(args, 0, name) ;
      PyTuple_SetItem(args, 1, tpl) ;
      PyTuple_SetItem(args, 2, stpl) ;
      PyTuple_SetItem(args, 3, mval) ;
      PyTuple_SetItem(args, 4, oval) ;
      if( (result = atomPut(NULL, args)) == NULL) {
	Py_DECREF(args) ; READERR ;
      }
      Py_DECREF(args) ;
      Py_DECREF(result) ;

      if( (dwargs = PyTuple_New(2)) == NULL ) return -1 ;
      name2 = Py_BuildValue("s", wordlist.words[0].word) ;
      PyTuple_SetItem(dwargs, 0, name2) ;

      if ( wordlist.nwords > 5 && wordlist.isnum[5] ) {
	/* convert Biso to Uiso */
	Uiso = wordlist.nums[5]/8./3.14159/3.14159 ;
	PyTuple_SetItem(dwargs, 1, Py_BuildValue("d", Uiso)) ;
      } else {
	PyTuple_SetItem(dwargs, 1, Py_BuildValue("d", 0.)) ;
      }
      /*
	NB this will set DW for all atoms of a given element
	ie not allowing for different site symmetries
	to do that would require Copying atomdefs to create diff sites
      */
      if( (result = atomDefDW(NULL, dwargs)) == NULL) {
	Py_DECREF(dwargs) ; READERR ;
      }
      Py_DECREF(dwargs) ;
      Py_DECREF(result) ;
      if( ! cel ) acel = 0 ;
      if( ! grp ) agrp = 0 ;
    }
  }
  fclose(fpt) ;
  if( (cel && !acel) || (grp && !agrp) ) badorder = 1 ;
  return 1 ;
}
static int readXTAL3D(char *fil, int ph)
{
  PyObject *tpl, *stpl, *mval, *oval ;
  PyObject *args, *dwargs, *name, *name2, *result ;
  FILE *fpt ;
  int nb, nd, ready ;
  int cel, grp, acel, agrp ;
  char buf[512] ;
  char *str, *end ;
  double dbls[6], s[3], Uiso ;

  cel = grp = 0 ;
  acel = agrp = 1 ;
  s[0] = s[1] = s[2] = 0. ;


  if( (fpt = fopen(fil, "r")) == NULL ) {
    PyErr_SetString(PyExc_ValueError, "readXTAL3D failed to open file") ;
    return -1 ;
  }
  clearAtomList() ;
  ready = 0 ;
  while(1) {
    if( ! ready ) if( ! fgets(buf, 511, fpt) ) break ;
    ready = 0 ;
    str = buf ;
    while( *str == ' ' ) str++ ;

    if( ! strncmp(str, "N", 1) ) {
      nd = strlen(str+1) ;
      if( strlen(str) > 511 ) str[511] = '\0' ;
      strcpy(titlebuf, str+5) ;
      tpl = Py_BuildValue("(s)", titlebuf) ;
      if( (result = title(NULL, tpl)) == NULL ) { Py_DECREF(tpl) ; READERR ; }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
    } else if( ! strncmp(str, "C", 1) ) {
      if ( lineparse(str+1, &wordlist, 0) < 0 ) READERR ;
      if ( wordlist.nwords > 0 ) {
	tpl = tupleFromDbls(wordlist.nwords, wordlist.nums) ;
	if( (result = lattice(NULL, tpl)) == NULL ) {
	  Py_DECREF(tpl) ; READERR ;
	}
	Py_DECREF(tpl) ;
	Py_DECREF(result) ;
      }
      cel = 1 ;
    } else if( ! strncmp(str, "S GRUP", 6) ) {
      tpl = Py_BuildValue("(s)", str+6) ;
      if( (result = spacegroup(NULL, tpl)) == NULL) {
	Py_DECREF(tpl) ; READERR ;
      }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
      grp = 1 ;
    } else if( ! strncmp(str, "A", 1) ) {
      if ( lineparse(str+1, &wordlist, 0) < 0 ) READERR ;
      /* word 1 is atomName, word 3-5 coords, 6 Biso 7 occ */
      if ( wordlist.nwords < 4 ) continue ;

      if( atomNumLookup(wordlist.words[0].word) <= 0 ) continue ;
      if( ! wordlist.isnum[1] ) continue ;
      if( ! wordlist.isnum[2] ) continue ;
      if( ! wordlist.isnum[3] ) continue ;

      if( (stpl = tupleFromDbls(3, s)) == NULL ) READERR ;
      if( (tpl = tupleFromDbls(3, wordlist.nums+1)) == NULL ) {
	Py_DECREF(stpl) ; READERR ;
      }
      name = Py_BuildValue("s", wordlist.words[0].word) ;
      mval = Py_BuildValue("d", 0.) ;
      if ( wordlist.nwords > 5 ) oval = Py_BuildValue("d", wordlist.nums[5]);
      else oval = Py_BuildValue("d", 1.) ;

      if( (args = PyTuple_New(5)) == NULL ) {
	Py_DECREF(stpl) ; Py_DECREF(tpl) ; Py_DECREF(name) ; Py_DECREF(mval) ;
	Py_DECREF(oval) ; READERR ;
      }
      PyTuple_SetItem(args, 0, name) ;
      PyTuple_SetItem(args, 1, tpl) ;
      PyTuple_SetItem(args, 2, stpl) ;
      PyTuple_SetItem(args, 3, mval) ;
      PyTuple_SetItem(args, 4, oval) ;
      if( (result = atomPut(NULL, args)) == NULL) {
	Py_DECREF(args) ; READERR ;
      }
      Py_DECREF(args) ;
      Py_DECREF(result) ;

      if( (dwargs = PyTuple_New(2)) == NULL ) READERR ;
      name2 = Py_BuildValue("s", wordlist.words[0].word) ;
      PyTuple_SetItem(dwargs, 0, name2) ;

      if ( wordlist.nwords > 4 && wordlist.isnum[4] ) {
	/* convert Biso to Uiso */
	Uiso = wordlist.nums[4]/8./3.14159/3.14159 ;
	PyTuple_SetItem(dwargs, 1, Py_BuildValue("d", Uiso)) ;
      } else {
	PyTuple_SetItem(dwargs, 1, Py_BuildValue("d", 0.)) ;
      }
      /*
	NB this will set DW for all atoms of a given element
	ie not allowing for different site symmetries
	to do that would require Copying atomdefs to create diff sites
      */
      if( (result = atomDefDW(NULL, dwargs)) == NULL) {
	Py_DECREF(dwargs) ; READERR ;
      }
      Py_DECREF(dwargs) ;
      Py_DECREF(result) ;
      if( ! cel ) acel = 0 ;
      if( ! grp ) agrp = 0 ;
    }
  }
  fclose(fpt) ;
  if( (cel && !acel) || (grp && !agrp) ) badorder = 1 ;
  return 1 ;
}
static int readPOWDERCELL(char *fil, int ph)
{
  PyObject *tpl, *stpl, *mval, *oval ;
  PyObject *args, *dwargs, *name, *name2, *result ;
  FILE *fpt ;
  int nb, nd, ready, ispc, iset ;
  int cel, grp, acel, agrp ;
  char buf[512] ;
  char *str, *end ;
  double dbls[6], s[3], Uiso ;

  cel = grp = 0 ;
  acel = agrp = 1 ;
  s[0] = s[1] = s[2] = 0. ;


  if( (fpt = fopen(fil, "r")) == NULL ) {
    PyErr_SetString(PyExc_ValueError, "readPOWDERCELL failed to open file") ;
    return -1 ;
  }
  clearAtomList() ;
  ready = 0 ;
  while(1) {
    if( ! ready ) if( ! fgets(buf, 511, fpt) ) break ;
    ready = 0 ;
    str = buf ;
    while( *str == ' ' ) str++ ;

    if( ! strncmp(str, "CELL", 4) ) {
      if ( lineparse(str+4, &wordlist, 0) < 0 ) READERR ;
      if ( wordlist.nwords > 0 ) {
	tpl = tupleFromDbls(wordlist.nwords, wordlist.nums) ;
	if( (result = lattice(NULL, tpl)) == NULL ) {
	  Py_DECREF(tpl) ; READERR ;
	}
	Py_DECREF(tpl) ;
	Py_DECREF(result) ;
      }
      cel = 1 ;
    } else if( ! strncmp(str, "RGNR", 4) ) {
      if ( (nd = sscanf(str+4, "%d %d", &ispc, &iset)) < 1 ) continue ;
      if ( nd == 2 ) tpl = Py_BuildValue("(ii)", ispc, iset) ;
      else tpl = Py_BuildValue("(i)", ispc) ;
      if( (result = spacegroup(NULL, tpl)) == NULL) {
	Py_DECREF(tpl) ; READERR ;
      }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
      grp = 1 ;
    } else {
      /* should be atom */
      if ( lineparse(str, &wordlist, 0) < 0 ) READERR ;
      /* word 1 is atomName, 2 atomType  word 3-5 coords, 6 occ 7 Biso */
      if ( wordlist.nwords < 5 ) continue ;
      if( atomNumLookup(wordlist.words[0].word) <= 0 ) continue ;
      if( ! wordlist.isnum[2] ) continue ;
      if( ! wordlist.isnum[3] ) continue ;
      if( ! wordlist.isnum[4] ) continue ;

      if( (tpl = tupleFromDbls(3, wordlist.nums+2)) == NULL ) READERR ;
      if( (stpl = tupleFromDbls(3, s)) == NULL ) {
	Py_DECREF(tpl) ; READERR ;
      }
      name = Py_BuildValue("s", wordlist.words[0].word) ;
      mval = Py_BuildValue("d", 0.) ;
      if ( wordlist.nwords > 5 ) oval = Py_BuildValue("d", wordlist.nums[5]);
      else oval = Py_BuildValue("d", 1.) ;

      if( (args = PyTuple_New(5)) == NULL ) {
	Py_DECREF(stpl) ; Py_DECREF(tpl) ; Py_DECREF(name) ; Py_DECREF(mval) ;
	Py_DECREF(oval) ; READERR ;
      }
      PyTuple_SetItem(args, 0, name) ;
      PyTuple_SetItem(args, 1, tpl) ;
      PyTuple_SetItem(args, 2, stpl) ;
      PyTuple_SetItem(args, 3, mval) ;
      PyTuple_SetItem(args, 4, oval) ;

      if( (result = atomPut(NULL, args)) == NULL) {
	Py_DECREF(args) ; READERR ;
      }
      Py_DECREF(args) ;
      Py_DECREF(result) ;

      if( (dwargs = PyTuple_New(2)) == NULL ) READERR ;
      name2 = Py_BuildValue("s", wordlist.words[0].word) ;
      PyTuple_SetItem(dwargs, 0, name2) ;
      if ( wordlist.nwords > 6 && wordlist.isnum[6] ) {
	/* convert Biso to Uiso */
	Uiso = wordlist.nums[6]/8./3.14159/3.14159 ;
	PyTuple_SetItem(dwargs, 1, Py_BuildValue("d", Uiso)) ;
      } else {
	PyTuple_SetItem(dwargs, 1, Py_BuildValue("d", 0.)) ;
      }
      /*
	NB this will set DW for all atoms of a given element
	ie not allowing for different site symmetries
	to do that would require Copying atomdefs to create diff sites
      */
      if( (result = atomDefDW(NULL, dwargs)) == NULL) {
	Py_DECREF(dwargs) ; READERR ;
      }
      Py_DECREF(dwargs) ;
      Py_DECREF(result) ;
      if( ! cel ) acel = 0 ;
      if( ! grp ) agrp = 0 ;
    }
  }
  fclose(fpt) ;
  if( (cel && !acel) || (grp && !agrp) ) badorder = 1 ;
  return 1 ;
}
static int readFULLPROF(char *fil, int ph)
{
  /* currently this will only read first phase in file */
  PyObject *tpl, *stpl, *mval, *oval ;
  PyObject *args, *dwargs, *name, *name2, *result ;
  PyObject *alst, *alist ;
  FILE *fpt ;
  int nb, nd, ready ;
  int cel, grp, acel, agrp, imult ;
  char buf[512] ;
  char *str, *end ;
  double dbls[6], s[3], Uiso, occ, mult ; ;

  cel = grp = 0 ;
  acel = agrp = 1 ;
  s[0] = s[1] = s[2] = 0. ;


  if( (fpt = fopen(fil, "r")) == NULL ) {
    PyErr_SetString(PyExc_ValueError, "readFULLPROF failed to open file") ;
    return -1 ;
  }
  clearAtomList() ;

  if( ! fgets(titlebuf, 511, fpt) ) READRET ;
  tpl = Py_BuildValue("(s)", titlebuf) ;
  if( (result = title(NULL, tpl)) == NULL ) { Py_DECREF(tpl) ; READERR ; }
  Py_DECREF(tpl) ;
  Py_DECREF(result) ;

  /* find the spacegroup line */
  while(1) {
    if( ! fgets(buf, 511, fpt) ) READRET ;
    if( (str = strstr(buf, "!Space group symbol")) != NULL ) {
      *str = '\0' ;
      tpl = Py_BuildValue("(s)", buf) ;
      if( (result = spacegroup(NULL, tpl)) == NULL) {
	Py_DECREF(tpl) ; READERR ;
      }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
      grp = 1 ;
      break ;
    }
  }

  /* next line is header labels for atoms ? is this always true */
  if( ! fgets(buf, 511, fpt) ) READRET ;
  /* now read atoms until comment line */
  nd = 1 ;
  while(1) {
    if( ! fgets(buf, 511, fpt) ) READRET ;
    if( buf[0] == '!' ) break ; /* end of atomslist */
    if ( lineparse(buf, &wordlist, 0) < 0 ) READERR ;
    /* word 1 is atomName, 2 atomtype word 3-5 coords, 6 Biso 7 mult*occ */
    if ( wordlist.nwords < 5 ) continue ;
    if( atomNumLookup(wordlist.words[0].word) <= 0 ) continue ;
    if( ! wordlist.isnum[2] ) continue ;
    if( ! wordlist.isnum[3] ) continue ;
    if( ! wordlist.isnum[4] ) continue ;
    if( (tpl = tupleFromDbls(3, wordlist.nums+2)) == NULL ) READERR ;
    if( (stpl = tupleFromDbls(3, s)) == NULL ) {
      Py_DECREF(tpl) ; READERR ;
    }
    name = Py_BuildValue("s", wordlist.words[0].word) ;
    mval = Py_BuildValue("d", 0.) ;
    oval = Py_BuildValue("d", 1.) ;
    if( (args = PyTuple_New(5)) == NULL ) {
      Py_DECREF(stpl) ; Py_DECREF(tpl) ; Py_DECREF(name) ; Py_DECREF(mval) ;
      Py_DECREF(oval) ; READERR ;
    }
    PyTuple_SetItem(args, 0, name) ;
    PyTuple_SetItem(args, 1, tpl) ;
    PyTuple_SetItem(args, 2, stpl) ;
    PyTuple_SetItem(args, 3, mval) ;
    PyTuple_SetItem(args, 4, oval) ;
    if( (result = atomPut(NULL, args)) == NULL) {
      Py_DECREF(args) ; READERR ;
    }
    Py_DECREF(args) ;
    Py_DECREF(result) ;

    if( (dwargs = PyTuple_New(2)) == NULL ) READERR ;
    name2 = Py_BuildValue("s", wordlist.words[0].word) ;
    PyTuple_SetItem(dwargs, 0, name2) ;
    if ( wordlist.nwords > 5 && wordlist.isnum[5] ) {
      /* convert Biso to Uiso */
      Uiso = wordlist.nums[5]/8./3.14159/3.14159 ;
      PyTuple_SetItem(dwargs, 1, Py_BuildValue("d", Uiso)) ;
    } else {
      PyTuple_SetItem(dwargs, 1, Py_BuildValue("d", 0.)) ;
    }
    /*
      NB this will set DW for all atoms of a given element
      ie not allowing for different site symmetries
      to do that would require Copying atomdefs to create diff sites
    */
    if( (result = atomDefDW(NULL, dwargs)) == NULL) {
      Py_DECREF(dwargs) ; READERR ;
    }
    Py_DECREF(dwargs) ;
    Py_DECREF(result) ;

    /* fullprof uses mult*occ */
    if ( wordlist.nwords > 6 && wordlist.isnum[6] ) {
      occ = wordlist.nums[6] ;
      /* FIX this as atomList no longer provides multiplicity directly */
      tpl = Py_BuildValue("(i)", nd) ;
      if( (alst = atomList(NULL, tpl)) == NULL ) { Py_DECREF(tpl) ; READERR ; }
      Py_DECREF(tpl) ;
      if( (alist = PyList_GetItem(alst, 0)) == NULL ) {
	Py_DECREF(alst) ; READERR ;
      }
      imult = PyList_Size(alist) ;
      if( imult < 1 ) {
	PyErr_SetString(PyExc_ValueError, "atomList returned empty list") ;
	Py_DECREF(alst) ; READERR ;
      }
      Py_DECREF(alst) ;
      mult = (double)imult ;
      if( mult > 1. ) occ /= mult ;
      oval = Py_BuildValue("(sd)", wordlist.words[0].word, occ);
      if( (result = editOcc(NULL, oval)) == NULL ) {
	Py_DECREF(oval) ; READERR ;
      }
      Py_DECREF(oval) ;
      Py_DECREF(result) ;
      if( PyErr_Occurred() ) {
	printf("readFULL err in multocc\n") ;
	PyErr_Clear() ;
      }
    }
    if( ! cel ) acel = 0 ;
    if( ! grp ) agrp = 0 ;
    /* read the codes line */
    if( ! fgets(buf, 511, fpt) ) READRET
    nd++ ;
  }


  /* look for the lattice line which has alpha beta ... */
  while(1) {
    if( ! fgets(buf, 511, fpt) ) READRET ;
    if( strstr(buf, "alpha") && strstr(buf, "beta") ) {
      if( ! fgets(buf, 511, fpt) ) READRET ;
      if ( lineparse(buf, &wordlist, 0) < 0 ) READERR ;
      if ( wordlist.nwords > 0 ) {
	tpl = tupleFromDbls(wordlist.nwords, wordlist.nums) ;
	if( (result = lattice(NULL, tpl)) == NULL ) {
	  Py_DECREF(tpl) ; READERR ;
	}
	Py_DECREF(tpl) ;
	Py_DECREF(result) ;
      }
      cel = 1 ;
      break ;
    }
  }

  fclose(fpt) ;
  if( (cel && !acel) || (grp && !agrp) ) badorder = 1 ;
  return 1 ;
}

static void parseCIFlabels(WordList *lbls, int *ic)
{
  int i ;
  for( i=0 ; i<6 ; i++ ) ic[i] = -1 ;
  for( i=0 ; i<lbls->nwords ; i++ ) {
    if( strstr(lbls->words[i].word, "_site_label") ) { ic[0] = i ;
    } else if( strstr(lbls->words[i].word, "_site_fract_x") ) { ic[1] = i ;
    } else if( strstr(lbls->words[i].word, "_site_fract_y") ) { ic[2] = i ;
    } else if( strstr(lbls->words[i].word, "_site_fract_z") ) { ic[3] = i ;
    } else if( strstr(lbls->words[i].word, "_site_occupan") ) { ic[4] = i ;
    } else if( strstr(lbls->words[i].word, "_site_B_iso") )   { ic[5] = i ;
    }
  }
}

static int addWord( char *txt, int ic, int ie, WordList *wordlist ) ;
static int readCIF(char *fil, int ph)
{
  /* can use phase number here ? */
  PyObject *tpl, *stpl, *mval, *oval ;
  PyObject *args, *dwargs, *name, *name2, *result ;
  PyObject *anisList, *editargs ;
  FILE *fpt ;
  int i, nl, nb, nd, ready, ispc, iset, phchk, atnum, Aatnum ;
  int labels, len, ic[6], loopstart ;
  int cel, grp, acel, agrp ;
  char styp, buf[512], namebuf[64], ibuf[8] ;
  char *str, *end, *first ;
  double dbls[6], r[3], s[3], Uiso ;
  double a[6], m[3], mmag ;
  static WordList lbls = { NULL, NULL, NULL, 0, 0, 0 } ;

  cel = grp = 0 ;
  acel = agrp = 1 ;
  s[0] = s[1] = s[2] = 0. ;
  a[3] = a[4] = a[5] = 90. ;


  if( (fpt = fopen(fil, "r")) == NULL ) {
    PyErr_SetString(PyExc_ValueError, "readCIF failed to open file") ;
    return -1 ;
  }
  clearAtomList() ;
  ready = 0 ;
  nl = 0 ;
  while(1) {
    if( ! ready ) if( ! fgets(buf, 511, fpt) ) break ;
    ready = 0 ;
    str = buf ;
    while( *str == ' ' ) str++ ;

    if ( lineparse(str, &wordlist, 0) < 0 ) READERR ;
    if ( wordlist.nwords < 1 ) continue ;
    first = wordlist.words[0].word ;

    if( ! strncmp(first, "data_", 5) ) {
      str = strstr(buf, "data_") ;
      str += 5 ;
      if( strlen(str) > 511 ) str[511] = '\0' ;
      tpl = Py_BuildValue("(s)", str) ;
      if( (result = title(NULL, tpl)) == NULL ) { Py_DECREF(tpl) ; READERR ; }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
    } else if( ! strncmp(first, "_cell_length", 12) ) {
      if ( wordlist.nwords < 2 ) continue ;
      i = first[13] - 97 ;
      if( i < 0 || i > 2 ) continue ;
      a[i] = wordlist.nums[1] ;
      nl += i ;
      if( nl < 3 ) continue ;
      tpl = tupleFromDbls(6, a) ;
      if( (result = lattice(NULL, tpl)) == NULL ) {
	Py_DECREF(tpl) ; READERR ;
      }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
      cel = 1 ;
    } else if( ! strncmp(first, "_cell_angle", 11) ) {
      if ( wordlist.nwords < 2 ) continue ;
      i = first[13] - 97 ;
      if( i == 6 ) i = 2 ;
      if( i < 0 || i > 2 ) continue ;
      a[i+3] = wordlist.nums[1] ;
      tpl = tupleFromDbls(6, a) ;
      if( (result = lattice(NULL, tpl)) == NULL ) {
	Py_DECREF(tpl) ; READERR ;
      }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
    } else if( ! strcmp(first, "_symmetry_space_group_name_H-M") ) {
      if ( wordlist.nwords < 2 ) continue ;
      while( (end = strchr(buf, '\'')) ) *end = ' ' ;
      while( (end = strchr(buf, '\"')) ) *end = ' ' ;
      if ( (str = strstr(buf, "_symmetry")) == NULL ) continue ;
      while( *str != ' ' && *str != '\0' ) str++ ;
      tpl = Py_BuildValue("(s)", str) ;
      if( (result = spacegroup(NULL, tpl)) == NULL) {
	Py_DECREF(tpl) ; READERR ;
      }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
      grp = 1 ;
    } else if( ! strncmp(first, "loop_", 5) ) {
      /* first is list of the datatypes followed by atom data */
      labels = 1 ;
      lbls.nwords = lbls.nnum = 0 ;
      loopstart = 1 ;
      while(1) {
	if( ! fgets(buf, 511, fpt) ) READRET ;
	str = stripstring(buf) ;
	if( (len = strlen(str)) < 1 ) break ; /* end of CIF block */
	if( labels ) {
	  if( strncmp(str, "_atom", 5) ) {
	    if( loopstart ) break ; /* must be a loop with first line _atom */
	    labels = 0 ;
	    parseCIFlabels(&lbls, ic) ;
	    if( ic[0] < 0 || ic[1] < 0 || ic[2] < 0 || ic[3] < 0 ) break;
	  } else {
	    addWord(str, 0, len, &lbls) ;
	    loopstart = 0 ;
	    continue ;
	  }
	}
	if ( lineparse(str, &wordlist, 0) < 0 ) READERR ;
	if ( wordlist.nwords < lbls.nwords ) continue ;
	if( atomNumLookup(wordlist.words[ic[0]].word) <= 0 ) continue ;
	strcpy(namebuf, wordlist.words[ic[0]].word) ;

	for( i=0 ; i<3 ; i++ ) r[i] = wordlist.nums[ic[i+1]] ;
	if( (tpl = tupleFromDbls(3, r)) == NULL ) READERR ;
	if( (stpl = tupleFromDbls(3, s)) == NULL ) {
	  Py_DECREF(tpl) ; READERR ;
	}
	if ( ic[4] >= 0 ) oval = Py_BuildValue("d", wordlist.nums[ic[4]]);
	else oval = Py_BuildValue("d", 1.) ;
	mval = Py_BuildValue("d", 0.) ;
	name = Py_BuildValue("s", namebuf) ;
	if( (args = PyTuple_New(5)) == NULL ) {
	  Py_DECREF(stpl) ; Py_DECREF(tpl) ; Py_DECREF(name) ;
	  Py_DECREF(mval) ; Py_DECREF(oval) ; READERR ;
	}
	PyTuple_SetItem(args, 0, name) ;
	PyTuple_SetItem(args, 1, tpl) ;
	PyTuple_SetItem(args, 2, stpl) ;
	PyTuple_SetItem(args, 3, mval) ;
	PyTuple_SetItem(args, 4, oval) ;
	if( (result = atomPut(NULL, args)) == NULL) {
	  Py_DECREF(args) ; READERR ;
	}
	Py_DECREF(args) ;
	Py_DECREF(result) ;
	if( ! cel ) acel = 0 ;
	if( ! grp ) agrp = 0 ;

	if( (dwargs = PyTuple_New(2)) == NULL ) READERR ;
	name2 = Py_BuildValue("s", namebuf) ;
	PyTuple_SetItem(dwargs, 0, name2) ;
	Uiso = 0. ;
	if ( ic[5] >= 0 ) {
	  /* CIF uses Biso */
	  Uiso = wordlist.nums[ic[5]]/8./3/14159/3.14159 ;
	}
	PyTuple_SetItem(dwargs, 1, Py_BuildValue("d", Uiso)) ;
	/*
	  NB this will set DW for all atoms of a given element
	  ie not allowing for different site symmetries
	  to do that would require Copying atomdefs to create diff sites
	*/
	if( (result = atomDefDW(NULL, dwargs)) == NULL) {
	  Py_DECREF(dwargs) ; READERR ;
	}
	Py_DECREF(dwargs) ;
	Py_DECREF(result) ;
      }
    }
  }

  fclose(fpt) ;
  if( (cel && !acel) || (grp && !agrp) ) badorder = 1 ;
  return 1 ;
}

static int readGSAS(char *fil, int ph)
{
  /* can use phase number here */
  PyObject *tpl, *stpl, *mval, *oval ;
  PyObject *args, *dwargs, *name, *name2, *result ;
  PyObject *anisList, *editargs ;
  FILE *fpt ;
  int i, nb, nd, ready, ispc, iset, phchk, atnum, Aatnum ;
  int cel, grp, acel, agrp ;
  char styp, buf[512], namebuf[64], ibuf[8] ;
  char *str, *end, *first ;
  double dbls[6], s[3], Uiso ;
  double a[6], m[3], mmag ;

  cel = grp = 0 ;
  acel = agrp = 1 ;
  s[0] = s[1] = s[2] = 0. ;
  a[3] = a[4] = a[5] = 90. ;

  if( (fpt = fopen(fil, "r")) == NULL ) {
    PyErr_SetString(PyExc_ValueError, "readGSAS failed to open file") ;
    return -1 ;
  }
  clearAtomList() ;
  ready = 0 ;
  while(1) {
    if( ! ready ) if( ! fgets(buf, 511, fpt) ) break ;
    ready = 0 ;
    str = buf ;
    while( *str == ' ' ) str++ ;

    /* look for CRS lines (crystal stucture) */
    if( strncmp(str, "CRS", 3) ) continue ;
    if( ph >= 0 && sscanf(str+3, "%d", &phchk) > 0 && phchk != ph ) continue ;
    else if( ph < 0 && sscanf(str+3, "%d", &phchk) > 0 ) ph = phchk ;
    /* use the first phase found if none is specified */
    while( *str != ' ' && *str != '\0' ) str++ ;
    if ( lineparse(str, &wordlist, 0) < 0 ) READERR ;
    if ( wordlist.nwords < 2 ) continue ;
    first = wordlist.words[0].word ;
    if( ! strcmp(first, "PNAM") ) {
      str = strstr(buf, "PNAM") ;
      while( *str != ' ' && *str != '\0' ) str++ ;
      if( strlen(str) > 511 ) str[511] = '\0' ;
      tpl = Py_BuildValue("(s)", str) ;
      if( (result = title(NULL, tpl)) == NULL ) { Py_DECREF(tpl) ; READERR ; }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
    } else if( ! strcmp(first, "ABC") ) {
      if ( wordlist.nwords < 4 ) continue ;
      for( i=0 ; i<3 ; i++ ) a[i] = wordlist.nums[i+1] ;
      tpl = tupleFromDbls(6, a) ;
      if( (result = lattice(NULL, tpl)) == NULL ) {
	Py_DECREF(tpl) ; READERR ;
      }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
      cel = 1 ;
    } else if( ! strcmp(first, "ANGLES") ) {
      if ( wordlist.nwords < 4 ) continue ;
      for( i=0 ; i<3 ; i++ ) a[i+3] = wordlist.nums[i+1] ;
      tpl = tupleFromDbls(6, a) ;
      if( (result = lattice(NULL, tpl)) == NULL ) {
	Py_DECREF(tpl) ; READERR ;
      }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
    } else if( ! strcmp(first, "SG") ) {
      if ( (str = strstr(buf, "SYM")) == NULL ) continue ;
      while( *str != ' ' && *str != '\0' ) str++ ;
      tpl = Py_BuildValue("(s)", str) ;
      if( (result = spacegroup(NULL, tpl)) == NULL) {
	Py_DECREF(tpl) ; READERR ;
      }
      Py_DECREF(tpl) ;
      Py_DECREF(result) ;
      grp = 1 ;
    } else if( ! strcmp(first, "AT") ) {
      /*
	atom
	parse after AT which should be
	nnnA  for atom number nnn and then  atomName x y z occLBL
	nnnB  for atom number nnn and then  Uiso I  OR  U11 U22 U33 U12 U13 U23
	nnnM  for atom number nnn and then  sx sy sz  =moment components along
	a b and c direct lattice
      */
      if ( (str = strstr(buf, "AT")) == NULL ) continue ;
      str += 2 ;
      if( strlen(str) < 8 ) continue ;
      if( sscanf(str, "%i%c", &atnum, &styp) < 2 ) continue ;
      while( *str == ' ' && *str != '\0' ) str++ ;
      while( *str != ' ' && *str != '\0' ) str++ ;
      /* this gets str to next word after atomNumberLineType */
      if( styp == 'A' ) {
	if ( lineparse(str, &wordlist, 0) < 0 ) READERR ;
	/* word 1 is atomName, word 2-4 coords, 5 occ */
	if ( wordlist.nwords < 4 ) continue ;
	/* append the atom number to the name */
	if( atomNumLookup(wordlist.words[0].word) <= 0 ) continue ;
	strcpy(namebuf, wordlist.words[0].word) ;
	sprintf(ibuf, "%d", atnum) ;
	strcat(namebuf, ibuf) ;
	if( (tpl = tupleFromDbls(3, wordlist.nums+1)) == NULL ) READERR ;
	if( (stpl = tupleFromDbls(3, s)) == NULL ) {
	  Py_DECREF(tpl) ; READERR ;
	}
	if ( wordlist.nwords > 5 ) oval = Py_BuildValue("d", wordlist.nums[5]);
	else oval = Py_BuildValue("d", 1.) ;
	mval = Py_BuildValue("d", 0.) ;
	name = Py_BuildValue("s", namebuf) ;
	if( (args = PyTuple_New(5)) == NULL ) {
	  Py_DECREF(stpl) ; Py_DECREF(tpl) ; Py_DECREF(name) ;
	  Py_DECREF(mval) ; Py_DECREF(oval) ; READERR ;
	}
	PyTuple_SetItem(args, 0, name) ;
	PyTuple_SetItem(args, 1, tpl) ;
	PyTuple_SetItem(args, 2, stpl) ;
	PyTuple_SetItem(args, 3, mval) ;
	PyTuple_SetItem(args, 4, oval) ;
	if( (result = atomPut(NULL, args)) == NULL) {
	  Py_DECREF(args) ; READERR ;
	}
	Py_DECREF(args) ;
	Py_DECREF(result) ;
	Aatnum = atnum ;
	if( ! cel ) acel = 0 ;
	if( ! grp ) agrp = 0 ;
      } else if( styp == 'B' ) {
	if ( atnum != Aatnum ) continue ;
	if ( lineparse(str, &wordlist, 0) < 0 ) READERR ;
	if ( wordlist.nwords < 2 ) continue ;
	if ( wordlist.nwords < 6 ) {
	  /* GSAS uses Uiso */
	  if( (dwargs = PyTuple_New(2)) == NULL ) READERR ;
	  name2 = Py_BuildValue("s", namebuf) ;
	  Uiso = wordlist.nums[0] ;
	  PyTuple_SetItem(dwargs, 0, name2) ;
	  PyTuple_SetItem(dwargs, 1, Py_BuildValue("d", Uiso)) ;
	} else if(wordlist.nwords >= 6) {
	  if( (dwargs = PyTuple_New(2)) == NULL ) READERR ;
	  if( (anisList = listFromDbls(6, wordlist.nums)) == NULL ) {
	    Py_DECREF(dwargs) ; READERR ;
	  }
	  name2 = Py_BuildValue("s", namebuf) ;
	  PyTuple_SetItem(dwargs, 0, name2) ;
	  PyTuple_SetItem(dwargs, 1, anisList) ;
	}
	/*
	  NB this will set DW for all atoms of a given element
	  ie not allowing for different site symmetries
	  to do that would require Copying atomdefs to create diff sites
	*/
	if( (result = atomDefDW(NULL, dwargs)) == NULL) {
	  Py_DECREF(dwargs) ; READERR ;
	}
	Py_DECREF(dwargs) ;
	Py_DECREF(result) ;
      } else if( styp == 'M' ) {
	/*
	  sx sy sz  here we can use the atom number to editDir editMom
	  NB GSAS specifies mom mag thru sx sy sz while this program
	  uses sx sy sz for direction only and moment is independent
	*/
	if ( lineparse(str, &wordlist, 0) < 0 ) READERR ;
	if ( wordlist.nnum < 3 || ! wordlist.isnum[0] || ! wordlist.isnum[1]
	     || ! wordlist.isnum[2] ) continue ;
	for( i=0 ; i<3 ; i++ ) s[i] = wordlist.nums[i] ;
	/* normalize the s and get the moment */
	mmag = 0. ;
	for( i=0 ; i<3 ; i++ ) {
	  m[i] = s[0]*auni[i] + s[1]*buni[i] + s[2]*cuni[i] ;
	  mmag += m[i]*m[i] ;
	}
	if( mmag > 0. ) {
	  mmag = sqrt(mmag) ;
	  for( i=0 ; i<3 ; i++ ) s[i] /= mmag ; 
	}
	if( (editargs = PyTuple_New(2)) == NULL ) READERR ;
	PyTuple_SetItem(editargs, 0, Py_BuildValue("i", atnum)) ;
	PyTuple_SetItem(editargs, 1, Py_BuildValue("d", mmag)) ;
	if( (result = editMom(NULL, editargs)) == NULL) {
	  Py_DECREF(editargs) ; READERR ;
	}
	Py_DECREF(editargs) ;
	Py_DECREF(result) ;
	if( (editargs = PyTuple_New(4)) == NULL ) READERR ;
	PyTuple_SetItem(editargs, 0, Py_BuildValue("i", atnum)) ;
	PyTuple_SetItem(editargs, 1, Py_BuildValue("d", s[0])) ;
	PyTuple_SetItem(editargs, 2, Py_BuildValue("d", s[1])) ;
	PyTuple_SetItem(editargs, 3, Py_BuildValue("d", s[2])) ;
	if( (result = editDir(NULL, editargs)) == NULL ) {
	  Py_DECREF(editargs) ; READERR ;
	}
	Py_DECREF(editargs) ;
	Py_DECREF(result) ;
      }
    }
  }

  fclose(fpt) ;
  if( (cel && !acel) || (grp && !agrp) ) badorder = 1 ;
  return 1 ;
}

static int readSHELX(char *fil, int ph)
{
  PyObject *tpl, *stpl, *mval, *oval, *args, *dwargs, *name, *name2 ;
  FILE *fpt ;
  int nb, nd, ready, ispc, iset ;
  int cel, grp, acel, agrp ;
  char buf[512] ;
  char *str, *end ;
  double dbls[6], s[3], occ, Uiso ;

  cel = grp = 0 ;
  acel = agrp = 1 ;
  s[0] = s[1] = s[2] = 0. ;
  if( (stpl = tupleFromDbls(3, s)) == NULL ) return -1 ;
  mval = Py_BuildValue("d", 0.) ;
  if( (args = PyTuple_New(5)) == NULL ) return -1 ;
  PyTuple_SetItem(args, 2, stpl) ;
  PyTuple_SetItem(args, 3, mval) ;
  if( (dwargs = PyTuple_New(2)) == NULL ) return -1 ;

  if( (fpt = fopen(fil, "r")) == NULL ) {
    PyErr_SetString(PyExc_ValueError, "readSHELX failed to open file") ;
    return -1 ;
  }
  clearAtomList() ;
  ready = 0 ;
  while(1) {
    if( ! ready ) if( ! fgets(buf, 511, fpt) ) break ;
    ready = 0 ;
    str = buf ;
    while( *str == ' ' ) str++ ;

    if( ! strncmp(str, "TITL", 4) ) {
      if( strlen(str) > 511 ) str[511] = '\0' ;
      if( title(NULL, Py_BuildValue("(s)", str+4)) == NULL ) READERR ;
    } else if( ! strncmp(str, "CELL", 4) ) {
      /* first number after CELL is the wavelength which we ignore */
      if ( lineparse(str+4, &wordlist, 0) < 0 ) READERR ;
      if ( wordlist.nwords > 1 ) {
	if( lattice(NULL, tupleFromDbls(wordlist.nwords-1, wordlist.nums+1))
	    == NULL ) READERR ;
      }
      cel = 1 ;
    } else if( ! strncmp(str, "SYMM", 4) ) {
      str += 4 ;
      if( (end = strchr(str, '(')) != NULL ) *end = '\0' ;
      if(spacegroup(NULL, Py_BuildValue("(s)", str)) == NULL) READERR ;
      grp = 1 ;
    } else if( ! strncmp(str, "UNIT", 4) ) {
      continue ;
    } else {
      /* try atom */
      if ( lineparse(str, &wordlist, 0) < 0 ) READERR ;
      /* word 1 is atomName, 2 atomTypeIndex  word 3-5 coords, 6 occ 7 Uiso */
      if ( wordlist.nwords < 5 ) continue ;
      name = Py_BuildValue("s", wordlist.words[0].word) ;
      name2 = Py_BuildValue("s", wordlist.words[0].word) ;
      if( atomNumLookup(wordlist.words[0].word) <= 0 ) continue ;
      if( ! wordlist.isnum[2] ) continue ;
      if( ! wordlist.isnum[3] ) continue ;
      if( ! wordlist.isnum[4] ) continue ;
      if( (tpl = tupleFromDbls(3, wordlist.nums+2)) == NULL ) READERR ;
      oval = Py_BuildValue("d", 1.) ;
      if ( wordlist.nwords > 5 ) {
	/* I cant figure why shelx convolutes occ like this */
	occ = wordlist.nums[5] ;
	occ = wordlist.nums[5] - (int)occ ;
	occ *= 2. ;
	if( occ <= 0. ) occ = 1. ;
	oval = Py_BuildValue("d", occ);
      }
      PyTuple_SetItem(args, 0, name) ;
      PyTuple_SetItem(args, 1, tpl) ;
      PyTuple_SetItem(args, 4, oval) ;
      PyTuple_SetItem(dwargs, 0, name2) ;
      if(atomPut(NULL, args) == NULL) READERR ;
      if ( wordlist.nwords > 6 && wordlist.isnum[6] ) {
	/* convert Biso to Uiso */
	Uiso = wordlist.nums[6] ;
	PyTuple_SetItem(dwargs, 1, Py_BuildValue("d", Uiso)) ;
      } else {
	PyTuple_SetItem(dwargs, 1, mval) ;
      }
      /*
	NB this will set DW for all atoms of a given element
	ie not allowing for different site symmetries
	to do that would require Copying atomdefs to create diff sites
      */
      if(atomDefDW(NULL, dwargs) == NULL) READERR ;
      if( ! cel ) acel = 0 ;
      if( ! grp ) agrp = 0 ;
    }
  }
  fclose(fpt) ;
  if( (cel && !acel) || (grp && !agrp) ) badorder = 1 ;
  return 1 ;
}
static int readDRAWXTL(char *fil, int ph)
{
  PyObject *tpl, *stpl, *mval, *oval, *args, *dwargs, *name, *name2 ;
  FILE *fpt ;
  int nb, nd, ready ;
  int cel, grp, acel, agrp ;
  char buf[512] ;
  char *str, *end ;
  double dbls[6], s[3], Uiso ;

  cel = grp = 0 ;
  acel = agrp = 1 ;
  s[0] = s[1] = s[2] = 0. ;
  if( (stpl = tupleFromDbls(3, s)) == NULL ) return -1 ;
  mval = Py_BuildValue("d", 0.) ;
  if( (args = PyTuple_New(5)) == NULL ) return -1 ;
  PyTuple_SetItem(args, 2, stpl) ;
  PyTuple_SetItem(args, 3, mval) ;
  if( (dwargs = PyTuple_New(2)) == NULL ) return -1 ;

  if( (fpt = fopen(fil, "r")) == NULL ) {
    PyErr_SetString(PyExc_ValueError, "readDRAWXTL failed to open file") ;
    return -1 ;
  }
  clearAtomList() ;
  ready = 0 ;
  while(1) {
    if( ! ready ) if( ! fgets(buf, 511, fpt) ) break ;
    ready = 0 ;
    str = buf ;
    while( *str == ' ' ) str++ ;

    if( ! strncmp(str, "titl", 4) ) {
      while( *str != ' ' && *str != '\0' ) str++ ;
      if( strlen(str) > 511 ) str[511] = '\0' ;
      if( title(NULL, Py_BuildValue("(s)", str)) == NULL ) READERR ;
    } else if( ! strcmp(str, "cell") ) {
      if ( lineparse(str+4, &wordlist, 0) < 0 ) READERR ;
      if ( wordlist.nwords > 0 ) {
	if( lattice(NULL, tupleFromDbls(wordlist.nwords, wordlist.nums))
	    == NULL ) READERR ;
      }
      cel = 1 ;
    } else if( ! strcmp(str, "spgp") || ! strcmp(str, "spgr")
	       || ! strcmp(str, "sgrp") ) {
      if(spacegroup(NULL, Py_BuildValue("(s)", str+4)) == NULL) READERR ;
      grp = 1 ;
    } else if( ! strcmp(str, "atom") ) {
      if ( lineparse(str+4, &wordlist, 0) < 0 ) READERR ;
      /* word 1 is atomName, word 2 number, word 3-5 coords */
      if ( wordlist.nwords < 5 ) continue ;
      name = Py_BuildValue("s", wordlist.words[0].word) ;
      name2 = Py_BuildValue("s", wordlist.words[0].word) ;
      if( atomNumLookup(wordlist.words[0].word) <= 0 ) continue ;
      if( ! wordlist.isnum[1] ) continue ;
      if( ! wordlist.isnum[2] ) continue ;
      if( ! wordlist.isnum[3] ) continue ;
      if( (tpl = tupleFromDbls(3, wordlist.nums+1)) == NULL ) READERR ;
      oval = Py_BuildValue("d", 1.) ;
      //if ( wordlist.nwords > 5 ) oval = Py_BuildValue("d", wordlist.nums[5]);
      PyTuple_SetItem(args, 0, name) ;
      PyTuple_SetItem(args, 1, tpl) ;
      PyTuple_SetItem(args, 4, oval) ;
      PyTuple_SetItem(dwargs, 0, name2) ;
      if(atomPut(NULL, args) == NULL) READERR ;
      PyTuple_SetItem(dwargs, 1, mval) ;
      /*
	NB this will set DW for all atoms of a given element
	ie not allowing for different site symmetries
	to do that would require Copying atomdefs to create diff sites
      */
      if(atomDefDW(NULL, dwargs) == NULL) READERR ;
      if( ! cel ) acel = 0 ;
      if( ! grp ) agrp = 0 ;
    }
    /*
      could also read arrow commands and convert to mag moment sx sy sz
      although DRAWxtal uses sx along direct-a sy along b* and sz perp to those
      i.e. orthogonal coord sys for sx,xy,sz
    */
  }
  fclose(fpt) ;
  if( (cel && !acel) || (grp && !agrp) ) badorder = 1 ;
  return 1 ;
}
static int readPYSQ(char *fil, int ph)
{
  /*
    NB dont use import which requires .py files
    and may need import statements within the file
    and only executes statements on first read
  */
  FILE *fpt ;
  char *str ;
  char buf[4096] ;
  if( (fpt = fopen(fil, "r")) == NULL ) {
    PyErr_SetString(PyExc_ValueError, "readPYSQ failed to open file") ;
    return -1 ;
  }
  clearAtomList() ;
  while(1) {
    if( ! fgets(buf, 4095, fpt) ) break ;
    if( buf[0] == '#' ) continue ;
    str = stripstring(buf) ;
    /* just pass this to Python interpreter using pysq dictionary */
    PyRun_String(str, Py_eval_input, pysqdict, pysqdict) ;
    if( PyErr_Occurred() ) return -1 ;
  }
  return 1 ;
}


static int readAtomDefFile(char *fil)
{
  PyObject *value ;
  FILE *fpt ;
  char buf[512] ;
  char prefix[8], name[32] ;
  double br, bi, m ;
  int i, atomnum, natom ;
  int ns, new ;
  int istat ;
  double bc[9] ;

  ATOMdef *atomDefPtr ;

  if( (fpt = fopen(fil, "r")) == NULL ) {
    PyErr_SetString(PyExc_IOError, "failed to open file for atom def") ;
    return -1 ;
  }
  natom = 0 ;
  while( fgets(buf, 512, fpt) ) {
    if( buf[0] == '#' ) continue ;

    if( strncmp(buf, "atomdef", 7) == 0 ) {
      value = PyRun_String(buf, Py_eval_input, pysqdict, pysqdict) ;
      if( PyErr_Occurred() != NULL ) PyErr_Clear() ;
      continue ;
    }

    ns = sscanf(buf, "%s %s %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		prefix, name, &br, &bi, &m,
		bc, bc+1, bc+2, bc+3, bc+4, bc+5, bc+6, bc+7, bc+8) ;
    if( ns < 3 ) continue ;
    if( ! isdigit(buf[0]) ) {
      atomnum = atomNumLookup(name) ;
    } else {
      atomnum = atoi(prefix) ;
    }

    if( ns < 4 ) bi = 0. ;
    if( ns < 5 ) m = 0. ;

    if ( (istat = atomLookupS(name, &value, &atomDefPtr)) < 0 ) {
      fclose(fpt) ;
      return -1 ;
    }
    if ( istat == 0 && (atomDefPtr = atomNewS(name)) == NULL ) {
      fclose(fpt) ;
      return -1 ;
    }
    if ( istat == 0 ) natom++ ;

    atomDefPtr->atomic = atomnum ;
    atomDefPtr->b.r = br ;
    atomDefPtr->b.i = bi ;
    atomDefPtr->m = m ;
    atomDefPtr->Eres = 0. ;
    atomDefPtr->HW = 0. ;
    for( i=0 ; i<9 ; i++ ) atomDefPtr->bcoef[i] = 0. ;
    if ( ns < 10 ) {
      atomDefPtr->nbcoef = 0 ;
    } else {
      atomDefPtr->nbcoef = ns - 5 ;
      for( i=0 ; i < ns - 5 ; i++ ) atomDefPtr->bcoef[i] = bc[i] ;
      if( ns == 10 ) { atomDefPtr->Eres = bc[3] ; atomDefPtr->HW = bc[4]/2. ; }
    }
  }
  fclose(fpt) ;
  return natom ;
}

static int readAtomDefFileFF(char *fil)
{
  PyObject *value ;
  PyObject *fflst, *newt ;
  PyObject *self ;
  FILE *fpt ;
  char buf[512] ;
  char name[32] ;
  double j2frac, cf[14] ;
  int i, natom, ns ;

  ATOMdef *atomDefPtr ;

  if( (fpt = fopen(fil, "r")) == NULL ) {
    PyErr_SetString(PyExc_IOError, "failed to open file to read atom def FF") ;
    return -1 ;
  }
  natom = 0 ;
  while( fgets(buf, 511, fpt) ) {
    if( buf[0] == '#' ) continue ;

    if( strncmp(buf, "atomdef", 7) == 0 ) {
      value = PyRun_String(buf, Py_eval_input, pysqdict, pysqdict) ;
      if( PyErr_Occurred() != NULL ) PyErr_Clear() ;
      else natom++ ;
      continue ;
    }

    ns = sscanf(buf,
  "%s %lf   %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf",
		name, &j2frac,
		cf,cf+1,cf+2,cf+3,cf+4,cf+5,cf+6,
		cf+7,cf+8,cf+9,cf+10,cf+11,cf+12,cf+13) ;
    if( ns < 3 ) continue ;

    /* build a list then tuple from the data and call atomdefFF */
    if( (fflst = PyList_New(ns)) == NULL ) { fclose(fpt) ; return -1 ; }
    PyList_SET_ITEM(fflst, 0, Py_BuildValue("s", name)) ;
    PyList_SET_ITEM(fflst, 1, Py_BuildValue("d", j2frac)) ;
    for( i=2 ; i<ns ; i++ )
      PyList_SET_ITEM(fflst, i, Py_BuildValue("d", cf[i-2])) ;
    newt = PyList_AsTuple(fflst) ;
    self = NULL ;
    if( ! atomDefFF(self, newt) ) { fclose(fpt) ; return -1 ; }
    natom++ ;
  }
  fclose(fpt) ;
  return natom ;
}


static int writeAtomDefFile(char *fil)
{
  PyObject *key, *value, *alist ;
  FILE *fpt ;
  char *name ;
  int natom ;
  int pos ;
  int i, j, np, row, col ;

  ATOMdef *a ;

  if( (fpt = fopen(fil, "w")) == NULL ) {
    PyErr_SetString(PyExc_IOError, "failed to open file for atom def") ;
    return -1 ;
  }

  natom = 0 ;
  pos = 0 ;
  if ( (alist = PyList_New(0)) == NULL ) return -1 ;
  while (PyDict_Next((PyObject*)atomDEFdict, &pos, &key, &value)) {
    a = (ATOMdef *) PyCObject_AsVoidPtr(value) ;
    if ( a == NULL ) {
      fclose(fpt) ;
      if ( PyErr_Occurred() ) return -1 ;
      PyErr_SetString(PyExc_ValueError,
		      "failed to convert to ATOMdef pointer") ;
      return -1 ;
    }
    name = PyString_AsString(key) ;

    fprintf(fpt, "atomdef('%s',%g,%g,%g)\n", name, a->b.r, a->b.i, a->m) ;
    if( a->ff != NULL ) {
      fprintf(fpt, "atomdefFF('%s',%g,[", name, a->ff->j2frac) ;
      for( i=0 ; i<a->ff->n ; i++ ) {
	if( i == a->ff->n - 1 ) fprintf(fpt, "%g", a->ff->coefs[i]) ;
	else fprintf(fpt, "%g,", a->ff->coefs[i]) ;
      }
      fprintf(fpt, "])\n") ;
    }
    if( a->dw != NULL ) {
      fprintf(fpt, "atomdefDW('%s',", name) ;
      if( a->dw->iso ) np = 1 ;
      else np = 9 ;
      for( i=0 ; i<np ; i++ ) {
	row = i/3 ;
	col = i%3 ;
	if( i == np-1 ) fprintf(fpt, "%g", a->dw->u[row][col]) ;
	else fprintf(fpt, "%g,", a->dw->u[row][col]) ;
      }
      fprintf(fpt, ")\n") ;
    }
    natom++ ;
  }
  fclose(fpt) ;
  return natom ;
}


static void calcAtomAllQ(int add, ATOMgroup *atomgroup, Qlist *qlst)
{
  /*
    add or subtract from phase sums for all Q in current Qlist
   */
  int j ;
  Qs *qs ;

  for( j=0 ; j<qlst->n ; j++ ) {
    qs = qlst->qlist + j ;
    calcAtomOneQ( add, atomgroup, qs ) ;
  }
}



#if defined(HAVE_DRAND48) && defined(NO_DECL_DRAND48)
extern double drand48 _ANSI_ARGS_((void));
#endif
#if defined(HAVE_SRAND48) && defined(NO_DECL_SRAND48)
extern void srand48 _ANSI_ARGS_((long int seed));
#endif

/*
 * Macros for testing floating-point values for certain special cases:
 *
 *	IS_NAN	Test for not-a-number by comparing a value against itself
 *	IF_INF	Test for infinity by comparing against the largest floating
 *		point value.
 */

#define IS_NAN(v) ((v) != (v))

#ifdef DBL_MAX
#   define IS_INF(v) (((v) > DBL_MAX) || ((v) < -DBL_MAX))
#else
#   define IS_INF(v) 0
#endif


#ifdef __STDC__
#endif /* __STDC__ */


/*---- vector functions -----------------------------------------------------*/

static void vecpro( double vec1[], double vec2[], double prod[] )
{
	prod[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1] ;
	prod[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2] ;
	prod[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0] ;
}

static double dotpro( double vec1[], double vec2[] )
{
	return (vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]) ;
}

static double vecmag( double vec1[] )
{
  double msq ;
  msq = dotpro(vec1, vec1) ;
  if( msq > 0. ) return ( sqrt(msq) ) ;
  return 0. ;
}
static double unitvec( double vec1[] )
{
  int i ;
  double msq ;
  msq = dotpro(vec1, vec1) ;
  if( msq > 0. ) {
    msq = sqrt(msq) ;
    for( i=0 ; i<3 ; i++ ) vec1[i] /= msq ;
  }
  return msq ;
}
static void sclpro( double scl, double vec[], double prod[] )
{
	int i ;
	for( i=0 ; i<3 ; ++i ) prod[i] = scl*vec[i] ;
}

static double recvec( double hkl[], double uni[] )
{
	int i ;
	double qm, qm2 ;

	for( i=0 ; i<3 ; ++i )
		uni[i] = hkl[0]*ast[i] + hkl[1]*bst[i] + hkl[2]*cst[i] ;
	qm2 = dotpro(uni,uni) ;
	if( qm2 <= 0. ) return (0.) ;
	qm = vecmag(uni) ;
	for( i=0 ; i<3 ; ++i )
		uni[i] /= qm ;
	return (qm) ;
}

static void matvec( double mat[3][3], double vec[3], double res[3] )
{
  int i, j ;
  for( i=0 ; i<3 ; i++ )
    {
      res[i] = 0. ;
      for( j=0 ; j<3 ; j++ ) res[i] += mat[i][j] * vec[j] ;
    }
  return ;
}

/* spacegroups */


static int copyPts(int n, Subpt *src, Subpt *dest)
{
  int i ;
  char *s ;
  for( i=0 ; i<n ; i++ ) {
    s = dest[i].ptstr ;
    dest[i] = src[i] ;
    dest[i].ptstr = s ;
    if( src[i].ptstr != NULL )
      /* again Subpt */
      if( ! setString(&(dest[i].ptstr),src[i].ptstr) ) return -1 ;
  }
  return 1 ;
}
static double parseCoef(char *cp, char *buf, int deflt)
{
  char *cs ;
  double sgn ;
  int num, den ;
  if( cp < buf ) return (double)deflt ;
  /*
    search backwards as long as we find digits or div sign
    then check the quit char for sign
  */
  cs = cp ;
  while( cs >= buf && (isdigit(*cs) || *cs == '/' || isspace(*cs)) ) cs-- ;
  sgn = 1. ;
  if( cs >= buf && *cs == '-' ) sgn = -1. ;
  cs++ ;
  num = deflt ;
  den = 1 ;
  while( cs <= cp && isspace(*cs) ) cs++ ;
  if( cs <= cp ) num = atoi(cs) ;
  while( cs <= cp && (isdigit(*cs) || isspace(*cs)) ) cs++ ;
  if( cs < cp && *cs == '/' ) {
    cs++ ;
    while( cs <= cp && isspace(*cs) ) cs++ ;
    den = atoi(cs) ;
  }
  return (sgn*(double)num/(double)den) ;
}

static int parseCoord(char *buf, int row, double m[3][3], double t[])
{
  int i ;
  char *cp ;
  for( i=0 ; i<3 ; i++ ) m[row][i] = 0. ;
  t[row] = 0. ;
  if( (cp = strchr(buf, 'x')) != NULL ) m[row][0] = parseCoef(cp-1, buf, 1) ;
  if( (cp = strchr(buf, 'y')) != NULL ) m[row][1] = parseCoef(cp-1, buf, 1) ;
  if( (cp = strchr(buf, 'z')) != NULL ) m[row][2] = parseCoef(cp-1, buf, 1) ;
  cp = buf + strlen(buf) ;
  t[row] = parseCoef(cp-1, buf, 0) ;
  return 1 ;
}

static int parseCoords(char *buf, int *npts, Subpt *pts)
{
  int i ;
  char *cp, *cp2, *pp, *pp2 ;
  char *sp[3] ;
  char csave ;
  Subpt *pt ;

  *npts = 0 ;
  cp = buf ;
  while( *cp != '\n' && *cp != '\0' ) {
    while( *cp != '(' && *cp != '\0' ) cp++ ;
    if( *cp != '(' ) return *npts ;
    cp2 = cp + 1 ;
    while( *cp2 != ')' && *cp2 != '\n' && *cp2 != '\0' ) cp2++ ;
    if( *cp2 != ')' ) return *npts ;
    csave = *(cp2+1) ;
    *(cp2+1) = '\0' ;
    pt = pts + (*npts) ;
    /* the pts already allocated storage for the string */
    strcpy(pt->ptstr, cp) ;
    *(cp2+1) = csave ;
    /* get the three comma sep strings between the parens */
    pp = cp + 1 ;
    while( isspace(*pp) ) pp++ ;
    pp2 = pp + 1 ;
    while( *pp2 != ',' ) pp2++ ;
    *pp2 = '\0' ;
    sp[0] = pp ;
    pp = pp2 + 1 ;
    while( isspace(*pp) ) pp++ ;
    pp2 = pp + 1 ;
    while( *pp2 != ',' ) pp2++ ;
    *pp2 = '\0' ;
    sp[1] = pp ;
    pp = pp2 + 1 ;
    while( isspace(*pp) ) pp++ ;
    pp2 = pp + 1 ;
    while( *pp2 != ')' ) pp2++ ;
    *pp2 = '\0' ;
    sp[2] = pp ;
    for( i=0 ; i<3 ; i++ ) {
      if( ! parseCoord(sp[i], i, pt->m, pt->t) ) return *npts ;
    }
    (*npts)++ ;
    cp = cp2 + 1 ;
  }
  return *npts ;
}
static int readWyck(char *fil, Spcgrp *grps, int Ngrps)
{
  static char location[] = "readWyck" ;
  /* NB assumes Ngrps storage for all groups is preallocated */
  FILE *fpt ;
  char buf[256] ;
  char *cp, *cp2 ;
  int i, j, nset ;
  int intN ;
  int nspc, nsub, nl, npts, np, ntrans ;
  Spcgrp *grp ;
  static Subpt *subpts = NULL ;

  if( (fpt = fopen(fil, "r")) == NULL ) return 0 ;

  /* allocate string storage for our Subpt buffers */
  if( subpts == NULL ) {
    subpts = (Subpt *)calloc(6, sizeof(Subpt)) ;
    if ( memExc(isNULL(subpts), location) ) return -1 ;
    for( i=0 ; i<6 ; i++ ) {
      subpts[i].ptstr = (char*)malloc(64*sizeof(char)) ;
      if ( memExc(isNULL(subpts[i].ptstr), location) ) return -1 ;
    }
  }
  intN = 0 ;
  nspc = 0 ;
  while( fgets(buf, 256, fpt) ) {
    /* first find a line beginning Wyck to start a grp */
    if( strncmp("Wyck", buf, 4) != 0 ) continue ;
    /* get the group number, symbol, settings and origin from this line */
    if( !(cp = strstr(buf, "Group")) ) return 0 ;
    cp += 6 ;
    intN = atoi(cp) ;
    if( intN < 1 || intN > 230 ) return 0 ;
    grp = grps + nspc ;
    grp->spcnum = intN ;
    /* move position past group number */
    while( isspace(*cp) ) cp++ ;
    while( isdigit(*cp) ) cp++ ;
    while( isspace(*cp) ) cp++ ;
    /* move past the opening paren for spcgrp symbol */
    cp++ ;
    i = 0 ;
    while( *cp != ')' && *cp != '\0' && *cp != '\n' ) {
      grp->symbol[i] = *cp ;
      i++ ;
      cp++ ;
    }
    grp->symbol[i] = '\0' ;

    /* move past the closing paren and look for [ */
    if( *cp == ')' ) cp++ ;
    while( *cp != '\0' && *cp != '\n' && *cp != '[' ) cp++ ;

    /* check for up to 6 setting numbers in remaining chars */
    grp->setting[0] = 1 ;
    for(j=1 ; j<6 ; j++) grp->setting[j] = 0 ;
    nset = 0 ;
    while( *cp != '\0' && *cp != '\n' && *cp != ']' ) {
      if( isdigit(*cp) ) {
	grp->setting[nset] = atoi(cp) ;
	nset++ ;
	if( nset > 5 ) break ;
	while( isdigit(*cp) ) cp++ ;
      } else {
	cp++ ;
      }
    }

    /* check for origin inside second set of square brackets */
    grp->origin[0] = '\0' ;
    while( *cp != '\0' && *cp != '\n' && *cp != '[' ) cp++ ;
    if( *cp == '[' ) {
      cp++ ;
      nset = 0 ;
      while( *cp != '\0' && *cp != '\n' && *cp != ']' ) {
	grp->origin[nset] = *cp ;
	cp++ ;
	nset++ ;
	if( nset >= 15 ) break ;
      }
      grp->origin[nset] = '\0' ;
    }

    /* read down past the line which says Coordinates */
    for( i=0 ; i<4 ; i++ ) if( ! fgets(buf, 256, fpt) ) return nspc ;
    /* is this line a list of translation vectors OR start of coords */
    grp->ntrans = 0 ;
    ntrans = 1 ;
    if( ! isdigit(buf[0]) ) {
      /* must be translation vectors line */
      if( parseCoords(buf, &(grp->ntrans), subpts) <= 0 )
	{
	  if( PyErr_Occurred() ) return -1 ;
	  return nspc ;
	}
      if( copyPts(grp->ntrans, subpts, grp->trans) < 0 ) return -1 ;
      ntrans = grp->ntrans ;
      if( ! fgets(buf, 256, fpt) ) return nspc ;
    }
    nsub = 0 ;
    grp->sub = NULL ;
    while( isdigit(buf[0]) ) {
      grp->nsub = nsub + 1 ;
      grp->sub = (Subgrp*)realloc(grp->sub, (nsub+1)*sizeof(Subgrp)) ;
      if ( memExc(isNULL(grp->sub), location) ) return -1 ;
      grp->sub[nsub].npts = atoi(buf)/ntrans ;
      grp->sub[nsub].pts =
	(Subpt *)calloc(grp->sub[nsub].npts, sizeof(Subpt)) ;
      if ( memExc(isNULL(grp->sub[nsub].pts), location) ) return -1 ;
      for( i=0 ; i<grp->sub[nsub].npts ; i++ )
	grp->sub[nsub].pts[i].ptstr = NULL ;
      cp = buf ;
      while( isdigit(*cp) ) cp++ ;
      while( isspace(*cp) ) cp++ ;
      grp->sub[nsub].Wyck = *cp ;
      grp->sub[nsub].index = intN ;
      cp++ ;
      while( isspace(*cp) ) cp++ ;
      cp2 = cp + 1 ;
      while( *cp2 != '\n' && *cp2 != '\0' ) cp2++ ;
      while( cp2 > cp && isspace(*(cp2-1)) ) cp2-- ;
      *cp2 = '\0' ;
      strcpy(grp->sub[nsub].sitesym, cp) ;
      /* calc number of lines of pts assuming up to 4 pts per line */
      nl = 1 + (grp->sub[nsub].npts - 1)/4 ;
      np = 0 ;
      for( i=0 ; i<nl ; i++ ) {
	if( ! fgets(buf, 256, fpt) ) return nspc ;
	if( parseCoords(buf, &npts, subpts) <= 0 )
	  {
	    if( PyErr_Occurred() ) return -1 ;
	    return nspc ;
	  }
	if( npts < 1 ) return nspc ;
	if( copyPts(npts, subpts, grp->sub[nsub].pts+np) < 0 ) return -1 ;
	np += npts ;
      }
      nsub++ ;
      if( ! fgets(buf, 256, fpt) ) return nspc ;
    } /* while isdigit buf[0] symmetry sub group multiplicity */
    nspc++ ;
    /* make sure group records are sep by at least one blank line */
  } /* while fgets line */
  return nspc ;
}

static int readAliases(char *fil, PyDictObject *d, PyDictObject *g)
{

  PyObject *avalue, *gvalue ;
  FILE *fpt ;
  char buf[256] ;
  char alias[128], sep[32], grp[128] ;
  int nalias ;

  if( (fpt = fopen(fil, "r")) == NULL ) return 0 ;

  nalias = 0 ;
  while( fgets(buf, 256, fpt) ) {
    if( buf[0] == '#' || strlen(buf) < 3 ) continue ;
    /* read an alias pair */
    if( sscanf(buf, "'%[^']'%[^']'%[^']", alias, sep, grp) < 2 ) continue ;
    if ((gvalue=PyDict_GetItemString((PyObject *)g, grp)) == NULL) continue ;
    if ((avalue=PyDict_GetItemString((PyObject *)d, alias)) != NULL) continue ;
    gvalue = Py_BuildValue("s", grp) ;
    if ( PyDict_SetItemString((PyObject *)d, alias, gvalue) ) continue ;
    nalias++ ;
  }
  fclose(fpt) ;
  return nalias ;
}


/* complex functions */


static DCMPLX Cadd( DCMPLX a, DCMPLX b )
{
	DCMPLX c ;
	c.r = a.r + b.r ;	c.i = a.i + b.i ;
	return c ;
}

static DCMPLX Csub( DCMPLX a, DCMPLX b )
{
	DCMPLX c ;
	c.r = a.r - b.r ;	c.i = a.i - b.i ;
	return c ;
}

static DCMPLX Cmul( DCMPLX a, DCMPLX b )
{
	DCMPLX c ;
	c.r = a.r*b.r - a.i*b.i ;	c.i = a.i*b.r + a.r*b.i ;
	return c ;
}

static DCMPLX Cdiv( DCMPLX a, DCMPLX b )
{
	DCMPLX c ;
	double r, den ;
	if( fabs(b.r) >= fabs(b.i) )
		{
		r = b.i/b.r ;			den = b.r + r*b.i ;
		c.r = (a.r + r*a.i)/den ;	c.i = (a.i - r*a.r)/den ;
		}
	else
		{
		r = b.r/b.i ;			den = b.i + r*b.r ;
		c.r = (r*a.r + a.i)/den ;	c.i = (r*a.i - a.r)/den ;
		}
	return c ;
}

static DCMPLX Cmplx( double r, double i )
{
	DCMPLX c ;
	c.r = r ;	c.i = i ;
	return c ;
}

static double Cabs( DCMPLX c )
{
	double x, y, m, temp ;
	x = fabs(c.r) ;	y = fabs(c.i) ;
	if	( x == 0.0 )	m = y ;
	else if	( y == 0.0 )	m = x ;
	else if ( x > y )	{ temp = y/x ;	m = x*sqrt(1. + temp*temp) ; }
	else			{ temp = x/y ;	m = y*sqrt(1. + temp*temp) ; }
	return m ;
}

static double Cabsq( DCMPLX c )
{
  return (c.r*c.r + c.i*c.i) ;
}

static DCMPLX Cconj( DCMPLX c )
{
	DCMPLX d ;
	d.r = c.r ;	d.i = -c.i ;
	return d ;
}
static DCMPLX Csqrt( DCMPLX c )
{
	DCMPLX d ;
	double x, y, w, r ;
	if( (c.r == 0.0 ) && (c.i == 0.0) )
		{ d.r = d.i = 0.0 ;	return d ; }
	else
		{
		x = fabs(c.r) ;	y = fabs(c.i) ;
		if( x >= y )
			{ r = y/x ; w = sqrt(x)*sqrt(0.5*(1.+sqrt(1.+r*r))) ; }
		else
			{ r = x/y ; w = sqrt(y)*sqrt(0.5*(r +sqrt(1.+r*r))) ; }
		if( c.r >= 0. )
			{ d.r = w ;	d.i = c.i/(2.*w) ; }
		else
			{ d.i = (c.i >= 0. ? w : -w) ; d.r = c.i/(2.*d.i) ; }
		return d ;
		}
}
static DCMPLX RCmul( double f, DCMPLX c )
{
	DCMPLX d ;
	d.r = f*c.r ;	d.i = f*c.i ;
	return d ;
}

static DCMPLX Cexp( DCMPLX c )
{
	DCMPLX d ;
	double mag ;
	mag = exp(c.r) ;
	d.r = mag*cos(c.i) ;
	d.i = mag*sin(c.i) ;
	return d ;
}
static DCMPLX RCexp( double r )	/* exp(ir) */
{
	DCMPLX d ;
	d.r = cos(r) ;
	d.i = sin(r) ;
	return d ;
}

static DCMPLX CRpow( DCMPLX c, double p )
{
	DCMPLX d ;
	double mag, phase ;
	mag = c.r*c.r + c.i*c.i ;
	phase = p*atan2(c.i, c.r) ;	/* -pi to pi */
	mag = pow(mag, p/2.) ;
	d.r = mag*cos(phase) ;
	d.i = mag*sin(phase) ;
	return d ;
}

static int memExc(int isnull, char *loc)
{
  static char unknown[] = "UNKNOWN location" ;
  if ( ! isnull ) return 0 ;
  if ( ! loc ) PyErr_SetString(PyExc_MemoryError, unknown) ;
  else PyErr_SetString(PyExc_MemoryError, loc) ;
  return 1 ;
}

static char *stralloc(const char *string)
{
  char *copy;
  //if (string == NULL) return NULL ;
  copy = malloc(strlen(string) + 1);
  if (copy == NULL)
    return NULL;
  strcpy(copy, string);
  return copy;
}
//static char *strdup(const char *string)
//{
// char *copy;
//if (string == NULL) return NULL ;
//  copy = malloc(strlen(string) + 1);
//if (copy == NULL)
//  return NULL;
//strcpy(copy, string);
//return copy;
//}

static int setString(char **s, const char *src)
{
  static char location[] = "setString" ;
  if( s == NULL ) return 1 ;
  if( *s != NULL ) { free(*s) ; *s = NULL ; }
  if( src == NULL ) return 1 ;
  *s = stralloc(src) ;
  if ( memExc(isNULL(*s), location) ) return 0 ;
  return 1 ;
}

static int recip()
{
  int i ;
  double dang[3], prod[3], dval    ;

  for( i=0 ; i<3 ; ++i )
    {
      dang[i] = 0. ;
      if( latt[i] <= 0. || angl[i] <= 0. ) return(0) ;
      else    dang[i] = DEGTORAD * angl[i] ;
    }

  /* x-axis along latt[0], bvec in x-y plane */
  /* find components of real space vectors in xtal Cartesian system */

  avec[0] = latt[0] ; avec[1] = 0. ; avec[2] = 0. ;
  bvec[0] = latt[1]*cos(dang[2])   ;
  bvec[1] = latt[1]*sin(dang[2])   ; bvec[2] = 0. ;
  if( bvec[1] == 0. ) return(0) ;
  cvec[0] = latt[2]*cos(dang[1])  ;
  cvec[1] = (latt[2]*latt[1]*cos(dang[0])-cvec[0]*bvec[0])/bvec[1] ;
  dval = latt[2]*latt[2] - cvec[0]*cvec[0] - cvec[1]*cvec[1] ;
  if( dval <= 0. ) return(0) ;
  cvec[2] = sqrt(dval) ;

  vecpro( bvec, cvec, prod ) ;
  dval = dotpro( avec, prod ) ;

  for( i=0 ; i<3 ; ++i ) ast[i] = TWOPI*prod[i]/dval ;
  vecpro( cvec, avec, prod ) ;
  for( i=0 ; i<3 ; ++i ) bst[i] = TWOPI*prod[i]/dval ;
  vecpro( avec, bvec, prod ) ;
  for( i=0 ; i<3 ; ++i ) cst[i] = TWOPI*prod[i]/dval ;

  for( i=0 ; i<3 ; i++ )
    {
      auni[i]=avec[i]/latt[0] ;
      buni[i]=bvec[i]/latt[1] ;
      cuni[i]=cvec[i]/latt[2] ;
    }

  rlat[0] = sqrt(ast[0]*ast[0]+ast[1]*ast[1]+ast[2]*ast[2]) ;
  rlat[1] = sqrt(bst[0]*bst[0]+bst[1]*bst[1]+bst[2]*bst[2]) ;
  rlat[2] = sqrt(cst[0]*cst[0]+cst[1]*cst[1]+cst[2]*cst[2]) ;
  for( i=0 ; i<3 ; ++i ) if(rlat[i] <= 0.) return(0) ;

  rang[0] = acos(dotpro(ast,bst)/(rlat[0]*rlat[1]))/DEGTORAD ;
  rang[1] = acos(dotpro(bst,cst)/(rlat[1]*rlat[2]))/DEGTORAD ;
  rang[2] = acos(dotpro(cst,ast)/(rlat[2]*rlat[0]))/DEGTORAD ;

  return(1) ;
}


static PyObject *MakeDblList(int n, double *dlst)
{
  int j ;
  PyObject *dList ;

  if ( dlst == NULL ) return PyList_New(0) ;
  if ( (dList = PyList_New(n)) == NULL ) return NULL ;

  for( j=0 ; j<n ; j++ ) {
    PyList_SET_ITEM(dList, j, Py_BuildValue("d", dlst[j])) ;
  }
  return dList ;
}


static PyObject *MakeDWList(DebyeWallerFactor *dw)
{
  int i, j ;
  PyObject *dList ;

  if ( dw == NULL ) return PyList_New(0) ;
  if ( dw->iso ) {
    if ( (dList = PyList_New(1)) == NULL ) return NULL ;
    PyList_SET_ITEM(dList, 0, Py_BuildValue("d", dw->u[0][0])) ;
    return dList ;
  }
  if ( (dList = PyList_New(9)) == NULL ) return NULL ;
  for( i=0 ; i<3 ; i++ ) {
    for( j=0 ; j<3 ; j++ )
      PyList_SET_ITEM(dList, 3*i + j, Py_BuildValue("d", dw->u[i][j])) ;
  }
  return dList ;
}


static void rhsetup(double hkl1[3], double hkl2[3], double hkl3[3])
{
  /*
    setup the rhc orientation vectors starting with hkl1 and if hkl2 is not
    parallel use it as the second one defining the scatt plane, else
    try hkl3 to define the scatt plane
    put the results in ip1.qun ip2.qun ip3.qun (ip3 will be out of plane)
   */
  ip1.q = recvec(hkl1, ip1.qun) ;
  ip2.q = recvec(hkl2, ip2.qun) ;
  vecpro( ip1.qun, ip2.qun, ip3.qun ) ;
  if( unitvec(ip3.qun) <= 0. )
    {
      ip2.q = recvec(hkl3, ip2.qun) ;      
      vecpro( ip1.qun, ip2.qun, ip3.qun ) ;
      unitvec(ip3.qun) ;
    }
  vecpro( ip3.qun, ip1.qun, ip2.qun ) ;
}


static void pbcalc()
{
  int i, jsm ;
  double pbmag, psm ;

  if( pbq != 0 ) {
    for( i=0 ; i<3 ; i++ ) pbx[i] = ip3.qun[i] ;
    return ;
  }

  pbmag = vecmag(pbz) ;
  if( pbmag <= 0. ) pbmag = 1. ;
  for( i=0 ; i<3 ; i++ ) pbz[i] /= pbmag ;
  jsm = 0 ;
  psm = fabs(pbz[0]) ;
  if( fabs(pbz[1]) < psm ) { jsm = 1 ; psm = fabs(pbz[1]) ; }
  if( fabs(pbz[2]) < psm ) { jsm = 2 ; psm = fabs(pbz[2]) ; }
  for( i=0 ; i<3 ; i++ ) pbx[i] = 0. ;
  pbx[jsm] = 1. ;
  vecpro(pbx, pbz, pby) ;
  vecpro(pby, pbz, pbx) ;
  /* pbx and pby are the other rhs axes in quantization system */
  
  /* construct the matrix transforms for pbx, pby, pbz so that */
  /* Mx(cosw, sinw, 1) = pbx(xtal coordinates) etc. */
  
  for( i=0 ; i<3 ; i++ )
    {
      trx[i][0] = pbx[0]*ip1.qun[i] + pbx[1]*ip2.qun[i] ;
      trx[i][1] = pbx[1]*ip1.qun[i] - pbx[0]*ip2.qun[i] ;
      trx[i][2] = pbx[2]*ip3.qun[i] ;
      try[i][0] = pby[0]*ip1.qun[i] + pby[1]*ip2.qun[i] ;
      try[i][1] = pby[1]*ip1.qun[i] - pby[0]*ip2.qun[i] ;
      try[i][2] = pby[2]*ip3.qun[i] ;
      trz[i][0] = pbz[0]*ip1.qun[i] + pbz[1]*ip2.qun[i] ;
      trz[i][1] = pbz[1]*ip1.qun[i] - pbz[0]*ip2.qun[i] ;
      trz[i][2] = pbz[2]*ip3.qun[i] ;
    }
}

static int atomNumLookup(char *buf)
{
  int i, nc ;

  nc = 0 ;
  while( (buf[nc]>=65 && buf[nc]<=90) || (buf[nc]>=97 && buf[nc]<=122) ) nc++ ;
  if( nc < 1 ) return 0 ;
  buf[0] = toupper(buf[0]) ;
  if( nc > 2 ) nc = 2 ;
  if( nc > 1 ) buf[1] = tolower(buf[1]) ;
  for( i=0 ; i<Nelem ; i++ ) {
    if( i==0 ) {
      if( nc == 1 && ((buf[0] == 'H') || (buf[0] == 'D') || (buf[0] == 'T')) )
	return 1 ;
      continue ;
    }
    if( strncmp(buf, ElemSymbols[i], nc) == 0 ) return (i+1) ;
  }
  /* on nc=2 failure try nc=1 */
  if( nc < 2 ) return 0 ;
  nc = 1 ;
  for( i=0 ; i<Nelem ; i++ ) {
    if( i==0 ) {
      if( nc == 1 && ((buf[0] == 'H') || (buf[0] == 'D') || (buf[0] == 'T')) )
	return 1 ;
      continue ;
    }
    if( strncmp(buf, ElemSymbols[i], nc) == 0 ) return (i+1) ;
  }
  return 0 ;
}

static int CopyAtomDef(ATOMdef *dest, ATOMdef *src)
{
  /* dont mess with name or description */
  int i ;
  static char location[] = "CopyAtomDef" ;
  if( src == NULL || dest == NULL ) return -1 ;
  if( dest->ff != NULL ) free(dest->ff) ;
  if( dest->dw != NULL ) free(dest->dw) ;
  dest->ff = NULL ; dest->dw = NULL ;
  dest->atomic = src->atomic ;
  dest->b.r = src->b.r ;
  dest->b.i = src->b.i ;
  dest->nbcoef = src->nbcoef ;
  for( i=0 ; i<9 ; i++ ) dest->bcoef[i] = src->bcoef[i] ;
  dest->m = src->m ;
  if( src->ff != NULL ) {
    dest->ff = (MagFormFactor*)calloc(1,sizeof(MagFormFactor)) ;
    if( memExc(isNULL(dest->ff), location) ) return 0 ;
    *dest->ff = *src->ff ;
  }
  if( src->dw != NULL ) {
    dest->dw = (DebyeWallerFactor*)calloc(1,sizeof(DebyeWallerFactor)) ;
    if( memExc(isNULL(dest->dw), location) ) return 0 ;
    *dest->dw = *src->dw ;
  }
  return 1 ;
}

static int atomLookupP(PyObject *key, PyObject **value, ATOMdef **aptr)
{
  if( key == NULL || aptr == NULL ) {
    PyErr_SetString(PyExc_ValueError, "null arg to atomLookupP") ;
    return -1 ;
  }
  *aptr = NULL ;
  if( (*value = PyDict_GetItem((PyObject *)atomDEFdict, key)) == NULL )
    return 0 ;
  *aptr = (ATOMdef *) PyCObject_AsVoidPtr(*value) ;
  if ( *aptr == NULL ) {
      if ( PyErr_Occurred() ) return 0 ;
      PyErr_SetString(PyExc_ValueError,
		      "failed to convert to ATOMdef pointer in atomLookupP") ;
      return 0 ;
  }
  return 1 ;
}
static int atomLookupS(char *key, PyObject **value, ATOMdef **aptr)
{
  if( key == NULL || aptr == NULL ) {
    PyErr_SetString(PyExc_ValueError, "null arg to atomLookupS") ;
    return -1 ;
  }
  *aptr = NULL ;
  if( (*value = PyDict_GetItemString((PyObject *)atomDEFdict, key)) == NULL )
    return 0 ;
  *aptr = (ATOMdef *) PyCObject_AsVoidPtr(*value) ;
  if ( *aptr == NULL ) {
      if ( PyErr_Occurred() ) return -1 ;
      PyErr_SetString(PyExc_ValueError,
		      "failed to convert to ATOMdef pointer in atomLookupP") ;
      return -1 ;
  }
  return 1 ;
}
static ATOMdef *atomNewP(PyObject *key)
{
  PyObject *value ;
  static char location[] = "atomNewP" ;
  char *name ;
  ATOMdef *aptr ;

  if( key == NULL ) {
    PyErr_SetString(PyExc_ValueError, "atomNewP called with NULL arg");
    return NULL ;
  }
  aptr = (ATOMdef *)calloc(1, sizeof(ATOMdef)) ;
  if ( memExc(isNULL(aptr), location) ) return NULL ;
  value = PyCObject_FromVoidPtr((void*)aptr, NULL) ;
  if ( PyDict_SetItem((PyObject *)atomDEFdict, key, value) ) {
    free(aptr) ;
    PyErr_SetString(PyExc_ValueError, "copy failed to put atom in dictionary");
    return NULL ;
  }
  if ( (name = PyString_AsString(key)) == NULL ) return NULL ;
  if( ! setString(&(aptr->name), name) ) return NULL ;
  aptr->atomic = atomNumLookup(name) ;
  ++nAtomDef ;
}
static ATOMdef *atomNewS(char *key)
{
  PyObject *value ;
  static char location[] = "atomNewS" ;
  char *name ;
  ATOMdef *aptr ;

  if( key == NULL ) {
    PyErr_SetString(PyExc_ValueError, "atomNewS called with NULL arg");
    return NULL ;
  }
  aptr = (ATOMdef *)calloc(1, sizeof(ATOMdef)) ;
  if ( memExc(isNULL(aptr), location) ) return NULL ;
  value = PyCObject_FromVoidPtr((void*)aptr, NULL) ;
  if ( PyDict_SetItemString((PyObject *)atomDEFdict, key, value) ) {
    free(aptr) ;
    PyErr_SetString(PyExc_ValueError, "copy failed to put atom in dictionary");
    return NULL ;
  }
  if( ! setString(&(aptr->name), key) ) return NULL ;
  aptr->atomic = atomNumLookup(key) ;
  ++nAtomDef ;
  return aptr ;
}

static int ParseFlagsFromArgs(int nflags, Flag *flags,
			      int narg, PyObject *args)
{
  int i, j, nstr, nchk, flen ;
  char *str ;
  PyObject *arg, *src ;
  /* args must alternate between keyword, int value */
  src = args ;
  if( narg > 0 ) {
    if ( (arg = PyTuple_GetItem(args, 0)) == NULL ) return 0 ;
    if ( PyList_Check(arg) ) {
      src = PyList_AsTuple(arg) ;
      narg = PyTuple_Size(src) ;
    }
  }
  for( i=0 ; i<narg ; i+=2 ) {
    if ( (arg = PyTuple_GetItem(src, i)) == NULL ) return 0 ;
    if ( ! PyString_Check(arg) ) {
      PyErr_SetString(PyExc_ValueError, "ParseFlags failed on keyword") ;
      return 0 ;
    }
    str = PyString_AsString(arg) ;
    nstr = strlen(str) ;
    if( nstr < 1 ) continue ;
    /* find the first match in flags */
    for( j=0 ; j<nflags ; j++ ) {
      nchk = nstr ;
      flen = strlen(flags[j].name) ;
      if( flen < 1 ) break ;
      if( flen < nstr ) nchk = flen ;
      if( ! strncmp(flags[j].name, str, nchk) ) {
	if ( (arg = PyTuple_GetItem(src, i+1)) == NULL ) return 0 ;	
	if ( ! PyInt_Check(arg) ) {
	  PyErr_SetString(PyExc_ValueError, "ParseFlags failed on int flag") ;
	  return 0 ;
	}
	flags[j].i = PyInt_AsLong(arg) ;
	break ;
      }
    }
  }
  return 1 ;
}

static PyObject *MakeFlagList(int n, Flag *flst)
{
  int j ;
  static char noname[] = "NOname" ;
  PyObject *sList ;

  if ( flst == NULL ) return PyList_New(0) ;
  if ( (sList = PyList_New(2*n)) == NULL ) return NULL ;

  for( j=0 ; j<n ; j++ ) {
    if ( flst[j].name == NULL ) {
      PyList_SET_ITEM(sList, 2*j, Py_BuildValue("s", noname)) ;
    } else {
      PyList_SET_ITEM(sList, 2*j, Py_BuildValue("s", flst[j].name)) ;
    }
    PyList_SET_ITEM(sList, 2*j + 1, Py_BuildValue("i", flst[j].i)) ;
  }
  return sList ;
}

static int sameAtomPos(myATOM *a, myATOM *b, double tol)
{
  int i ;
  double diff, mag ;
  mag = 0. ;
  for( i=0 ; i<3 ; i++ ) {
    diff = a->r[i] - b->r[i] ;
    mag += diff*diff ;
  }
  if( mag > 0 ) mag = sqrt(mag) ;
  if( mag <= tol ) return 1 ;
  return 0 ;
}
static int pruneAtoms()
{
  int i, j, k, l, ndel ;
  double toler ;
  ATOMgroup *gi, *gk ;
  myATOM *aj, *al ;

  /*
    for each atomgroups atoms check if any atomgroup
    farther down in the list has an atom with same coordinates
    to within the positionTolerance
  */

  for( i=0 ; i<atomslist.n ; i++ ) atomslist.atomslist[i].omit = 0 ;
  toler = pow(10., (double)(*atomtolerance)) ;
  for( i=0 ; i<atomslist.n - 1 ; i++ ) {
    gi = atomslist.atomslist + i ;
    for( j=0 ; j<gi->natoms ; j++ ) {
      aj = gi->atoms + j ;
      for( k=i+1 ; k<atomslist.n ; k++ ) {
	gk = atomslist.atomslist + k ;
	for( l=0 ; l<gk->natoms ; l++ ) {
	  al = gk->atoms + l ;
	  if( sameAtomPos(aj, al, toler) ) gk->omit = 1 ;
	}
      }
    }
  }
  ndel = omitAtoms() ;
  return ndel ;
}

static int omitAtoms()
{
  int i, ndel, nlist ;
  /*
    now truncate the atomslist
    first delete storage for atomgroups to be deleted
  */

  ndel = 0 ;
  for( i=0 ; i<atomslist.n ; i++ ) {
    if( atomslist.atomslist[i].omit == 1 ) {
      if ( atomslist.atomslist[i].name != NULL ) {
	free(atomslist.atomslist[i].name) ;
	atomslist.atomslist[i].name = NULL ;
      }
      free(atomslist.atomslist[i].atoms) ;
      ndel++ ;
    }
  }
  /* now rebuild the list from ones not omited */
  nlist = 0 ;
  for( i=0 ; i<atomslist.n ; i++ ) {
    if( atomslist.atomslist[i].omit == 1 ) continue ;
    *(atomslist.atomslist+nlist) = *(atomslist.atomslist+i) ;
    nlist++ ;
  }
  /* NULL storage for zombies left at end */
  for( i=nlist ; i<atomslist.n ; i++ ) {
    atomslist.atomslist[i].mult = 0 ;
    atomslist.atomslist[i].natoms = 0 ;
    atomslist.atomslist[i].atoms = NULL ;
    atomslist.atomslist[i].name = NULL ;
  }
  atomslist.n = nlist ;
  return ndel ;
}


static int addWord( char *txt, int ic, int ie, WordList *wordlist )
{
  static char location[] = "addWord" ;
  /*
    grab a word from txt delimited by indices ic ie
    index ie is just after word end in txt
    add the word to wordlist
  */
  int i, j, m, n, na, len, len1 ;
  Word *word ;

  if ( wordlist == NULL ) return 0 ;
  m = wordlist->nwords ;
  n = m + 1 ;
  na = wordlist->nalloc ;
  if( n > na ) {
    wordlist->nalloc = 0 ;
    wordlist->words =
      (Word*)realloc((char*)wordlist->words, n*sizeof(Word)) ;
    if ( memExc(isNULL(wordlist->words), location) ) return -1 ;
    wordlist->nums = (double*)realloc((char*)wordlist->nums, n*sizeof(double));
    if ( memExc(isNULL(wordlist->nums), location) ) return -1 ;
    wordlist->isnum = (int*)realloc((char*)wordlist->isnum, n*sizeof(int)) ;
    if ( memExc(isNULL(wordlist->isnum), location) ) return -1 ;
    for( i=na ; i<n ; i++ ) {
      wordlist->words[i].word = NULL ;
      wordlist->words[i].nalloc = 0 ;
    }
    wordlist->nalloc = n ;
  }
  wordlist->nwords = n ;
  len = ie - ic ;
  len1 = len + 1 ;
  word = wordlist->words + m ;
  if( len1 > word->nalloc ) {
    word->nalloc = 0 ;
    word->word = (char*)realloc((char*)word->word, len1*sizeof(char)) ;
    if ( memExc(isNULL(wordlist->words), location) ) return -1 ;
    word->nalloc = len1 ;
  }
  for( i=ic, j=0 ; i<ie ; i++, j++ ) {
    word->word[j] = txt[i] ;
  }
  word->word[len] = '\0' ;
  word->size = len ;
  return 1 ;
}

static void readNumbers(WordList *wordlist)
{
  int i ;
  wordlist->nnum = 0 ;
  for( i=0 ; i<wordlist->nwords ; i++ ) {
    if( sscanf(wordlist->words[i].word, "%lf", wordlist->nums+i) > 0 ) {
      wordlist->isnum[i] = 1 ;
      wordlist->nnum++ ;
    } else {
      wordlist->isnum[i] = 0 ;
    }
  }
}

static int lineparse(char *line, WordList *wordlist, int parenflag)
{
  /*
    line parsing optimized for reading all number columns
    return -1 for paren error
  */
  int j, k ;
  int ic, ie ;
  int startparen, parenok, iparen ;

  /* ' " { [ (    ' " } ] ) */
  static char paren[] = { 39, 34, 123, 91, 40 } ;
  static char endparen[] = { 39, 34, 125, 93, 41 } ;
  static int nparen = 5 ;

  if( wordlist == NULL || line == NULL ) {
    PyErr_SetString(PyExc_ValueError, "lineparse called with NULL args") ;
    return -1 ;
  }

  wordlist->nwords = 0 ;
  parenok = 0 ;
  iparen = 0 ;

  /* try a quick whitespace separation */

  ic = 0 ;
  while( line[ic] != '\n' && line[ic] != '\0' ) {
    /* skip white */
    while( line[ic] == ' ' || line[ic] == '\t' ) ic++ ;
    if( line[ic] == '\n' || line[ic] == '\0' ) break ;
    /* start word */
    ie = ic + 1 ;
    while( line[ie] != ' ' && line[ie] != '\t' && 
	   line[ie] != '\n' && line[ie] != '\0' ) ie++ ;
    /* now add this word */
    if ( addWord( line, ic, ie, wordlist ) < 0 ) return -1 ;
    ic = ie ;
  }

  /* now try quick read numbers */
  readNumbers(wordlist) ;
  if( wordlist->nnum == wordlist->nwords ) return 1 ;

  /*
    parenflag:
    0 no paren/quotes checking
    1 check but dont remove blanks
    2 check and remove blanks
  */

  wordlist->nwords = 0 ;
  wordlist->nnum = 0 ;
  ic = 0 ;

  while( line[ic] != '\n' && line[ic] != '\0' ) {
    /* remove starting whitespace */
    while( line[ic] == ' ' || line[ic] == '\t') ic++ ;
    if( line[ic] == '\n' || line[ic] == '\0' ) break ;

    /* check for grouping */
    startparen = 0 ;
    if( parenflag > 0 ) {
      for( j=0 ; j<nparen ; j++ ) {
	if( line[ic] == paren[j] ) {
	  iparen = j ;
	  parenok = 0 ;
	  startparen = 1 ;
	  ic++ ;
	  break ;
	}
      }
    }

    if( startparen ) {
      /* grouped word, note we dont incl grouping symbols in parsed word */
      ie = ic ;
      while( line[ie] != '\n' && line[ie] != '\0' ) {
	if( line[ie] == endparen[iparen] ) { 
	  parenok = 1 ;
	  break ;
	}
	if( parenflag>1 && (line[ie] == ' ' || line[ie] == '\t') ) {
	  k = ie ;
	    do {
	      line[k] = line[k+1] ;
	      k++ ;
	    } while( line[k] != '\0' ) ;
	  ie-- ;
	}
	ie++ ;
      }

      if( ! parenok ) {
	wordlist->nwords-- ;
	return -1 ;
      }

    } else {
      /* regular word start */
      ie = ic + 1 ;
      while( line[ie] != ' ' && line[ie] != '\t' &&
	     line[ie] != '\n' && line[ie] != '\0' ) ie++ ;
    }

    /* now add this word */
    if ( addWord( line, ic, ie, wordlist ) < 0 ) return -1 ;
    ic = ie ;
    if( startparen ) ic++ ;
  }
  readNumbers(wordlist) ;
  return 1 ;
}

static int dblsFromTuple(PyObject *tpl, int n, double *d)
{
  PyObject *value ;
  int i, nt ;
  if( ! PyTuple_Check(tpl) ) {
    PyErr_SetString(PyExc_ValueError, "dblsFromTuple arg is NOT tuple") ;
    return -1 ;
  }
  nt = PyTuple_Size(tpl) ;
  if( nt > n ) nt = n ;
  for( i=0 ; i<nt ; i++ ) {
    if( (value = PyTuple_GetItem(tpl, i)) == NULL ) return -1 ;
    if( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError,
		      "dblsFromTuple tuple args must be numeric") ;
      return -1 ;
    }
    d[i] = PyFloat_AsDouble(value) ;
  }
  return nt ;
}
static int dblsFromList(PyObject *lst, int n, double *d)
{
  PyObject *value ;
  int i, nt ;
  if( ! PyList_Check(lst) ) {
    PyErr_SetString(PyExc_ValueError, "dblsFromList arg is NOT list") ;
    return -1 ;
  }
  nt = PyList_Size(lst) ;
  if( nt > n ) nt = n ;
  for( i=0 ; i<nt ; i++ ) {
    if( (value = PyList_GetItem(lst, i)) == NULL ) return -1 ;
    if( ! PyFloat_Check(value) && ! PyInt_Check(value) ) {
      PyErr_SetString(PyExc_ValueError,
		      "dblsFromList list args must be numeric") ;
      return -1 ;
    }
    d[i] = PyFloat_AsDouble(value) ;
  }
  return nt ;
}

static PyObject *tupleFromDbls(int n, double *d)
{
  int i ;
  PyObject *tpl, *value ;
  if( (tpl = PyTuple_New(n)) == NULL ) return NULL ;
  for( i=0 ; i<n ; i++ ) {
    if( (value = Py_BuildValue("d", d[i])) == NULL ) return NULL ;
    PyTuple_SetItem(tpl, i, value) ;
  }
  return tpl ;
}
static PyObject *listFromDbls(int n, double *d)
{
  int i ;
  PyObject *lst, *value ;
  if( (lst = PyList_New(n)) == NULL ) return NULL ;
  for( i=0 ; i<n ; i++ ) {
    if( (value = Py_BuildValue("d", d[i])) == NULL ) return NULL ;
    PyList_SetItem(lst, i, value) ;
  }
  return lst ;
}
static PyObject *atomdefTuple(char *name, ATOMdef *atomDefPtr)
{
  PyObject *atom, *Clist ;
  if( atomDefPtr->nbcoef > 0 ) {
    if( (Clist = listFromDbls(atomDefPtr->nbcoef, atomDefPtr->bcoef))
	== NULL ) return NULL ;
    if( (atom = Py_BuildValue("sidddO", name, atomDefPtr->atomic,
			      atomDefPtr->b.r, atomDefPtr->b.i,
			      atomDefPtr->m, Clist)) == NULL )
      return NULL ;
  } else if( (atom = Py_BuildValue("siddd", name, atomDefPtr->atomic,
				   atomDefPtr->b.r, atomDefPtr->b.i,
				   atomDefPtr->m)) == NULL ) {
    return NULL ;
  }
  return atom ;
}

static char *stripstring(char *s)
{
  char *ss ;
  if( s == NULL ) return s ;
  while( *s != '\0' && (*s == ' ' || *s == '\t')) s++ ;
  ss = s ;
  while( *s != '\0' ) s++ ;
  s-- ;
  while( s > ss && (*s <= 32) ) s-- ;
  s++ ;
  *s = '\0' ;
  return ss ;
}

static void clearAtomList()
{
  int i ;
  for( i=0 ; i<atomslist.n ; i++ ) atomslist.atomslist[i].omit = 1 ;
  omitAtoms() ;
}

static int checkLattice(int fix)
{
  int i, spc, ok, hexok, mono, ia, ib ;
  if( ! *usespacegroup ) return 1 ;
  ok = 1 ;
  spc = spcgrps[curSpaceGrpN].spcnum ;
  if( spc >= 195 ) {
    /* cubic a=b=c alp=bet=gam=90 */
    if( latt[1] != latt[0] ) { if(fix) latt[1] = latt[0] ; ok = 0 ; }
    if( latt[2] != latt[0] ) { if(fix) latt[2] = latt[0] ; ok = 0 ; }
    if( angl[0] != 90. ) { if(fix) angl[0] = 90. ; ok = 0 ; }
    if( angl[1] != 90. ) { if(fix) angl[1] = 90. ; ok = 0 ; }
    if( angl[2] != 90. ) { if(fix) angl[2] = 90. ; ok = 0 ; }
  } else if( spc >= 168 ) {
    /* hexagonal 2 latts= with angl betw = 60 or 120 other 2 angls 90 */
    if( latt[0] == latt[1] && latt[0] != latt[2] ) {
      if( angl[2] != 60. && angl[2] != 120. ) {
		if(fix) angl[2] = 120. ; ok = 0 ;
	  }
	  if( angl[0] != 90. ) { if(fix) angl[0] = 90. ; ok = 0 ; }
	  if( angl[1] != 90. ) { if(fix) angl[1] = 90. ; ok = 0 ; }
	} else if( latt[0] == latt[2] && latt[0] != latt[1] ) {
      if( angl[1] != 60. && angl[1] != 120. ) {
	    if(fix) angl[1] = 120. ; ok = 0 ;
      }
      if( angl[0] != 90. ) { if(fix) angl[0] = 90. ; ok = 0 ; }
      if( angl[2] != 90. ) { if(fix) angl[2] = 90. ; ok = 0 ; }
    } else if( latt[1] == latt[2] && latt[0] != latt[1] ) {
      if( angl[0] != 60. && angl[0] != 120. ) {
	    if(fix) angl[0] = 120. ; ok = 0 ;
      }
      if( angl[1] != 90. ) { if(fix) angl[1] = 90. ; ok = 0 ; }
      if( angl[2] != 90. ) { if(fix) angl[2] = 90. ; ok = 0 ; }
    } else if( latt[0] == latt[1] && latt[0] == latt[2] ) {
      /* hex with all lattice params = */
      hexok = 0 ;
      for( i=0 ; i<3 ; i++ ) {
	    if( angl[i] == 60. || angl[i] == 120. ) {
	      if( hexok ) { if(fix) angl[i] = 90. ; ok = 0 ; }
	      else { hexok = 1 ; }
	    } else if( angl[i] != 90. ) {
	      if(fix) angl[i] = 90. ; ok = 0 ;
	    }
      }
      if( ! hexok ) { if(fix) angl[2] = 120. ; ok = 0 ; }
    } else {
      /* all lattice params are different so look for hex angle */
      ok = 0 ;
      hexok = 0 ;
      for( i=0 ; i<3 ; i++ ) {
	    if( angl[i] == 60. || angl[i] == 120. ) {
	      if( hexok ) { if(fix) angl[i] = 90. ; }
	      else { hexok = i + 1 ; }
	    } else if( angl[i] != 90. ) {
	      if(fix) angl[i] = 90. ;
	    }
      }
      if( hexok ) {
	hexok-- ;
	ia = (hexok+1)%3 ;
	ib = (hexok+2)%3 ;
	if(fix) latt[ib] = latt[ia] ;
      } else {
	if(fix) angl[2] = 120. ;
	if(fix) latt[1] = latt[0] ;
      }
    }
  } else if( spc >= 143 ) {
    /* trigonal a=b=c  alp=bet=gam
    if( latt[1] != latt[0] ) { if(fix) latt[1] = latt[0] ; ok = 0 ; }
    if( latt[2] != latt[0] ) { if(fix) latt[2] = latt[0] ; ok = 0 ; }
    if( angl[1] != angl[0] ) { if(fix) angl[1] = angl[0] ; ok = 0 ; }
    if( angl[2] != angl[0] ) { if(fix) angl[2] = angl[0] ; ok = 0 ; } */
  } else if( spc >= 75 ) {
    /* tetragonal two latt == all 90 angls */
    for( i=0 ; i<3 ; i++ )
      if( angl[i] != 90. ) { if(fix) angl[i] = 90. ; ok = 0 ; }
    if( latt[1] != latt[0] && latt[2] != latt[0] && latt[1] != latt[2] ) {
      if(fix) latt[1] = latt[0] ; ok = 0 ;
    }
  } else if( spc >= 16 ) {
    /* orthorhombic all angles 90 */
    for( i=0 ; i<3 ; i++ )
      if( angl[i] != 90. ) { if(fix) angl[i] = 90. ; ok = 0 ; }
  } else if( spc >= 3 ) {
    /* monoclinic 2 angles = 90 */
    mono = 0 ;
    for( i=0 ; i<3 ; i++ ) {
      if( angl[i] != 90. ) {
	if( mono ) { if(fix) angl[i] = 90. ; ok = 0 ; }
	else { mono = 1 ; }
      }
    }
  }
  /* no conditions on triclinic spc 1 0r 2 */
  return ok ;
}

static DCMPLX bGSAS(DCMPLX b0, double E, int n, double *b)
{
  static int j[3] = { 3, 6, 8 } ;
  int i, imax ;
  double Ei, Eisq, Wsq, den ; 
  double cr[3], ci[3] ;
  DCMPLX bE ;

  if( n < 5 ) return b0 ;
  bE.r = b[0] ; /* nonresonant bReal */
  bE.i = 0. ;
  Wsq = b[4]*b[4] ;
  cr[0] = b[1] ;
  ci[0] = -b[2] ;
  if( n > 5 ) { cr[1] = b[1]*b[5] ; ci[1] = -b[2]*b[5] ; }
  if( n > 7 ) { cr[2] = b[1]*b[7] ; ci[2] = -b[2]*b[7] ; }
  imax = (n-1)/2 - 1 ;
  for( i=0 ; i < (n <= 9 ? imax : 3) ; i++ ) {
    Ei = E - b[j[i]] ;
    Eisq = Ei*Ei ;
    den = Eisq + Wsq ;
    /*
      I think the formula in GSAS manual is wrong
      and Wsq is added not subtracted in the denominator
    */
    bE.r += cr[i]*Ei/den ;
    bE.i += ci[i]/den ;
  }
  return bE ;
}
static DCMPLX bresonant(double k, DCMPLX b0, double Eres, double HW)
{
  /*
    single resonance fit for b wavelength dependence as in M&L.
    N.B. b = -f(scattering amplitude)
    Given b at 1.8 Angstroms = 25.25 meV  k = 3.4906585  and
    resonance energy Eres and HalfWidth HW
    return b at input energy E
    using resonant real and imag parts of b as M&Lovesey Appendix A,
    br = (HWn/k) (E - Eres)/den
    bi = -(HWn/k)(HW + HWn)/den
    den = (E - Eres)^2 + (HW + HWn)^2
    HWn/HW @ Eres ~= 0.006  HWn ~ k  E = Dn k^2 OR k = sqrt(E/Dn)
    Dn = 2.072141789 mev A^2 = 2.072141789 * (10^-8)^2 mev cm^2
    so k[cm-1] = sqrt(E[meV]/2.072141789)*10^12*10^-4
    and 1/k [cm] = 10^4*sqrt(2.07214/E)*10^-12 [cm]
    Let HWn = HWn0 (k/k0). Since we are going to fudge HWn0 anyway
    we can arbitrarily pick k0
    Now
    br = (HWn0/k0) (E - Eres)/den
    bi = (HWn0/k0) (HW + HWn0(k/k0))/den
    den = (E - Eres)^2 + (HW + HWn0(k/k0))^2
    Use 1/k0 [cm] = 10^4 * sqrt(2.07214/En0) * 10^-12 [cm]
    factor Eres^2 out of den
    br =  (HWn0/Eres) 10^4 sqrt(2.07214/En0) 10^-12 [cm] (E/Eres - 1)/den'
    bi = -(HWn0/Eres) 10^4 sqrt(2.07214/En0) 10^-12 [cm] (HW + HWn)/Eres/den'
    den' = (E/Eres - 1)^2 + ((HW + HWn)/Eres)^2
    So for simplicity take k0 = 1 A^-1 so En0 = 2.07214 and then
    br =  (HWn1/Eres) 10^4 10^-12 [cm] (E/Eres - 1)/den'
    bi =  (HWn1/Eres) 10^4 10^-12 [cm] (HW + HWn1*k)/Eres/den'
    den' = (E/Eres - 1)^2 + ((HW + HWn1*k)/Eres)^2
    given Eres and HW
    fudge HWn1 to get the imaginary part to agree at 25.25 meV
    At 25.25 mev lambda = 1.8 Ang  k = 2Pi/lambda = 3.4906585
    solving
    HW1 = -B + B sqrt(1 + bi [(25.25-Eres)^2 + HW^2]/3.49/B^2/(10^4 - e)
    where B = HW (10^4 - 2e)/(10^4 - e)/(2*3.49)
    and e = 3.49 * bi  bi in units 10^-12 cm
    or
    HW1 = -B + B sqrt(1 + bi (4*3.49/10^4)((1-10^-4e)/(1-2*10^-4e)^2)/L)
    L = HW^2/((25.25 - Er)^2 + HW^2)
   */
  DCMPLX b ;
  double E, DE, H, HWn1, HWn, Ht, DE2, L, e, f, bi, e1, e2, r, r2, den ;
  double b0res, b0non ;
  static double E0 = 25.24908 ;
  static double k0 = 3.4906585 ;
  static double Dn = 2.072141789 ;

  if( Eres <= 0. ) return b0 ;
  bi = b0.i ;
  e = bi*k0 ;
  e1 = 1. - e/10000. ;
  e2 = 1. - 2.*e/10000. ;
  H = HW*e2/e1/(2.*k0) ;
  r = HW/Eres ;
  r2 = r*r ;
  DE = E0/Eres - 1. ;
  L = 10000.*r2/(r2 + DE*DE) ;
  f = 4.*k0*(e1/e2/e2)/L ;
  HWn1 = -H + H*sqrt(1 + bi*f) ;
  /*
    HWn1 is supposed to be HWn at k=1 E=Dn very low neutron energy 
    Now we can calculate the real and imaginary parts due to resonance at E 
    b0.r = non-res + res at k=1.8 A 
    We need to calc the non-res real part 
    by calculating the resonant real part at 1.8
  */
  f = 10000.*HWn1/Eres ;

  HWn = HWn1 * k0 ;
  Ht = (HW + HWn)/Eres ;
  DE = E0/Eres - 1. ;
  den = DE*DE + Ht*Ht ;
  b0res = f + DE/den ;
  b0non = b0.r - b0res ;

  E = Dn*k*k ;
  HWn = HWn1 * k ;
  Ht = (HW + HWn)/Eres ;
  DE = E/Eres - 1. ;
  den = DE*DE + Ht*Ht ;

  b.i = -f * Ht/den ;
  b.r =  b0non + f * DE/den ;
  return b ;
}

static PyObject *MakeAtomList(char *name, int isub, myATOM *a)
{
  PyObject *groupList ;
  if ( (groupList = PyList_New(10)) == NULL ) return NULL ;
  PyList_SET_ITEM(groupList, 0,Py_BuildValue("s", name));
  PyList_SET_ITEM(groupList, 1, Py_BuildValue("i", isub)) ;
  PyList_SET_ITEM(groupList, 2, Py_BuildValue("d", a->r[0])) ;
  PyList_SET_ITEM(groupList, 3, Py_BuildValue("d", a->r[1])) ;
  PyList_SET_ITEM(groupList, 4, Py_BuildValue("d", a->r[2])) ;
  PyList_SET_ITEM(groupList, 5, Py_BuildValue("d", a->s[0])) ;
  PyList_SET_ITEM(groupList, 6, Py_BuildValue("d", a->s[1])) ;
  PyList_SET_ITEM(groupList, 7, Py_BuildValue("d", a->s[2])) ;
  PyList_SET_ITEM(groupList, 8, Py_BuildValue("d", a->m)) ;
  PyList_SET_ITEM(groupList, 9, Py_BuildValue("d", a->c)) ;
  return groupList ;
}
