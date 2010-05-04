/*
  declarations for use of libPB

  in PBindata and PBoutdata the cross-sections flipper states are
  pp = ++ = PoffAoff
  mm = -- = PonAon
  pm = +- = PonAoff
  mp = -+ = PoffAon
  mm = -- = PoffAoff

  in PBindata:
  Ei Ef are incident and exit energies in meV for each datapoint of up to 4 CS msrd
  Cpp Cmm Cpm Cmp are the 4 msrd count-rates for each datapoint.
  Epp Emm Epm Emp are the stnd dev of the msrd count-rates
  tpp ... are the time stamps in seconds

  in PBoutdata:
  Spp Smm ... are the extracted cross-sections
  Epp Emm ... are the stnd dev of the extracted cross-sections

  See the example .CTRL file to help understand what to put into
  the flags for CountsEnable CountsAdd Sconstrain Spp ...
*/

typedef struct {
  double *Ei, *Ef ;
  double *Cpp, *Cmm, *Cpm, *Cmp ;
  double *Epp, *Emm, *Epm, *Emp ;
  unsigned long *tpp, *tmm, *tpm, *tmp ;
} PBindata ;

typedef struct {
  double *Spp, *Smm, *Spm, *Smp ;
  double *Epp, *Emm, *Epm, *Emp ;
  double *R ;
} PBoutdata ;

typedef struct {
  int MonitorCorrect ;
  int PolMonitorCorrect ;
  int MonoSelect ;  /* 0 is unfiltered PG, 1 is 2cmPGfilter E=14.7 */
  int Debug ;
  int SimFlux ;
  int SimDeviate ;
  int NoNegativeCS ;
  int HalfPolarized ;

  int CountsEnable[4] ;
  int CountsAdd1[4] ;
  int CountsAdd2[4] ;
  int Sconstrain[4] ;
  double Spp[4], Smm[4], Spm[4], Smp[4] ;
} PBflags ;


/* entrypoints */

int PBcorrectData(char *CellFile, PBflags *flgs,
		  int npts, PBindata *in, PBoutdata *out) ;

int PBsim(char *filename) ;

int PBreadflags(char *filename) ;

int PBreaddata(char *filename) ;
