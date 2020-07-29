#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _hhaxon_reg();
extern void _hhsoma_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," hhaxon.mod");
fprintf(stderr," hhsoma.mod");
fprintf(stderr, "\n");
    }
_hhaxon_reg();
_hhsoma_reg();
}
