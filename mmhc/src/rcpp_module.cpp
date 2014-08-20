#include "mmhc.h"

RCPP_MODULE(mmhc) {
    class_<MMHC>("MMHC")
        .default_constructor()
		.constructor<SEXP>()
		.method("pc", &MMHC::GetPC, "returns the PC set from MMPC()")
		.method("adjMat", &MMHC::GetGraph, "returns the adjacency from MMHC()")
		.method("score", &MMHC::GetScore, "return the score of the final graph")
		.method("mmpc", &MMHC::mmpc, "executes the MMPC() function")
		.method("mmhc", &MMHC::mmhc, "executes the MMHC() function")
	;
}
