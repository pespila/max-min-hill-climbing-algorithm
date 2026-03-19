#include "mmhc.h"

// The RCPP_MODULE sets which elements and methods (including the constructor) of the class can be accessed
RCPP_MODULE(mmhc) {
    class_<MMHC>("MMHC")
        .default_constructor()
		.constructor<SEXP>()
        .method("mat", &MMHC::GetMat, "returns the matrix from MMPC()")
		.method("pc", &MMHC::GetPC, "returns the PC set from MMPC()")
		.method("adjMat", &MMHC::GetGraph, "returns the adjacency from MMHC()")
		.method("score", &MMHC::GetScore, "return the score of the final graph")
		.method("mmpc", &MMHC::mmpc, "executes the MMPC() function")
		.method("mmhc", &MMHC::mmhc, "executes the MMHC() function")
	;
}
