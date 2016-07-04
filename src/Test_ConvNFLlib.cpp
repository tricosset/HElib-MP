
#define __TEST_CONVNFLLIB__

#include <NTL/ZZ.h>
#include "FHEContext.h"
#include "FHE.h"
#include "Ctxt.h"
#include <time.h>
#include <ctime>
#include <fstream>
#include "Util.h"
#include "DoubleCRT.h"
#include "NTL/ZZ_pX.h"
#include <chrono>
#include <valgrind/callgrind.h>
#include <sys/time.h>
#include "Test_Params.hpp" // Parameters

int main() {
	SetSeed(ZZ(0));

	struct timeval tbeg, tend;
	double texe = 0;


	/*
	 * Instantiation
	 */
	cout << endl
		 << "***************************" << endl
		 << "* Test ConvFrom/To NFLlib *" << endl
#ifdef __TEST_ECDSA__
		 << "*        ECDSA 256        *" << endl
#endif
#ifdef __TEST_RSA__
		 << "*        RSA 2048         *" << endl
#endif
		 << "***************************" << endl;
	cout << "   Parameters" << endl
		 << "---------------------------" << endl
	     << "  p:         " << plaintextModulus << endl
	     << "  m:         " << m 				<< endl
	     << "  wndw:      " << WNDW 			<< endl
	     << "  hght:      " << hght 			<< endl
	     << "  nPrms:     " << nPrms 			<< endl
	     << "  nDgts:     " << nDgts 			<< endl;
	FHEcontext context(m, plaintextModulus, g_GFm);
	buildModChain(context, lvl, nDgts, nHlfPrmsByLvl);


	/*
	 * Keys Generation
	 */
	cout << "===========================" << endl
		 << "   KeyGen"                   << endl;

	FHESecKey secretKey(context);
	const FHEPubKey& publicKey = secretKey;
	secretKey.GenSecKey(512, plaintextModulus); // A Hamming-weight-w secret key


	/*
	 * String to DoubleCRT
	 */
	cout << "===========================" << endl
		 << "   " << nb_encrypt << " String to DoubleCRT" << endl
		 << "---------------------------" << endl;
	Ctxt* c[25];
	gettimeofday(&tbeg,NULL);
	for (unsigned i = 0; i < 25; i++) {
		c[i] = new Ctxt(publicKey, "c0_out", "c1_out", plaintextModulus);
	}
	gettimeofday(&tend,NULL);
	texe = ((double)(tend.tv_sec-tbeg.tv_sec)) + ((double)(tend.tv_usec-tbeg.tv_usec))/1000000.;
	cout << "  log2 q:    " << context.logOfProduct(c[0]->getPrimeSet())/log(2) << endl
		 << "  Time:      " << texe/25.*((double)nb_encrypt) << " s" << std::endl
	     << "  Size:      " << context.logOfProduct(c[0]->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM()*((double)nb_encrypt) << " Mb" << std::endl
	     << "  Flow rate: " << (context.logOfProduct(c[0]->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM())/(texe/25.) << " Mbps" << std::endl;


	/*
	 * Reduction to last q
	 */
	for (unsigned i = 0; i < 25; i++) {
		c[i]->modDownToLevel(lastBaseLvl);
	}


	/*
	 * DoubleCRT to String
	 */
	cout << "===========================" << endl
		 << "   " << nb_decrypt*2 << " DoubleCRT to String" << endl
		 << "---------------------------" << endl;
	ofstream file0("newc0_out", ios::out | ios::trunc);
	ofstream file1("newc1_out", ios::out | ios::trunc);
	ZZX tmp0, tmp1;
	gettimeofday(&tbeg,NULL);
	for (unsigned i = 0; i < 25; i++) {
		c[i]->parts[0].toPoly(tmp0, true);
		c[i]->parts[1].toPoly(tmp1, true);
	}
	gettimeofday(&tend,NULL);
	texe = ((double)(tend.tv_sec-tbeg.tv_sec)) + ((double)(tend.tv_usec-tbeg.tv_usec))/1000000.;
	cout << "  log2 q:    " << context.logOfProduct(c[0]->getPrimeSet())/log(2) << endl
		 << "  Time:      " << texe/25.*((double)nb_decrypt) << " s" << endl
	     << "  Size:      " << context.logOfProduct(c[0]->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM()*2.*((double)nb_decrypt) << " Mb" << std::endl
	     << "  Flow rate: " << (context.logOfProduct(c[0]->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM())*2./(texe/25.) << " Mbps" << std::endl;
	cout << "***************************" << endl;

	file0 << tmp0;
	file1 << tmp1;
	file0.close();
	file1.close();
	return 0;
}
