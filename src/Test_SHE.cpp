
#define __TEST_SHE_512__

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
#include <sys/time.h>
#include "Test_Params.hpp"


int main() {
	SetSeed(ZZ(0));
	struct timeval tbeg, tend;
	double texe = 0;

	/*
	 * Instantiation
	 */
	cout << endl
		 << "***************************" << endl
		 << "*        Test SWHE        *" << endl
		 << "***************************" << endl;
	cout << "   Parameters"               << endl
	     << "---------------------------" << endl
	     << "  p:           " << plaintextModulus << endl
	     << "  m:           " << m 				  << endl
	     << "  depth:       " << lvl 			  << endl
	     << "  nPrms:       " << nPrms 			  << endl
	     << "  nDgts:       " << nDgts 			  << endl;
	FHEcontext context(m, plaintextModulus);
	buildModChain(context, lvl, nDgts, nHlfPrmsByLvl);

	/*
	 * Keys Generation
	 */
	cout << "===========================" << endl
	     << "   KeyGen"                   << endl;
	FHESecKey secretKey(context);
	const FHEPubKey& publicKey = secretKey;
	secretKey.GenSecKey(32, plaintextModulus); // A Hamming-weight-w secret key
	ZZX p = to_ZZX(plaintextModulus-1);
	Ctxt *c = new Ctxt(publicKey);

	/*
	 * Encryptions
	 */
	cout << "===========================" << endl
	     << "   Encrypt"                  << endl
	     << "---------------------------" << endl;
	gettimeofday(&tbeg,NULL);
	for(long i = 0; i < repeti; i++){
		secretKey.Encrypt(*c, p, plaintextModulus);
		Debug(cout << "." << flush);
	}
	Debug(cout << endl);
	gettimeofday(&tend,NULL);
	texe = ((double)(tend.tv_sec-tbeg.tv_sec)) + ((double)(tend.tv_usec-tbeg.tv_usec))/1000000.;
	cout << "  log2 q:      " << context.logOfProduct(c->getPrimeSet())/log(2) << endl
		 << "  Time:        " << texe/(double)repeti << " s" << endl
	     << "  Size:        " << context.logOfProduct(c->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM()*2 << " Mb" << std::endl
	     << "  Rate:        " << (context.logOfProduct(c->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM())*2/(texe/((double)(repeti))) << " Mbps" << std::endl;


	/*
	 * Multiplications
	 */
	cout << "===========================" << endl
	     << "   " << lvl << " Mul. (1 per lvl)"       << endl
	     << "---------------------------" << endl;
	gettimeofday(&tbeg,NULL);
	for (int i=0; i<lvl; i++) {
		c->multiplyBy(*c);
		Debug(cout << "." << flush);
	}
	Debug(cout << endl);
	gettimeofday(&tend,NULL);
	texe = ((double)(tend.tv_sec-tbeg.tv_sec)) + ((double)(tend.tv_usec-tbeg.tv_usec))/1000000.;
	cout << "  Time:        " << texe << " s" << std::endl;
	cout << "  Avg:         " << texe/lvl << " s" << std::endl;


	/*
	 * Squashing
	 */
	c->modDownToLevel(c->findBaseLevel());
	c->cleanUp();


	/*
	 * Decryptions
	 */
	cout << "===========================" << endl
	     << "   Decrypt"                  << endl
	     << "---------------------------" << endl;
	gettimeofday(&tbeg,NULL);
	for(unsigned i = 0; i < repeti; i++) {
		secretKey.Decrypt(p, *c);
		Debug(cout << "." << flush);
	}
	Debug(cout << endl);
	gettimeofday(&tend,NULL);
	texe = ((double)(tend.tv_sec-tbeg.tv_sec)) + ((double)(tend.tv_usec-tbeg.tv_usec))/1000000.;
	cout << "  log2 q:      " << context.logOfProduct(c->getPrimeSet())/log(2) << endl
		 << "  Correctness: " << ((p[0]==to_ZZ(1))?"true":"false") << endl
         << "  Time:        " << texe/(double)repeti << " s" << std::endl
	     << "  Size:        " << context.logOfProduct(c->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM()*2. << " Mb" << std::endl
	     << "  Rate:        " << (context.logOfProduct(c->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM())*((double)repeti)*2./texe << " Mbps" << std::endl
         << "===========================" << endl;
}
