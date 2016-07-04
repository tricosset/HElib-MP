
#define __TEST_RSA_2048__

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
		 << "*      Test RSA 2048      *" << endl
		 << "***************************" << endl;
	cout << "   Parameters"               << endl
	     << "---------------------------" << endl
	     << "  p:             " << plaintextModulus << endl
	     << "  m:             " << m 				  << endl
	     << "  depth:         " << lvl 				  << endl
	     << "  wndw:          " << WNDW 			  << endl
	     << "  hght:          " << hght 			  << endl
	     << "  nPrms:         " << nPrms 			  << endl
	     << "  nDgts:         " << nDgts 			  << endl;
	FHEcontext context(m, plaintextModulus);
	buildModChain(context, lvl, nDgts, nHlfPrmsByLvl);
   cout << "  phi(m):        " << context.zMStar.getPhiM() << endl;

	/*
	 * Keys Generation
	 */
	cout << "===========================" << endl
	     << "   KeyGen"                   << endl;
	FHESecKey secretKey(context);
	const FHEPubKey& publicKey = secretKey;
	secretKey.GenSecKey(512, plaintextModulus); // A Hamming-weight-w secret key
	ZZX p = to_ZZX(plaintextModulus-1);
	Ctxt *c[2048/WNDW];
	for (unsigned i = 0; i < 2048/WNDW; i++)
	{
		c[i] = new Ctxt(publicKey);
	}


	/*
	 * Encryptions
	 */
	cout << "==========================="    << endl
	     << "   " << 2048/WNDW << " Encrypt" << endl
	     << "---------------------------"    << endl;
	gettimeofday(&tbeg,NULL);
	for(long i = 0; i < 2048/WNDW; i++){
		secretKey.Encrypt(*c[i], p, plaintextModulus);
		Debug(cout << "." << flush);
	}
	Debug(cout << endl);
	gettimeofday(&tend,NULL);
	texe = ((double)(tend.tv_sec-tbeg.tv_sec)) + ((double)(tend.tv_usec-tbeg.tv_usec))/1000000.;
	cout << "  log2 q:      " << context.logOfProduct(c[0]->getPrimeSet())/log(2) << endl
		 << "  Time:        " << texe << " s" << endl
		 << "  Avg:         " << texe/((double)(2048/WNDW)) << " s" << endl
	     << "  Size:        " << context.logOfProduct(c[0]->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM()*(double)(2048/WNDW) << " Mb" << std::endl
	     << "  Rate:        " << (context.logOfProduct(c[0]->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM())/(texe/((double)(2048/WNDW))) << " Mbps" << std::endl;


	/*
	 * Multiplications
	 */
	cout << "===========================" << endl
	     << "   ";
	     for (unsigned i = 0; i < hght-1; i++) {
	    	 cout << ((2048/(2*WNDW))>>i) << "+";
	     }
	cout << ((2048/(2*WNDW))>>(hght-1)) << " Mul." << endl
	     << "---------------------------"          << endl;
	gettimeofday(&tbeg,NULL);
	long k;
	for (k=2; k<(1<<(hght+1)); k<<=1) {
		for (unsigned i=0; i<2048/WNDW; i+=k) {
			c[i]->multiplyBy(*c[i+k/2]);
			Debug(cout << "." << flush);
			Debug(for(unsigned j=0; j<k-1; j++)cout << " " << flush);
		}
		Debug(cout << endl);
	}
	gettimeofday(&tend,NULL);
	texe = ((double)(tend.tv_sec-tbeg.tv_sec)) + ((double)(tend.tv_usec-tbeg.tv_usec))/1000000.;
	cout << "  Time:        " << texe << " s" << std::endl
	     << "  Avg:         " << texe/(double)((2048/WNDW)-1) << " s" << endl;

	/*
	 * Squashing
	 */
	for(unsigned i = 0; i < 2048/WNDW; i+=k/2){
		c[i]->modDownToLevel(c[i]->findBaseLevel());
		c[i]->cleanUp();
	}


	/*
	 * Decryptions
	 */
	cout << "==========================="          << endl
	     << "   " << 2048/WNDW/(k/2) << " Decrypt" << endl
	     << "---------------------------"          << endl;
	gettimeofday(&tbeg,NULL);
	for(unsigned i = 0; i < 2048/WNDW; i+=k/2) {
		secretKey.Decrypt(p, *c[i]);
		Debug(cout << "." << flush);
		Debug(for(unsigned j=0; j<k/2-1; j++)cout << " " << flush);
	}
	Debug(cout << endl);
	gettimeofday(&tend,NULL);
	texe = ((double)(tend.tv_sec-tbeg.tv_sec)) + ((double)(tend.tv_usec-tbeg.tv_usec))/1000000.;
	cout << "  log2 q:      " << context.logOfProduct(c[0]->getPrimeSet())/log(2) << endl
		 << "  Correctness: " << ((p[0]==to_ZZ(1))?"true":"false") << endl
         << "  Time:        " << texe << " s" << std::endl
	     << "  Avg:         " << texe/(double)(2048/WNDW/(k/2)) << " s" << endl
	     << "  Size:        " << context.logOfProduct(c[0]->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM()*((double)((2048/WNDW)/(k/2)))*2. << " Mb" << std::endl
	     << "  Rate:        " << (context.logOfProduct(c[0]->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM())*((double)((2048/WNDW)/(k/2)))*2./texe << " Mbps" << std::endl
         << "===========================" << endl;


	return 0;
}
