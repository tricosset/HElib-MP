
//#define __TEST_ECC_P256__
#define __AUTO_TEST_ECC_P256__

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
#include "elliptic_curve.hpp"

int main() {
#ifdef __AUTO_TEST_ECC_P256__
	for (unsigned wndw = 8; wndw < 256; wndw*=2) {
		unsigned hght_max = 0;
		while ((256/wndw)>>(++hght_max) > 0);
		for (unsigned hght = 1; hght < hght_max; hght++) {
			unsigned lvl = hght;
			unsigned nPrms = (lvl+1)*nHlfPrmsByLvl/2;
#endif
			SetSeed(ZZ(0));
			struct timeval tbeg, tend;
			double texe = 0;

			ZZX pG_x = to_ZZX(G_x);
			ZZX pG_y = to_ZZX(G_y);
			ZZX pG_z = to_ZZX(G_z);

			Ctxt* cG_x[256];
			Ctxt* cG_y[256];
			Ctxt* cG_z[256];

			ECPoint* eG[256];

			/*
			 *  Instantiation
			 */
			cout << endl
				 << "***************************" << endl
				 << "*     Test ECDSA 256      *" << endl
				 << "***************************" << endl;
			cout << "   Parameters"               << endl
			     << "---------------------------" << endl
			     << "  p:         " << plaintextModulus << endl
			     << "  m:         " << m 				<< endl
			     << "  depth:     " << lvl		     	<< endl
			     << "  wndw:      " << wndw 			<< endl
			     << "  hght:      " << hght 			<< endl
			     << "  nPrms:     " << nPrms 			<< endl
			     << "  nDgts:     " << nDgts 			<< endl;

			FHEcontext context(m, plaintextModulus);
			/*
			 * Ifndef NO_HALF_SIZE_PRIME then nPrimes = (nLevels+1)*nPrimesByLvl/2
			 */
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
			 * Encryption
			 */
			for(long i = 0; i < 256/wndw; i++) {
				cG_x[i] = new Ctxt(publicKey);
				cG_y[i] = new Ctxt(publicKey);
				cG_z[i] = new Ctxt(publicKey);
			}

			cout << "===========================" << endl
			     << "   " << 256/wndw << "*(x,y) Encrypt" << endl
			     << "---------------------------" << endl;
			gettimeofday(&tbeg,NULL);
			for(long i = 0; i < 256/wndw; i++) {
				secretKey.Encrypt(*cG_x[i], pG_x, plaintextModulus);
				secretKey.Encrypt(*cG_y[i], pG_y, plaintextModulus);
				Debug(cout << "." << flush);
			}
			gettimeofday(&tend,NULL);
			Debug(cout << endl);
			texe = ((double)(tend.tv_sec-tbeg.tv_sec)) + ((double)(tend.tv_usec-tbeg.tv_usec))/1000000.;
			cout << "  log2 q:    " << context.logOfProduct(cG_x[0]->getPrimeSet())/log(2) << endl
				 << "  Time:      " << texe << " s" << endl
				 << "  Avg:       " << texe/(256/wndw) << " s" << endl
			     << "  Size:      " << context.logOfProduct(cG_x[0]->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM()*((double)(256/wndw)*2) << " Mb" << endl
			     << "  Rate:      " << (context.logOfProduct(cG_x[0]->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM())*((double)(256/wndw)*2)/(texe) << " Mbps" << endl;


			/*
			 * Z-coordinates Encryption
			 */
			for(long i = 0; i < 256/wndw; i++) {
				secretKey.Encrypt(*cG_z[i], pG_z, plaintextModulus);
			}


			/*
			 * Precomputation
			 */
			ECPrecomputationEncrypted precomputation_encrypted(context);
			for(long i = 0; i < 256/wndw; i++) {
				eG[i] = new ECPoint(*cG_x[i], *cG_y[i], *cG_z[i], publicKey);
			}


			/*
			 * Scalar Multiplication
			 */
			cout << "===========================" << endl
			     << "   ";
			     for (unsigned i = 0; i < hght-1; i++) {
			    	 cout << ((256/(2*wndw))>>i) << "+";
			     }
			cout << ((256/(2*wndw))>>(hght-1)) << " Add." << endl
			     << "---------------------------" << endl;
			gettimeofday(&tbeg,NULL);
			long k;
			for (k=2; k<1<<(hght+1); k<<=1) {
				for (long i=0; i<256/wndw; i+=k) {
					ec_addition(*eG[i], *eG[i], *eG[i+k/2], precomputation_encrypted, publicKey);
					Debug(cout << "." << flush);
					Debug(for(unsigned j=0; j<k-1; j++)cout << " " << flush);
				}
				Debug(cout << endl);
			}
			gettimeofday(&tend,NULL);
			texe = ((double)(tend.tv_sec-tbeg.tv_sec)) + ((double)(tend.tv_usec-tbeg.tv_usec))/1000000.;
			cout << "  Time:      " << texe << " s" << std::endl
			     << "  Avg:       " << texe/(double)((256/wndw)-1) << " s" << endl;


			/*
			 * Squashing
			 */
			for(long i = 0; i < 256/wndw; i+=k/2) {
				eG[i]->X->modDownToLevel(eG[i]->X->findBaseLevel());
				eG[i]->Y->modDownToLevel(eG[i]->Y->findBaseLevel());
				eG[i]->Z->modDownToLevel(eG[i]->Z->findBaseLevel());

				eG[i]->X->cleanUp();
				eG[i]->Y->cleanUp();
				eG[i]->Z->cleanUp();
			}


			/*
			 * Decryption
			 */
			cout << "===========================" << endl
			     << "   " << 256/wndw/(k/2) << "*(x,y,z) Decrypt" << endl
			     << "---------------------------" << endl;
			gettimeofday(&tbeg,NULL);
			for(long i = 0; i < 256/wndw; i+=k/2) {
				secretKey.Decrypt(pG_x, *eG[i]->X);
				secretKey.Decrypt(pG_y, *eG[i]->Y);
				secretKey.Decrypt(pG_z, *eG[i]->Z);
				Debug(cout << "." << flush);
				Debug(for(unsigned j=0; j<k/2-1; j++)cout << " " << flush);
			}
			gettimeofday(&tend,NULL);
			Debug(cout << endl);
			texe = ((double)(tend.tv_sec-tbeg.tv_sec)) + ((double)(tend.tv_usec-tbeg.tv_usec))/1000000.;

			ZZ mG_x = pG_x[0];
			ZZ mG_y = pG_y[0];
			ZZ mG_z = pG_z[0];
			ZZ invZ = InvMod(mG_z, plaintextModulus);

			mG_x *= invZ;
			mG_x %= plaintextModulus;

			mG_y *= invZ;
			mG_y %= plaintextModulus;

			cout << "  log2 q:    " << context.logOfProduct(eG[0]->X->getPrimeSet())/log(2) << endl
				 << "  " << (1<<hght) << "G_x:     " << mG_x << endl
				 << "  " << (1<<hght) << "G_y:     " << mG_y << endl
		         << "  Time:      " << texe << " s" << endl
				 << "  Avg:       " << texe/(256/wndw/(k/2)) << " s" << endl
			     << "  Size:      " << context.logOfProduct(eG[0]->X->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM()*((double)(256/wndw/(k/2)))*6. << " Mb" << endl
			     << "  Rate:      " << (context.logOfProduct(eG[0]->X->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM())*((double)(256/wndw/(k/2)))*6./texe << " Mbps" << endl
			     << "===========================" << endl;
#ifdef __AUTO_TEST_ECC_P256__
		}
	}
#endif
	return 0;
}
