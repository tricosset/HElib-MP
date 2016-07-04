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
#include "elliptic_curve.hpp"
#include "Test_Params.hpp"

int main() {

	  SetSeed(ZZ(0));
	  long n = 1;

	  struct timeval tbeg, tend;
	  double texe = 0;

	  FHEcontext context(m, plaintextModulus, g_GFm);
	  buildModChain(context, 0, nDgts, 296);

		cout << "===========================" << endl
		     << "   Test ECDSA "              << endl
		     << "===========================" << endl;
		cout << "Parameters: " << endl
		     << "\tm: " << m << endl
		     << "\tlog2 q: " << context.logOfProduct(context.ctxtPrimes)/log(2) << endl
		     << "\tp: " << plaintextModulus << endl
		     << "\tnDgts: " << nDgts << endl;

		cout << "===========================" << endl
		     << "   KeyGen"                   << endl
		     << "===========================" << endl;
	  FHESecKey secretKey(context);
	  const FHEPubKey& publicKey = secretKey;
	  secretKey.GenSecKey(32, plaintextModulus); // A Hamming-weight-w secret key

      ECPrecomputation precomputation_encrypted;
		Ctxt cG_x(publicKey);
		Ctxt cG_y(publicKey);
		Ctxt cG_z(publicKey);
		Ctxt cG_t(publicKey);

		ZZX pG_x = to_ZZX(G_x);
		ZZX pG_y = to_ZZX(G_y);
		ZZX pG_z = to_ZZX(G_z);
		ZZX pG_t = to_ZZX(G_t);

		gettimeofday(&tbeg,NULL);
		for(long i = 0; i < n; i++){
			  CALLGRIND_START_INSTRUMENTATION;
			secretKey.Encrypt(cG_x, pG_x, plaintextModulus);
			  CALLGRIND_STOP_INSTRUMENTATION;
			  CALLGRIND_DUMP_STATS;
			secretKey.Encrypt(cG_y, pG_y, plaintextModulus);
			secretKey.Encrypt(cG_z, pG_z, plaintextModulus);
			secretKey.Encrypt(cG_t, pG_t, plaintextModulus);
		}
		gettimeofday(&tend,NULL);
		texe = ((double)(tend.tv_sec-tbeg.tv_sec)) + ((double)(tend.tv_usec-tbeg.tv_usec))/1000000.;
		cout << "===========================" << endl;
		std::cout  << "Chiffrement avec la cle secrete:" << std::endl
				   << "\tTemps: " << texe/(4.*(double)n) << " s" << std::endl
				   << "\tTaille: " << context.logOfProduct(cG_x.getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM() << " Mb" << std::endl
				   << "\tDebit: " << (context.logOfProduct(cG_x.getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM())/(texe/(4.*(double)n)) << " Mbps" << std::endl;
		cout << "===========================" << endl
		     << "   Homom. 16G"               << endl
		     << "===========================" << endl;
		ECPoint eG(cG_x, cG_y, cG_z, cG_t, publicKey);
		ECPoint eGpG(publicKey);

	  gettimeofday(&tbeg,NULL);

	  ec_addition(eGpG,   eG,   eG, precomputation_encrypted, publicKey); //  2G
	  ec_addition(eGpG, eGpG, eGpG, precomputation_encrypted, publicKey); //  4G
	  ec_addition(eGpG, eGpG, eGpG, precomputation_encrypted, publicKey); //  8G
	  ec_addition(eGpG, eGpG, eGpG, precomputation_encrypted, publicKey); // 16G

	  gettimeofday(&tend,NULL);
	  texe = ((double)(tend.tv_sec-tbeg.tv_sec)) + ((double)(tend.tv_usec-tbeg.tv_usec))/1000000.;
	  std::cout << texe << " s" << std::endl;

	  gettimeofday(&tbeg,NULL);
	  for(long i = 0; i < n; i++){
			secretKey.Decrypt(pG_x, *eGpG.X);
			secretKey.Decrypt(pG_y, *eGpG.Y);
			secretKey.Decrypt(pG_z, *eGpG.Z);
			secretKey.Decrypt(pG_t, *eGpG.T);
	  }
		gettimeofday(&tend,NULL);

		ZZ mG_x = pG_x[0];
		ZZ mG_y = pG_y[0];
		ZZ mG_z = pG_z[0];
		ZZ mG_t = pG_t[0];

		ZZ invZ = InvMod(mG_z, plaintextModulus);
		mG_x *= invZ;
		mG_x %= plaintextModulus;
		mG_y *= invZ;
		mG_y %= plaintextModulus;

		std::cout << "\t"  << mG_x << std::endl;
		std::cout << "\t"  << mG_y << std::endl;

		texe = ((double)(tend.tv_sec-tbeg.tv_sec)) + ((double)(tend.tv_usec-tbeg.tv_usec))/1000000.;
		cout << "===========================" << endl;
		std::cout  << "Dechiffrement:" << std::endl
				   << "\tTemps: " << texe/(4.*(double)n) << " s" << std::endl
				   << "\tTaille: " << context.logOfProduct(eGpG.X->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM() << " Mb" << std::endl
				   << "\tDebit: " << (context.logOfProduct(eGpG.X->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM())/(texe/(4.*(double)n)) << " Mbps" << std::endl;
		cout << "===========================" << endl;

		std::cout << "eGpG.X useless primes: " << eGpG.X->uselessPrimes(secretKey, pG_x) << std::endl;
		std::cout << "eGpG.Y useless primes: " << eGpG.Y->uselessPrimes(secretKey, pG_y) << std::endl;
		std::cout << "eGpG.Z useless primes: " << eGpG.Z->uselessPrimes(secretKey, pG_z) << std::endl;
		std::cout << "eGpG.T useless primes: " << eGpG.T->uselessPrimes(secretKey, pG_t) << std::endl;

		cout << "===========================" << endl;
  return 0;
}
