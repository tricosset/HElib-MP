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

#define ECDSA_256
//#define RSA_2048

#ifdef ECDSA_256

unsigned m = 1031;
unsigned g_GFm = 14;
unsigned nDgts = 3;
unsigned nPrimesByLvl = 29;
unsigned lvl = 9;

// Curve P256
ZZ plaintextModulus = conv<ZZ>("115792089210356248762697446949407573530086143415290314195533631308867097853951");
ZZ coeff_a = plaintextModulus-3;
ZZ coeff_b = conv<ZZ>("41058363725152142129326129780047268409114441015993725554835256314039467401291");
ZZ order_G = conv<ZZ>("115792089210356248762697446949407573529996955224135760342422259061068512044369");
ZZ G_x = conv<ZZ>("48439561293906451759052585252797914202762949526041747995844080717082404635286");
ZZ G_y = conv<ZZ>("36134250956749795798585127919587881956611106672985015071877198253568414405109");
ZZ G_z = conv<ZZ>("1");
ZZ G_t = conv<ZZ>("69187469364232031836548821531971153808731075654725806004116076052366752432012");

#endif // ECDSA_256

#ifdef RSA_2048

unsigned m = 1031;
unsigned g_GFm = 14;
unsigned nDgts = 3;
unsigned nPrimesByLvl = 200;
unsigned lvl = 4;

ZZ plaintextModulus = conv<ZZ>(
"323170060713110073007148766886699519604441026697154840321303454275246551\
388678908931972014115229134636887179609218980194941195591504909210950881\
523864482831206308773673009960917501977503896521067960576383840675682767\
922186426197561618380943384761704705816458520363050428875758915410658086\
075523991239303855219143333896683424206849747865645694948561760353263220\
580778056593310261927084603141502585928641771167259436037184618573575983\
511523016459044036976132332872312271256847108202097251571017269313234696\
785425806566979350459972683529986382155251663894373355436021354332296046\
45318478604952148193555853611059596231637");

#endif // RSA_2048

int main() {
	SetSeed(ZZ(0));

	FHEcontext context(m, plaintextModulus, g_GFm);
	buildModChain(context, lvl, nDgts, nPrimesByLvl);

	cout << "***************************" << endl
		 << "* Test Import From NFLlib *" << endl
		 << "***************************" << endl;
	cout << "   Parameters" << endl
		 << "---------------------------" << endl
		 << "  m: " << m << endl
		 << "  log2 q: " << context.logOfProduct(context.ctxtPrimes)/log(2) << endl
		 << "  p: " << plaintextModulus << endl
		 << "  nDgts: " << nDgts << endl;

	cout << "===========================" << endl
		 << "   KeyGen"                   << endl;

	FHESecKey secretKey(context);
	const FHEPubKey& publicKey = secretKey;
	secretKey.GenSecKey(512, plaintextModulus); // A Hamming-weight-w secret key

	struct timeval tbeg, tend;
	double texe = 0;
	double tsum = 0;
	long   nexe = 0;
	cout << "===========================" << endl
		 << "   String to DoubleCRT" << endl
		 << "---------------------------" << endl;
	Ctxt* c[50];
	gettimeofday(&tbeg,NULL);
	for (unsigned i = 0; i < 50; i++) {
		c[i] = new Ctxt(publicKey, "c0_out", "c1_out", plaintextModulus);
	}
	gettimeofday(&tend,NULL);
	texe = ((double)(tend.tv_sec-tbeg.tv_sec)) + ((double)(tend.tv_usec-tbeg.tv_usec))/1000000.;
	cout << "  Time: " << texe/50. << " s" << std::endl
	     << "  Size: " << context.logOfProduct(c[0]->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM() << " Mb" << std::endl
	     << "  Flow rate: " << (context.logOfProduct(c[0]->getPrimeSet())/log(2)/1000000.*context.zMStar.getPhiM())/(texe/50.) << " Mbps" << std::endl;
	cout << "===========================" << endl;
	return 0;
}
