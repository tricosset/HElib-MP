#ifndef ELLIPTIC_CURVE_HPP
#define ELLIPTIC_CURVE_HPP

struct ECPoint
{
	Ctxt *X;
	Ctxt *Y;
	Ctxt *Z;
	Ctxt *T;
	ECPoint(const FHEPubKey& publicKey)
	{
		X = new Ctxt(publicKey);
		Y = new Ctxt(publicKey);
		Z = new Ctxt(publicKey);
		T = new Ctxt(publicKey);
	}
	ECPoint(const Ctxt& x, const Ctxt& y, const Ctxt& z, const Ctxt& t, const FHEPubKey& publicKey)
	{
		X = new Ctxt(publicKey);
		Y = new Ctxt(publicKey);
		Z = new Ctxt(publicKey);
		T = new Ctxt(publicKey);

		*X = x;
		*Y = y;
		*Z = z;
		*T = t;
	}
};

#if 0

class ECPrecomputation
{
public:
	ZZ two;
	ZZ two_a;
	ZZ three_a;
	ZZ two_b;
	ZZ twelve_b;
	ZZ four_b;
	ZZ a_sq;
	ZZ three_a_sq;
	ZZ a;
	ZZ two_a_sq;
	ZZ four_a_b;
	ZZ a_3_8_b_2;

	ECPrecomputation()
	{
		two = to_ZZ(2);
		two_a = to_ZZ(2*coeff_a);
		three_a = to_ZZ(3*coeff_a);
		two_b = to_ZZ(2*coeff_b);
		twelve_b = to_ZZ(12*coeff_b);
		four_b = to_ZZ(4*coeff_b);
		a_sq = to_ZZ(coeff_a*coeff_a);
		three_a_sq = to_ZZ(3*coeff_a*coeff_a);
		a = to_ZZ(coeff_a);
		two_a_sq = to_ZZ(2*coeff_a*coeff_a);
		four_a_b = to_ZZ(4*coeff_a*coeff_b);
		a_3_8_b_2 = to_ZZ(coeff_a*coeff_a*coeff_a+8*coeff_b*coeff_b);
	}
};

// C = A+B
void ec_addition(ECPoint& C, const ECPoint& A, const ECPoint& B, const ECPrecomputation& precomputation, const FHEPubKey& publicKey)
{
	Ctxt X1X2 = *A.X, X1Z2 = *A.X, X1T2 = *A.X;
	Ctxt Y1Y2 = *A.Y, Y1Z2 = *A.Y, Y1T2 = *A.Y;
	Ctxt Z1X2 = *A.Z, Z1Y2 = *A.Z, Z1Z2 = *A.Z;
	Ctxt T1X2 = *A.T, T1Y2 = *A.T, T1T2 = *A.T;

	X1X2.multiplyBy(*B.X);
	X1Z2.multiplyBy(*B.Z);
	X1T2.multiplyBy(*B.T);

	Y1Y2.multiplyBy(*B.Y);
	Y1Z2.multiplyBy(*B.Z);
	Y1T2.multiplyBy(*B.T);

	Z1X2.multiplyBy(*B.X);
	Z1Y2.multiplyBy(*B.Y);
	Z1Z2.multiplyBy(*B.Z);

	T1X2.multiplyBy(*B.X);
	T1Y2.multiplyBy(*B.Y);
	T1T2.multiplyBy(*B.T);

	Ctxt tmp(publicKey);

	Ctxt F(publicKey);
	F = T1T2;
		tmp = X1X2;
		tmp.multByConstant(precomputation.two_a);
	F -= tmp;
		tmp = X1Z2;
		tmp += Z1X2;
		tmp.multByConstant(precomputation.four_b);
	F -= tmp;
		tmp = Z1Z2;
		tmp.multByConstant(precomputation.a_sq);
	F += tmp;

	Ctxt H(publicKey);
	H = X1T2;
	H += T1X2;
		tmp = Y1Y2;
		tmp.multByConstant(precomputation.two);
	H += tmp;
		tmp = X1Z2;
		tmp += Z1X2;
		tmp.multByConstant(precomputation.a);
	H += tmp;
		tmp = Z1Z2;
		tmp.multByConstant(precomputation.two_b);
	H += tmp;
	
	Ctxt GA1(publicKey);
	GA1 = T1T2;
		tmp = X1X2;
		tmp.multByConstant(precomputation.two_a);
	GA1 += tmp;
		tmp = Z1Z2;
		tmp.multByConstant(precomputation.three_a_sq);
	GA1 -= tmp;

	Ctxt GA2(publicKey);
	GA2 = GA1;
	
		tmp = X1Z2;
		tmp.multByConstant(precomputation.four_b);
	GA1 += tmp;
		tmp = Z1X2;
		tmp.multByConstant(precomputation.twelve_b);
	GA1 += tmp;

		tmp = Z1X2;
		tmp.multByConstant(precomputation.four_b);
	GA2 += tmp;
		tmp = X1Z2;
		tmp.multByConstant(precomputation.twelve_b);
	GA2 += tmp;

	Ctxt GB1(publicKey);
		tmp = X1T2;
		tmp.multByConstant(precomputation.four_b);
	GB1 = tmp;
		tmp = X1X2;
		tmp.multByConstant(precomputation.two_a_sq);
	GB1 -= tmp;
		tmp = X1Z2;
		tmp.multByConstant(precomputation.four_a_b);
	GB1 -= tmp;
		tmp = T1T2;
		tmp.multByConstant(precomputation.three_a);
	GB1 += tmp;
		tmp = Z1Z2;
		tmp.multByConstant(precomputation.a_3_8_b_2);
	GB1 -= tmp;

	Ctxt GB2(publicKey);
		tmp = T1X2;
		tmp.multByConstant(precomputation.four_b);
	GB2 = tmp;
		tmp = X1X2;
		tmp.multByConstant(precomputation.two_a_sq);
	GB2 -= tmp;
		tmp = Z1X2;
		tmp.multByConstant(precomputation.four_a_b);
	GB2 -= tmp;
		tmp = T1T2;
		tmp.multByConstant(precomputation.three_a);
	GB2 += tmp;
		tmp = Z1Z2;
		tmp.multByConstant(precomputation.a_3_8_b_2);
	GB2 -= tmp;

	Ctxt G1(publicKey);
	G1 = T1Y2;
	G1.multiplyBy(GA1);

	Ctxt G2(publicKey);
	G2 = Z1Y2;
	G2.multiplyBy(GB1);

	Ctxt G3(publicKey);
	G3 = Y1T2;
	G3.multiplyBy(GA2);

	Ctxt G4(publicKey);
	G4 = Y1Z2;
	G4.multiplyBy(GB2);

	*C.X = F;
	C.X->multiplyBy(H);
	
	*C.Y = G1;
	*C.Y += G2;
	*C.Y += G3;
	*C.Y += G4;

	*C.Z = H;
	C.Z->multiplyBy(H);

	*C.T = F;
	C.T->multiplyBy(F);
}
#else

class ECPrecomputation
{
public:
	ZZ two;
	ZZ two_a;
	ZZ three_a;
	ZZ two_b;
	ZZ twelve_b;
	ZZ four_b;
	ZZ a_sq;
	ZZ three_a_sq;
	ZZ a;
	ZZ two_a_sq;
	ZZ four_a_b;
	ZZ a_3_8_b_2;

	ECPrecomputation()
	{
		two = to_ZZ(2);
		two_a = to_ZZ(2*coeff_a);
		three_a = to_ZZ(3*coeff_a);
		two_b = to_ZZ(2*coeff_b);
		twelve_b = to_ZZ(12*coeff_b);
		four_b = to_ZZ(4*coeff_b);
		a_sq = to_ZZ(coeff_a*coeff_a);
		three_a_sq = to_ZZ(3*coeff_a*coeff_a);
		a = to_ZZ(coeff_a);
		two_a_sq = to_ZZ(2*coeff_a*coeff_a);
		four_a_b = to_ZZ(4*coeff_a*coeff_b);
		a_3_8_b_2 = to_ZZ(coeff_a*coeff_a*coeff_a+8*coeff_b*coeff_b);
	}
};

class ECPrecomputationEncrypted
{
public:
	DoubleCRT two;
	DoubleCRT two_a;
	DoubleCRT three_a;
	DoubleCRT two_b;
	DoubleCRT twelve_b;
	DoubleCRT four_b;
	DoubleCRT a_sq;
	DoubleCRT three_a_sq;
	DoubleCRT a;
	DoubleCRT two_a_sq;
	DoubleCRT four_a_b;
	DoubleCRT a_3_8_b_2;

	ECPrecomputationEncrypted(const FHEcontext& context):
	two(context),
	two_a(context),
	three_a(context),
	two_b(context),
	twelve_b(context),
	four_b(context),
	a_sq(context),
	three_a_sq(context),
	a(context),
	two_a_sq(context),
	four_a_b(context),
	a_3_8_b_2(context)
	{
		two = DoubleCRT(ZZX(2), (context));
		two_a = DoubleCRT(ZZX(2*coeff_a), (context));
		three_a = DoubleCRT(ZZX(3*coeff_a), (context));
		two_b = DoubleCRT(ZZX(2*coeff_b), (context));
		twelve_b = DoubleCRT(ZZX(12*coeff_b), (context));
		four_b = DoubleCRT(ZZX(4*coeff_b), (context));
		a_sq = DoubleCRT(ZZX(coeff_a*coeff_a), (context));
		three_a_sq = DoubleCRT(ZZX(3*coeff_a*coeff_a), (context));
		a = DoubleCRT(ZZX(coeff_a), (context));
		two_a_sq = DoubleCRT(ZZX(2*coeff_a*coeff_a), (context));
		four_a_b = DoubleCRT(ZZX(4*coeff_a*coeff_b), (context));
		a_3_8_b_2 = DoubleCRT(ZZX(coeff_a*coeff_a*coeff_a+8*coeff_b*coeff_b), (context));
	}
};

// C = A+B
void ec_addition(ECPoint& C, const ECPoint& A, const ECPoint& B, const ECPrecomputation& precomputation, const FHEPubKey& publicKey)
{
	Ctxt X1X2 = *A.X, X1Z2 = *A.X, X1T2 = *A.X;
	Ctxt Y1Y2 = *A.Y, Y1Z2 = *A.Y, Y1T2 = *A.Y;
	Ctxt Z1X2 = *A.Z, Z1Y2 = *A.Z, Z1Z2 = *A.Z;
	Ctxt T1X2 = *A.T, T1Y2 = *A.T, T1T2 = *A.T;

	X1X2.multiplyBy(*B.X);
	X1Z2.multiplyBy(*B.Z);
	X1T2.multiplyBy(*B.T);

	Y1Y2.multiplyBy(*B.Y);
	// Y1Z2.multiplyBy(*B.Z);
	// Y1T2.multiplyBy(*B.T);

	Z1X2.multiplyBy(*B.X);
	// Z1Y2.multiplyBy(*B.Y);
	Z1Z2.multiplyBy(*B.Z);

	T1X2.multiplyBy(*B.X);
	// T1Y2.multiplyBy(*B.Y);
	T1T2.multiplyBy(*B.T);

	Ctxt tmp(publicKey);

	*C.X = T1T2;
		tmp = X1X2;
		tmp.multByConstant(precomputation.two_a);
	*C.X -= tmp;
		tmp = X1Z2;
		tmp += Z1X2;
		tmp.multByConstant(precomputation.four_b);
	*C.X -= tmp;
		tmp = Z1Z2;
		tmp.multByConstant(precomputation.a_sq);
	*C.X += tmp;

	*C.Z = X1T2;
	*C.Z += T1X2;
		tmp = Y1Y2;
		tmp.multByConstant(precomputation.two);
	*C.Z += tmp;
		tmp = X1Z2;
		tmp += Z1X2;
		tmp.multByConstant(precomputation.a);
	*C.Z += tmp;
		tmp = Z1Z2;
		tmp.multByConstant(precomputation.two_b);
	*C.Z += tmp;
	
	Ctxt GA1(publicKey);
	Ctxt GA2(publicKey);
	GA1 = T1T2;
		tmp = X1X2;
		tmp.multByConstant(precomputation.two_a);
	GA1 += tmp;
		tmp = Z1Z2;
		tmp.multByConstant(precomputation.three_a_sq);
	GA1 -= tmp;

	GA2 = GA1;
	
		tmp = X1Z2;
		tmp.multByConstant(precomputation.four_b);
	GA1 += tmp;
		tmp = Z1X2;
		tmp.multByConstant(precomputation.twelve_b);
	GA1 += tmp;

		tmp = Z1X2;
		tmp.multByConstant(precomputation.four_b);
	GA2 += tmp;
		tmp = X1Z2;
		tmp.multByConstant(precomputation.twelve_b);
	GA2 += tmp;

	Ctxt GB1(publicKey);
	Ctxt GB2(publicKey);

		tmp = T1T2;
		tmp.multByConstant(precomputation.three_a);
	GB1 = tmp;
		tmp = X1X2;
		tmp.multByConstant(precomputation.two_a_sq);
	GB1 -= tmp;
		tmp = Z1Z2;
		tmp.multByConstant(precomputation.a_3_8_b_2);
	GB1 -= tmp;

	GB2 = GB1;

		tmp = X1T2;
		tmp.multByConstant(precomputation.four_b);
	GB1 += tmp;
		tmp = X1Z2;
		tmp.multByConstant(precomputation.four_a_b);
	GB1 -= tmp;

		tmp = T1X2;
		tmp.multByConstant(precomputation.four_b);
	GB2 += tmp;
		tmp = Z1X2;
		tmp.multByConstant(precomputation.four_a_b);
	GB2 -= tmp;


	// T1Y2.multiplyBy(GA1);
	// Z1Y2.multiplyBy(GB1);
	// Y1T2.multiplyBy(GA2);
	// Y1Z2.multiplyBy(GB2);

	T1Y2.multiplyBy2(*B.Y, GA1);
	Z1Y2.multiplyBy2(*B.Y, GB1);
	Y1T2.multiplyBy2(*B.T, GA2);
	Y1Z2.multiplyBy2(*B.Z, GB2);

	*C.Y = T1Y2;
	*C.Y += Z1Y2;
	*C.Y += Y1T2;
	*C.Y += Y1Z2;

	*C.T = *C.X;
	C.T->multiplyBy(*C.T);
	C.X->multiplyBy(*C.Z);
	C.Z->multiplyBy(*C.Z);
}
#endif

#if 0

p = 2^256-2^224+2^192+2^96-1
proof.arithmetic(False) # turn off primality checking
F = GF(p)
n = 115792089210356248762697446949407573529996955224135760342422259061068512044369
a = F(-3)
b = F(41058363725152142129326129780047268409114441015993725554835256314039467401291)
E = EllipticCurve([a, b])
G = E(48439561293906451759052585252797914202762949526041747995844080717082404635286,36134250956749795798585127919587881956611106672985015071877198253568414405109)
G.set_order(115792089210356248762697446949407573529996955224135760342422259061068512044369)

#endif

#endif
