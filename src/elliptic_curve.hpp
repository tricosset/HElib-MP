#ifndef ELLIPTIC_CURVE_HPP
#define ELLIPTIC_CURVE_HPP

struct ECPoint
{
	Ctxt *X;
	Ctxt *Y;
	Ctxt *Z;
	ECPoint(const FHEPubKey& publicKey)
	{
		X = new Ctxt(publicKey);
		Y = new Ctxt(publicKey);
		Z = new Ctxt(publicKey);
	}
	ECPoint(const Ctxt& x, const Ctxt& y, const Ctxt& z, const FHEPubKey& publicKey)
	{
		X = new Ctxt(publicKey);
		Y = new Ctxt(publicKey);
		Z = new Ctxt(publicKey);

		*X = x;
		*Y = y;
		*Z = z;
	}
};

class ECPrecomputationEncrypted
{
public:
	DoubleCRT three;
	DoubleCRT b;
	ECPrecomputationEncrypted(const FHEcontext& context):
	three(context),
	b(context)
	{
		three = DoubleCRT(ZZX(3), (context));
		b = DoubleCRT(ZZX(coeff_b), (context));
	}
};

// C = A+B
void ec_addition(ECPoint& C, const ECPoint& A, const ECPoint& B, const ECPrecomputationEncrypted& precomputation, const FHEPubKey& publicKey)
{
	Ctxt X1X2 = *A.X, X1Y2 = *A.X, X1Z2 = *A.X;
	Ctxt Y1X2 = *A.Y, Y1Y2 = *A.Y, Y1Z2 = *A.Y;
	Ctxt Z1X2 = *A.Z, Z1Y2 = *A.Z, Z1Z2 = *A.Z;

	X1X2.multiplyBy(*B.X);
	X1Y2.multiplyBy(*B.Y);
	X1Z2.multiplyBy(*B.Z);

	Y1X2.multiplyBy(*B.X);
	Y1Y2.multiplyBy(*B.Y);
	Y1Z2.multiplyBy(*B.Z);

	Z1X2.multiplyBy(*B.X);
	Z1Y2.multiplyBy(*B.Y);
	Z1Z2.multiplyBy(*B.Z);

	Ctxt T0(publicKey), T1(publicKey), T2(publicKey);
	Ctxt T3(publicKey), T4(publicKey), T5(publicKey);

	Ctxt tmp(publicKey);

	T0 = X1Y2;
		T0 += Y1X2;
	T1 = Y1Y2;
	T2 = X1Z2;
		T2 += Z1X2;
		tmp = Z1Z2;
		tmp.multByConstant(precomputation.b);
		T2 -= tmp;
		T2.multByConstant(precomputation.three);
	T3 = Y1Z2;
		T3 += Z1Y2;
	T4 = X1Z2;
		T4 += Z1X2;
		T4.multByConstant(precomputation.b);
		tmp = Z1Z2;
		tmp.multByConstant(precomputation.three);
		T4 -= X1X2;
		T4 -= tmp;
	T5 = X1X2;
		T5 -= Z1Z2;
		T5.multByConstant(precomputation.three);

	*C.X = T1;
		*C.X += T2;
		C.X->multiplyBy(T0);
		tmp = T3;
		tmp.multiplyBy(T4);
		tmp.multByConstant(precomputation.three);
		*C.X -= tmp;

	*C.Y = T1;
		*C.Y -= T2;
		tmp = T1;
		tmp += T2;
		C.Y->multiplyBy(tmp);
		tmp = T4;
		tmp.multiplyBy(T5);
		tmp.multByConstant(precomputation.three);
		*C.Y += tmp;

	*C.Z = T1;
		*C.Z -= T2;
		C.Z->multiplyBy(T3);
		tmp = T0;
		tmp.multiplyBy(T5);
		*C.Z += tmp;
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
