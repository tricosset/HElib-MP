/* Copyright (C) 2012,2013 IBM Corp.
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
/* CModulus.cpp - supports forward and backward length-m FFT transformations
 *
 * This is a wrapper around the bluesteinFFT routines, for one modulus q.
 *
 * On initialization, it initizlies NTL's zz_pContext for this q
 * and computes a 2m-th root of unity r mod q and also r^{-1} mod q.
 * Thereafter this class provides FFT and iFFT routines that converts between
 * time & frequency domains. Some tables are computed the first time that
 * each directions is called, which are then used in subsequent computations.
 * 
 * The "time domain" polynomials are represented as ZZX, which are reduced
 * modulo Phi_m(X). The "frequency domain" are jusr vectors of integers
 * (vec_long), that store only the evaluation in primitive m-th
 * roots of unity.
 */

#include "CModulus.h"
#include "timing.h"


// It is assumed that m,q,context, and root are already set. If root is set
// to zero, it will be computed by the compRoots() method. Then rInv is
// computed as the inverse of root.


zz_pContext BuildContext(long p, long maxroot) {
   if (maxroot <= CalcMaxRoot(p))
      return zz_pContext(INIT_USER_FFT, p);
   else
      return zz_pContext(p, maxroot);
}


// Constructor: it is assumed that zms is already set with m>1
// If q == 0, then the current context is used
Cmodulus::Cmodulus(const PAlgebra &zms, long qq, long rt)
{
  assert(zms.getM()>1);
  bool explicitModulus = true;

  if (qq == 0) {
    q = zz_p::modulus();
    explicitModulus = false;
  }
  else
    q = qq;

  zMStar = &zms;
  root = rt;

  long mm;
  mm = zms.getM();
  m_inv = InvMod(mm, q);

  zz_pBak bak; 

  if (explicitModulus) {
    bak.save(); // backup the current modulus
    context = BuildContext(q, NextPowerOfTwo(zms.getM()) + 1);
    context.restore();       // set NTL's current modulus to q
  }
  else
    context.save();

  if (root==0) { // Find a 2m-th root of unity modulo q, if not given
    zz_p rtp;
    long e = 2*zms.getM();
    FindPrimitiveRoot(rtp,e); // NTL routine, relative to current modulus
    if (rtp==0) // sanity check
      Error("Cmod::compRoots(): no 2m'th roots of unity mod q");
    root = rep(rtp);
  }
  rInv = InvMod(root,q); // set rInv = root^{-1} mod q

  // Allocate memory (relative to current modulus that was defined above).
  // These objects will be initialized when anyone calls FFT/iFFT.

  zz_pX phimx_poly;
  conv(phimx_poly, zms.getPhimX());

  powers.set_ptr(new zz_pX);
  Rb.set_ptr(new fftRep);
  ipowers.set_ptr(new zz_pX);
  iRb.set_ptr(new fftRep);
  phimx.set_ptr(new zz_pXModulus1(zms.getM(), phimx_poly));

  BluesteinInit(mm, conv<zz_p>(root), *powers, powers_aux, *Rb);
  BluesteinInit(mm, conv<zz_p>(rInv), *ipowers, ipowers_aux, *iRb);
}

Cmodulus& Cmodulus::operator=(const Cmodulus &other)
{
  if (this == &other) return *this;

  zMStar  =  other.zMStar; // Yes, really copy this pointer
  q       = other.q;
  m_inv   = other.m_inv;

  context = other.context;
  zz_pBak bak; bak.save(); // backup the current modulus
  context.restore();       // Set NTL's current modulus to q

  // NOTE: newer versions of NTL allow fftRep's and zz_pXModulus's to be copied
  // "out of context" (versions after 7.0.*). However, those copies
  // are not intended to allow copies out of one context into another,
  // so we still need to use copied_ptr's (but not context restoration).
  // All of this is fairly academic, as I don't think we really
  // copy FHEcontexts around anywhere. Also, it would be cleaner
  // to make the vector in FHEcontext be a vector of copied_ptr<Cmodulus>

  root = other.root;
  rInv = other.rInv;

  powers_aux = other.powers_aux;
  ipowers_aux = other.ipowers_aux;

  // copy data, not pointers in these fields
  powers = other.powers;
  Rb = other.Rb;
  ipowers = other.ipowers;
  iRb = other.iRb;
  phimx = other.phimx;



  return *this;
}

void Cmodulus::FFT(vec_long &y, const ZZX& x) const
{
  FHE_TIMER_START;
  zz_pBak bak; bak.save();
  context.restore();
  zz_p rt;
  zz_pX& tmp = Cmodulus::getScratch_zz_pX();

  conv(tmp,x);      // convert input to zpx format
  conv(rt, root);  // convert root to zp format

  BluesteinFFT(tmp, getM(), rt, *powers, powers_aux, *Rb); // call the FFT routine

  // copy the result to the output vector y, keeping only the
  // entries corresponding to primitive roots of unity
  y.SetLength(zMStar->getPhiM());
  long i,j;
  long m = getM();
  for (i=j=0; i<m; i++)
    if (zMStar->inZmStar(i)) y[j++] = rep(coeff(tmp,i));
}


void Cmodulus::iFFT(zz_pX &x, const vec_long& y)const
{
  FHE_TIMER_START;
  zz_pBak bak; bak.save();
  context.restore();
  zz_p rt;

  long m = getM();

  // convert input to zpx format, initializing only the coeffs i s.t. (i,m)=1
  x.rep.SetLength(m);
  long i,j;
  for (i=j=0; i<m; i++)
    if (zMStar->inZmStar(i)) x.rep[i].LoopHole() = y[j++]; // DIRT: y[j] already reduced
  x.normalize();
  conv(rt, rInv);  // convert rInv to zp format

  BluesteinFFT(x, m, rt, *ipowers, ipowers_aux, *iRb); // call the FFT routine

  // reduce the result mod (Phi_m(X),q) and copy to the output polynomial x
  { FHE_NTIMER_START(iFFT_division);
    rem(x, x, *phimx); // out %= (Phi_m(X),q)
  }

  // normalize
  zz_p mm_inv;
  conv(mm_inv, m_inv);
  x *= mm_inv; 
}


zz_pX& Cmodulus::getScratch_zz_pX() 
{
   NTL_THREAD_LOCAL static zz_pX scratch;
   return scratch;
}


fftRep& Cmodulus::getScratch_fftRep(long k)
{
  NTL_THREAD_LOCAL static fftRep rep;
  NTL_THREAD_LOCAL static long MaxK[4] = {-1, -1, -1, -1};

  long NumPrimes = zz_pInfo->NumPrimes;

  for (long i = 0; i < NumPrimes; i++) {
    if (k > MaxK[i]) {
      rep.tbl[i].SetLength(1L << k);
      MaxK[i] = k;
    }
  }

  rep.NumPrimes = NumPrimes;
  rep.k = rep.MaxK = k;

  return rep;
}




