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

#include "Ctxt.h"
#include "FHE.h"
#include "timing.h"

#ifndef BIG_P
// Dummy encryption, this procedure just encodes the plaintext in a Ctxt object
void Ctxt::DummyEncrypt(const ZZX& ptxt, double size)
{
  if (size < 0.0) {
    size = ((double) context.zMStar.getPhiM()) * ptxtSpace*ptxtSpace /12.0;
  }
  noiseVar = size;
  primeSet = context.ctxtPrimes;

  // A single part, with the plaintext as data and handle pointing to 1

  long f = rem(context.productOfPrimes(context.ctxtPrimes),ptxtSpace);
  if (f == 1) { // scale by constant
    DoubleCRT dcrt(ptxt, context, primeSet);  
    parts.assign(1, CtxtPart(dcrt));
  } else {
    ZZX tmp;
    MulMod(tmp, ptxt, f, ptxtSpace, /*positive=*/false);
    DoubleCRT dcrt(tmp, context, primeSet);  
    parts.assign(1, CtxtPart(dcrt));
  }
}
#else
// Dummy encryption, this procedure just encodes the plaintext in a Ctxt object
void Ctxt::DummyEncrypt(const ZZX& ptxt, double size)
{
  if (size < 0.0) {
	size = ((double) context.zMStar.getPhiM()) * to_double(ptxtSpace)*to_double(ptxtSpace) /12.0;
  }
  noiseVar = size;
  primeSet = context.ctxtPrimes;

  // A single part, with the plaintext as data and handle pointing to 1

  ZZ f;
  rem(f, context.productOfPrimes(context.ctxtPrimes),ptxtSpace);
  if (f == 1) { // scale by constant
    DoubleCRT dcrt(ptxt, context, primeSet);
    parts.assign(1, CtxtPart(dcrt));
  } else {
    ZZX tmp;
    MulMod(tmp, ptxt, f, ptxtSpace, /*positive=*/false);
    DoubleCRT dcrt(tmp, context, primeSet);
    parts.assign(1, CtxtPart(dcrt));
  }
}
#endif


// Sanity-check: Check that prime-set is "valid", i.e. that it 
// contains either all the special primes or none of them
bool Ctxt::verifyPrimeSet() const
{
  IndexSet s = primeSet & context.specialPrimes; // special primes in primeSet
  if (!empty(s) && s!=context.specialPrimes) return false;

  s = primeSet / s;                              // ctxt primes in primeSet
  return (s.isInterval() && s.first()<=1 && !empty(s));
}

#ifndef BIG_P
//! @brief How many levels in the "base-set" for that ciphertext
long Ctxt::findBaseLevel() const 
{
  IndexSet s;
  findBaseSet(s);

  if (context.containsSmallPrime()) {
    if (s.contains(context.ctxtPrimes.first()))
      return 2*card(s) -1; // 1st prime is half size
    else
      return 2*card(s);
  }
  else return card(s);     // one prime per level
}
#else
//! @brief How many levels in the "base-set" for that ciphertext
long Ctxt::findBaseLevel() const
{
	  IndexSet s;
	  findBaseSet(s);
  if (context.containsSmallPrime()) {
    if (s.contains(context.ctxtPrimes.first()))
      return 2*card(s) -1; // 1st prime is half size
    else
      return 2*card(s);
  }
  else return card(s);     // one prime per level
}
#endif

bool CtxtPart::operator==(const CtxtPart& other) const
{
  if (((DoubleCRT&)*this)!=((DoubleCRT&)other)) return false;
  
  return (skHandle==other.skHandle);
}

// Checking equality between ciphertexts. This routine performs a
// "shallow" check, comparing only pointers to ciphertext parts.
bool Ctxt::equalsTo(const Ctxt& other, bool comparePkeys) const
{
  if (&context != &other.context) return false;
  if (comparePkeys && &pubKey != &other.pubKey) return false;

  if (parts.size() != other.parts.size()) return false;
  for (size_t i=0; i<parts.size(); i++)
    if (parts[i] != other.parts[i]) return false;

  if (primeSet != other.primeSet) return false;
  if (ptxtSpace != other.ptxtSpace) return false;

  // compare noiseVar, ignoring small deviations
  if (noiseVar == 0.0) return (other.noiseVar == 0.0);
  xdouble ratio = other.noiseVar / noiseVar;
  return (ratio>0.9 && ratio<1.1);
}

#ifndef BIG_P
// Constructor
Ctxt::Ctxt(const FHEPubKey& newPubKey, long newPtxtSpace):
  context(newPubKey.getContext()), pubKey(newPubKey), ptxtSpace(newPtxtSpace),
  noiseVar(to_xdouble(0.0))
{
  if (ptxtSpace<=0) ptxtSpace = pubKey.getPtxtSpace();
  else assert (GCD(ptxtSpace, pubKey.getPtxtSpace()) > 1); // sanity check
  primeSet=context.ctxtPrimes;
}
#else
// Constructor
Ctxt::Ctxt(const FHEPubKey& newPubKey, ZZ newPtxtSpace):
  context(newPubKey.getContext()), pubKey(newPubKey), ptxtSpace(newPtxtSpace),
  noiseVar(to_xdouble(0.0))
{
  if (ptxtSpace<=0) ptxtSpace = pubKey.getPtxtSpace();
  else assert (GCD(ptxtSpace, pubKey.getPtxtSpace()) > 1); // sanity check
  primeSet=context.ctxtPrimes;
}
#endif

#ifdef BIG_P
Ctxt::Ctxt(const FHEPubKey& newPubKey, const char* file_c0, const char* file_c1, ZZ newPtxtSpace): // constructor
  context(newPubKey.getContext()), pubKey(newPubKey), ptxtSpace(newPtxtSpace),
  noiseVar(to_xdouble(0.0))
{
  if (ptxtSpace<=0) ptxtSpace = pubKey.getPtxtSpace();
  else assert (GCD(ptxtSpace, pubKey.getPtxtSpace()) > 1); // sanity check
  primeSet=context.ctxtPrimes;
  pubKey.NFLlib2HElib(*this, file_c0, file_c1, ptxtSpace);
}
#endif

// Constructor
Ctxt::Ctxt(ZeroCtxtLike_type, const Ctxt& ctxt):
  context(ctxt.getPubKey().getContext()), pubKey(ctxt.getPubKey()), 
  ptxtSpace(ctxt.getPtxtSpace()),
  noiseVar(to_xdouble(0.0))
{
  // same body as previous constructor
  if (ptxtSpace<=0) ptxtSpace = pubKey.getPtxtSpace();
  else assert (GCD(ptxtSpace, pubKey.getPtxtSpace()) > 1); // sanity check
  primeSet=context.ctxtPrimes;
}


// A private assignment method that does not check equality of context or
// public key, this needed for example when we copy the pubEncrKey member
// between different public keys.
Ctxt& Ctxt::privateAssign(const Ctxt& other)
{
  if (this == &other) return *this; // both point to the same object

  parts = other.parts;
  primeSet = other.primeSet;
  ptxtSpace = other.ptxtSpace;
  noiseVar  = other.noiseVar;
  return *this;
}

// Ciphertext maintenance

// mod-switch up to add the primes in s \setminus primeSet, after this call we
// have s<=primeSet. s must contain either all special primes or none of them.
void Ctxt::modUpToSet(const IndexSet &s)
{
#ifdef VERBOSE
	std::cout << "Ctxt::modUpToSet" << std::endl
			  << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl
			  ;
#endif
  //  FHE_TIMER_START;
  IndexSet setDiff = s/primeSet; // set minus (primes in s but not in primeSet)
  if (empty(setDiff)) return;    // nothing to do, no primes are added

  // scale up all the parts to use also the primes in setDiff
  double f = 0.0;
  for (long i=0; i<lsize(parts); i++) {
    // addPrimesAndScale returns the logarithm of the product of added primes,
    // all calls should return the same value = log(prod. of primes in setDiff)
    f = parts[i].addPrimesAndScale(setDiff);
  }

  // The variance estimate grows by a factor of exp(f)^2 = exp(2f)
  noiseVar *= xexp(2*f);
#ifdef VERBOSE
	std::cout << "\tnoiseVar *= xexp(2*f);" << std::endl
			  << "\tf: " << f << std::endl
			  << "\txexp(2*f): " << xexp(2*f) << std::endl
			  << "\tnoiseVar: " << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif

  primeSet.insert(setDiff); // add setDiff to primeSet
  assert(verifyPrimeSet()); // sanity-check: ensure primeSet is still valid
  //  FHE_TIMER_STOP;
}

#ifndef BIG_P
// mod-switch down to primeSet \intersect s, after this call we have
// primeSet<=s. s must contain either all special primes or none of them.
void Ctxt::modDownToSet(const IndexSet &s)
{
  IndexSet intersection = primeSet & s;
  //  assert(!empty(intersection));       // some primes must be left
  if (empty(intersection)) {
    cerr << "modDownToSet called from "<<primeSet<<" to "<<s<<endl;
    exit(1);
  }
  if (intersection==primeSet) return; // nothing to do, removing no primes
  FHE_TIMER_START;

  IndexSet setDiff = primeSet / intersection; // set-minus

  // Scale down all the parts: use either a simple "drop down" (just removing
  // primes, i.e., reducing the ctxt modulo the samaller modulus), or a "real
  // modulus switching" with rounding, basically whichever yeilds smaller
  // noise. Recall that we keep the invariant that a ciphertext mod Q is
  // decrypted to Q*m (mod p), so if we just "drop down" we still need to
  // multiply by (Q^{-1} mod p).

  // Get an estimate for the added noise term for modulus switching
  xdouble addedNoiseVar = modSwitchAddedNoiseVar();
  if (noiseVar*ptxtSpace*ptxtSpace < addedNoiseVar) {     // just "drop down"
    long prodInv = InvMod(rem(context.productOfPrimes(setDiff),ptxtSpace), ptxtSpace);

    for (size_t i=0; i<parts.size(); i++) {
      parts[i].removePrimes(setDiff);         // remove the primes not in s
      parts[i] *= prodInv;
      // WARNING: the following line is written just so to prevent overflow
      noiseVar = noiseVar*prodInv*prodInv;
    }
    //    cerr << "DEGENERATE DROP\n";
  } 
  else {                                      // do real mod switching
    for (size_t i=0; i<parts.size(); i++) 
      parts[i].scaleDownToSet(intersection, ptxtSpace);

    // update the noise estimate
    double f = context.logOfProduct(setDiff);
    noiseVar /= xexp(2*f);
    noiseVar += addedNoiseVar;
  }
  primeSet.remove(setDiff); // remove the primes not in s
  assert(verifyPrimeSet()); // sanity-check: ensure primeSet is still valid
  FHE_TIMER_STOP;
}
#else
// mod-switch down to primeSet \intersect s, after this call we have
// primeSet<=s. s must contain either all special primes or none of them.
void Ctxt::modDownToSet(const IndexSet &s)
{
#ifdef VERBOSE
	std::cout << "Ctxt::modDownToSet:" << std::endl;
	  std::cout << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif

	  IndexSet intersection = primeSet & s;
	  //  assert(!empty(intersection));       // some primes must be left
	  if (empty(intersection)) {
	    cerr << "modDownToSet called from "<<primeSet<<" to "<<s<<endl;
	    exit(1);
	  }
	  if (intersection==primeSet) return; // nothing to do, removing no primes
	  FHE_TIMER_START;

	  IndexSet setDiff = primeSet / intersection; // set-minus

	  // Scale down all the parts: use either a simple "drop down" (just removing
	  // primes, i.e., reducing the ctxt modulo the samaller modulus), or a "real
	  // modulus switching" with rounding, basically whichever yeilds smaller
	  // noise. Recall that we keep the invariant that a ciphertext mod Q is
	  // decrypted to Q*m (mod p), so if we just "drop down" we still need to
	  // multiply by (Q^{-1} mod p).

	  // Get an estimate for the added noise term for modulus switching
	  xdouble addedNoiseVar = modSwitchAddedNoiseVar();
#ifdef VERBOSE
	  std::cout << "\taddedNoiseVar: " << addedNoiseVar << std::endl;
#endif
	  if (noiseVar*to_xdouble(ptxtSpace)*to_xdouble(ptxtSpace) < addedNoiseVar) {     // just "drop down"
#ifdef VERBOSE
		  std::cout << "\tnoiseVar*ptxtSpace*ptxtSpace < addedNoiseVar" << std::endl;
#endif

	    ZZ prodInv;
	    InvMod(prodInv, context.productOfPrimes(setDiff)%ptxtSpace, ptxtSpace);
#ifdef VERBOSE
		  std::cout << "\tprodInv: " << prodInv << std::endl;
#endif

	    for (size_t i=0; i<parts.size(); i++) {
	      parts[i].removePrimes(setDiff);         // remove the primes not in s
	      parts[i] *= prodInv;
	      // WARNING: the following line is written just so to prevent overflow
	      noiseVar = noiseVar*to_xdouble(prodInv)*to_xdouble(prodInv);
#ifdef VERBOSE
	      std::cout << "\tnoiseVar = noiseVar*to_xdouble(prodInv)*to_xdouble(prodInv);" << std::endl
	    		    << "\tto_xdouble(prodInv): " << to_xdouble(prodInv) << std::endl
	    		    << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif

	    }
	    //    cerr << "DEGENERATE DROP\n";
	  }
	  else {                                      // do real mod switching
#ifdef VERBOSE
		  std::cout << "\tnoiseVar*ptxtSpace*ptxtSpace >= addedNoiseVar" << std::endl;
#endif
	    for (size_t i=0; i<parts.size(); i++)
	      parts[i].scaleDownToSet(intersection, ptxtSpace);

	    // update the noise estimate
	    double f = context.logOfProduct(setDiff);
	    noiseVar /= xexp(2*f);
	    noiseVar += addedNoiseVar;
#ifdef VERBOSE
	      std::cout << "\tdouble f = context.logOfProduct(setDiff);" << std::endl
	    		    << "\tnoiseVar /= xexp(2*f);" << std::endl
	    		    << "\tnoiseVar += addedNoiseVar;" << std::endl
	    		    << "\tf: " << f << std::endl
	    		    << "\txexp(2*f): " << xexp(2*f) << std::endl
	    		    << "\taddedNoiseVar: " << addedNoiseVar << std::endl
	    		    << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif
	  }
	  primeSet.remove(setDiff); // remove the primes not in s
	  assert(verifyPrimeSet()); // sanity-check: ensure primeSet is still valid
	  FHE_TIMER_STOP;
}
#endif

// Modulus-switching down
void Ctxt::modDownToLevel(long lvl)
{
  long currentLvl;
  IndexSet targetSet;
  IndexSet currentSet = primeSet & context.ctxtPrimes;
  if (context.containsSmallPrime()) {
    currentLvl = 2*card(currentSet);
    if (currentSet.contains(0))
      currentLvl--;  // first prime is half the size

    if (lvl & 1) {   // odd level, includes the half-size prime
      targetSet = IndexSet(0,(lvl-1)/2);
    } else {
      targetSet = IndexSet(1,lvl/2);
    }
  }
  else {
    currentLvl = card(currentSet);
    targetSet = IndexSet(0,lvl-1);    // one prime per level
  }

  // If target is not below the current level, nothing to do
  if (lvl >= currentLvl && currentSet==primeSet) return;

  if (lvl >= currentLvl) { // just remove the special primes
    targetSet = currentSet;
  }

  // sanity-check: interval does not contain special primes
  assert(targetSet.disjointFrom(context.specialPrimes));

  // may need to mod-UP to include the smallest prime
  if (targetSet.contains(0) && !currentSet.contains(0))
    modUpToSet(targetSet); // adds the primes in targetSet / primeSet

  modDownToSet(targetSet); // removes the primes in primeSet / targetSet
}

void Ctxt::blindCtxt(const ZZX& poly)
{
  Ctxt tmp(pubKey);
  pubKey.Encrypt(tmp,poly,ptxtSpace,/*highNoise=*/true);
  *this += tmp;
  // FIXME: This implementation is not optimized, the levels in the
  //  modulus chain should be handled much better.
}

#ifndef BIG_P
// Reduce plaintext space to a divisor of the original plaintext space
void Ctxt::reducePtxtSpace(long newPtxtSpace)
{
  long g = GCD(ptxtSpace, newPtxtSpace);
  assert (g>1);
  ptxtSpace = g;
}
#else
// Reduce plaintext space to a divisor of the original plaintext space
void Ctxt::reducePtxtSpace(ZZ newPtxtSpace)
{
  ZZ g;
  GCD(g, ptxtSpace, newPtxtSpace);
  assert (g>1);
  ptxtSpace = g;
}
#endif

#ifndef BIG_P
// key-switch to (1,s_i), s_i is the base key with index keyID. If
// keyID<0 then re-linearize to any key for which a switching matrix exists
void Ctxt::reLinearize(long keyID)
{
  // Special case: if *this is empty or already re-linearized then do nothing
  if (this->isEmpty() || this->inCanonicalForm(keyID)) return;

  FHE_TIMER_START;
  this->reduce();

  // To relinearize, the primeSet must be disjoint from the special primes
  if (!primeSet.disjointFrom(context.specialPrimes)) {
    modDownToSet(primeSet / context.specialPrimes);
    // cout << "<<< special primes\n";
  }

  long g = ptxtSpace;
  Ctxt tmp(pubKey, ptxtSpace); // an empty ciphertext, same plaintext space

  double logProd = context.logOfProduct(context.specialPrimes);
  tmp.noiseVar = noiseVar * xexp(2*logProd); // The noise after mod-UP

  for (size_t i=0; i<parts.size(); i++) {
    CtxtPart& part  = parts[i];

    // For a part relative to 1 or base,  only scale and add
    if (part.skHandle.isOne() || part.skHandle.isBase(keyID)) {
      part.addPrimesAndScale(context.specialPrimes);
      tmp.addPart(part, /*matchPrimeSet=*/true);
      continue;
    }
    // Look for a key-switching matrix to re-linearize this part
    const KeySwitch& W = (keyID>=0)? 
      pubKey.getKeySWmatrix(part.skHandle,keyID) :
      pubKey.getAnyKeySWmatrix(part.skHandle);

    assert(W.toKeyID>=0);       // verify that a switching matrix exists

    g = GCD(W.ptxtSpace, g);    // verify that the plaintext spaces match
    assert (g>1);
    tmp.ptxtSpace = g;

    tmp.keySwitchPart(part, W); // switch this part & update noiseVar
  }
  *this = tmp;

#if 0
  // alternative strategy: get rid of special primes
  // after each reLin...
  if (!primeSet.disjointFrom(context.specialPrimes)) {
    modDownToSet(primeSet / context.specialPrimes);
    // cout << ">>> special primes\n";
  }
#endif

  FHE_TIMER_STOP;
}
#else
// key-switch to (1,s_i), s_i is the base key with index keyID. If
// keyID<0 then re-linearize to any key for which a switching matrix exists
void Ctxt::reLinearize(long keyID)
{
#ifdef VERBOSE
	std::cout << "Ctxt::reLinearize" << std::endl
			  << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif
  // Special case: if *this is empty or already re-linearized then do nothing
  if (this->isEmpty() || this->inCanonicalForm(keyID)) return;

  FHE_TIMER_START;
  this->reduce();

  // To relinearize, the primeSet must be disjoint from the special primes
  if (!primeSet.disjointFrom(context.specialPrimes)) {
	modDownToSet(primeSet / context.specialPrimes);
	// cout << "<<< special primes\n";
  }

  ZZ g = ptxtSpace;
  Ctxt tmp(pubKey, ptxtSpace); // an empty ciphertext, same plaintext space

  double logProd = context.logOfProduct(context.specialPrimes);

  tmp.noiseVar = noiseVar * xexp(2*logProd); // The noise after mod-UP
#ifdef VERBOSE
	std::cout << "\ttmp.noiseVar = noiseVar * xexp(2*logProd);" << std::endl
			  << "\tdouble logProd = context.logOfProduct(context.specialPrimes);" << std::endl
			  << "\tlogProd: " << logProd << std::endl
			  << "\txexp(2*logProd): " << xexp(2*logProd) << std::endl
			  << "\ttmp.noiseVar: " << tmp.noiseVar << " (" << log(tmp.noiseVar)/log(2)/2 << ")" << std::endl;
#endif
  for (size_t i=0; i<parts.size(); i++) {
    CtxtPart& part  = parts[i];

    // For a part relative to 1 or base,  only scale and add
    if (part.skHandle.isOne() || part.skHandle.isBase(keyID)) {
      part.addPrimesAndScale(context.specialPrimes);
      tmp.addPart(part, /*matchPrimeSet=*/true);
      continue;
    }
    // Look for a key-switching matrix to re-linearize this part
    const KeySwitch& W = (keyID>=0)?
      pubKey.getKeySWmatrix(part.skHandle,keyID) :
      pubKey.getAnyKeySWmatrix(part.skHandle);

    assert(W.toKeyID>=0);       // verify that a switching matrix exists

    g = GCD(W.ptxtSpace, g);    // verify that the plaintext spaces match
    assert (g>1);
    tmp.ptxtSpace = g;

    tmp.keySwitchPart(part, W); // switch this part & update noiseVar
  }
  *this = tmp;
#ifdef VERBOSE
std::cout << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif
#if 0
  // alternative strategy: get rid of special primes
  // after each reLin...
  if (!primeSet.disjointFrom(context.specialPrimes)) {
    modDownToSet(primeSet / context.specialPrimes);
    // cout << ">>> special primes\n";
  }
#endif

  FHE_TIMER_STOP;
}
#endif

void Ctxt::cleanUp()
{
  reLinearize();
  reduce();
  if (!primeSet.disjointFrom(context.specialPrimes)) {
    modDownToSet(primeSet / context.specialPrimes);
  }
}

#ifndef BIG_P
// Takes as arguments a ciphertext-part p relative to s' and a key-switching
// matrix W = W[s'->s], uses W to switch p relative to (1,s), and adds the
// result to *this.
// It is assumed that the part p does not include any of the special primes,
// and that if *this is not an empty ciphertext then its primeSet is
// p.getIndexSet() \union context.specialPrimes
void Ctxt::keySwitchPart(const CtxtPart& p, const KeySwitch& W)
{
  // no special primes in the input part
  assert(context.specialPrimes.disjointFrom(p.getIndexSet()));

  // For parts p that point to 1 or s, only scale and add
  if (p.skHandle.isOne() || p.skHandle.isBase(W.toKeyID)) { 
    CtxtPart pp = p;
    pp.addPrimesAndScale(context.specialPrimes);
    addPart(pp, /*matchPrimeSet=*/true);
    return; 
  }

  // some sanity checks
  assert(W.fromKey == p.skHandle);  // the handles must match

  // Compute the number of digits that we need and the esitmated added noise
  // from switching this ciphertext part.
  long pSpace = W.ptxtSpace;
  long nDigits = 0;
  xdouble addedNoise = to_xdouble(0.0);
  double sizeLeft = context.logOfProduct(p.getIndexSet());
  for (size_t i=0; i<context.digits.size() && sizeLeft>0.0; i++) {    
    nDigits++;

    double digitSize = context.logOfProduct(context.digits[i]);
    if (sizeLeft<digitSize) digitSize=sizeLeft; // need only part of this digit

    // Added noise due to this digit is phi(m) * sigma^2 * pSpace^2 * |Di|^2/4, 
    // where |Di| is the magnitude of the digit

    // WARNING: the following line is written just so to prevent overflow
    addedNoise += to_xdouble(context.zMStar.getPhiM()) * pSpace*pSpace
      * xexp(2*digitSize) * context.stdev*context.stdev / 4.0;

    sizeLeft -= digitSize;
  }

  // Sanity-check: make sure that the added noise is not more than the special
  // primes can handle: After dividing the added noise by the product of all
  // the special primes, it should be smaller than the added noise term due
  // to modulus switching, i.e., keyWeight * phi(m) * pSpace^2 / 12

  long keyWeight = pubKey.getSKeyWeight(p.skHandle.getSecretKeyID());
  double phim = context.zMStar.getPhiM();
  double logModSwitchNoise = log((double)keyWeight) 
    +2*log((double)pSpace) +log(phim) -log(12.0);
  double logKeySwitchNoise = log(addedNoise) 
    -2*context.logOfProduct(context.specialPrimes);
  assert(logKeySwitchNoise < logModSwitchNoise);

  // Break the ciphertext part into digits, if needed, and scale up these
  // digits using the special primes. This is the most expensive operation
  // during homormophic evaluation, so it should be thoroughly optimized.

  vector<DoubleCRT> polyDigits;
  p.breakIntoDigits(polyDigits, nDigits);

  // Finally we multiply the vector of digits by the key-switching matrix

  // An object to hold the pseudorandom ai's, note that it must be defined
  // with the maximum number of levels, else the PRG will go out of synch.
  // FIXME: This is a bug waiting to happen.
  DoubleCRT ai(context);

  // Set the first ai using the seed, subsequent ai's (if any) will
  // use the evolving RNG state (NOTE: this is not thread-safe)

  RandomState state;
  SetSeed(W.prgSeed);

  // Add the columns in, one by one
  DoubleCRT tmp(context, IndexSet::emptySet());
  
  for (unsigned long i=0; i<polyDigits.size(); i++) {
    ai.randomize();
    tmp = polyDigits[i];
  
    // The operations below all use the IndexSet of tmp
  
    // add part*a[i] with a handle pointing to base of W.toKeyID
    tmp.Mul(ai,  /*matchIndexSet=*/false);
    addPart(tmp, SKHandle(1,1,W.toKeyID), /*matchPrimeSet=*/true);
  
    // add part*b[i] with a handle pointing to one
    polyDigits[i].Mul(W.b[i], /*matchIndexSet=*/false);
    addPart(polyDigits[i], SKHandle(), /*matchPrimeSet=*/true);
  }
  noiseVar += addedNoise;
} // restore random state upon destruction of the RandomState, see NumbTh.h
#else
// Takes as arguments a ciphertext-part p relative to s' and a key-switching
// matrix W = W[s'->s], uses W to switch p relative to (1,s), and adds the
// result to *this.
// It is assumed that the part p does not include any of the special primes,
// and that if *this is not an empty ciphertext then its primeSet is
// p.getIndexSet() \union context.specialPrimes
void Ctxt::keySwitchPart(const CtxtPart& p, const KeySwitch& W)
{
#ifdef VERBOSE
	std::cout << "Ctxt::keySwitchPart" << std::endl;
#endif

  // no special primes in the input part
  assert(context.specialPrimes.disjointFrom(p.getIndexSet()));

  // For parts p that point to 1 or s, only scale and add
  if (p.skHandle.isOne() || p.skHandle.isBase(W.toKeyID)) {
	CtxtPart pp = p;
	pp.addPrimesAndScale(context.specialPrimes);
	addPart(pp, /*matchPrimeSet=*/true);
	return;
  }

  // some sanity checks
  assert(W.fromKey == p.skHandle);  // the handles must match

  // Compute the number of digits that we need and the esitmated added noise
  // from switching this ciphertext part.
  ZZ pSpace = W.ptxtSpace;
  long nDigits = 0;
  xdouble addedNoise = to_xdouble(0.0);

  double sizeLeft = context.logOfProduct(p.getIndexSet());
  for (size_t i=0; i<context.digits.size() && sizeLeft>0.0; i++) {
	nDigits++;

	double digitSize = context.logOfProduct(context.digits[i]);

	if (sizeLeft<digitSize) digitSize=sizeLeft; // need only part of this digit

	// Added noise due to this digit is phi(m) * sigma^2 * pSpace^2 * |Di|^2/4,
	// where |Di| is the magnitude of the digit

	// WARNING: the following line is written just so to prevent overflow
	addedNoise += to_xdouble(context.zMStar.getPhiM()) * to_xdouble(pSpace)*to_xdouble(pSpace)
	  * xexp(2*digitSize) * context.stdev*context.stdev / 4.0;

#ifdef VERBOSE
	std::cout << "\taddedNoise += to_xdouble(context.zMStar.getPhiM()) * to_xdouble(pSpace)*to_xdouble(pSpace) * xexp(2*digitSize) * context.stdev*context.stdev / 4.0;" << std::endl
			  << "\tto_xdouble(context.zMStar.getPhiM()): " << to_xdouble(context.zMStar.getPhiM()) << std::endl
			  << "\tto_xdouble(pSpace): " << to_xdouble(pSpace) << std::endl
			  << "\tdigitSize: " << digitSize << std::endl
			  << "\txexp(2*digitSize): " << xexp(2*digitSize) << std::endl
			  << "\tcontext.stdev: " << context.stdev << std::endl
			  << "\taddedNoise: " << addedNoise << std::endl;
#endif

	sizeLeft -= digitSize;
  }

  // Sanity-check: make sure that the added noise is not more than the special
  // primes can handle: After dividing the added noise by the product of all
  // the special primes, it should be smaller than the added noise term due
  // to modulus switching, i.e., keyWeight * phi(m) * pSpace^2 / 12

  long keyWeight = pubKey.getSKeyWeight(p.skHandle.getSecretKeyID());
  double phim = context.zMStar.getPhiM();
  double logModSwitchNoise = log((double)keyWeight)
	+2*log(pSpace) +log(phim) -log(12.0);
  double logKeySwitchNoise = log(addedNoise)
	-2*context.logOfProduct(context.specialPrimes);
  assert(logKeySwitchNoise < logModSwitchNoise);

  // Break the ciphertext part into digits, if needed, and scale up these
  // digits using the special primes. This is the most expensive operation
  // during homormophic evaluation, so it should be thoroughly optimized.

  vector<DoubleCRT> polyDigits;
  p.breakIntoDigits(polyDigits, nDigits);

  // Finally we multiply the vector of digits by the key-switching matrix

  // An object to hold the pseudorandom ai's, note that it must be defined
  // with the maximum number of levels, else the PRG will go out of synch.
  // FIXME: This is a bug waiting to happen.
  DoubleCRT ai(context);

  // Set the first ai using the seed, subsequent ai's (if any) will
  // use the evolving RNG state (NOTE: this is not thread-safe)

  RandomState state;
  SetSeed(W.prgSeed);

  // Add the columns in, one by one
  DoubleCRT tmp(context, IndexSet::emptySet());

  for (unsigned long i=0; i<polyDigits.size(); i++) {
	ai.randomize();
	tmp = polyDigits[i];

	// The operations below all use the IndexSet of tmp

	// add part*a[i] with a handle pointing to base of W.toKeyID
	tmp.Mul(ai,  /*matchIndexSet=*/false);
	addPart(tmp, SKHandle(1,1,W.toKeyID), /*matchPrimeSet=*/true);

	// add part*b[i] with a handle pointing to one
	polyDigits[i].Mul(W.b[i], /*matchIndexSet=*/false);
	addPart(polyDigits[i], SKHandle(), /*matchPrimeSet=*/true);
  }
  noiseVar += addedNoise;
#ifdef VERBOSE
  std::cout << "\tnoiseVar += addedNoise;" << std::endl
		    << "\taddedNoise: " << addedNoise << std::endl
		    << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif
} // restore random state upon destruction of the RandomState, see NumbTh.h
#endif

// Find the IndexSet such that modDown to that set of primes makes the
// additive term due to rounding into the dominant noise term 
void Ctxt::findBaseSet(IndexSet& s) const
{
#ifdef VERBOSE
	std::cout << "Ctxt::findBaseSet" << std::endl;
#endif
  if (getNoiseVar()<=0.0) { // an empty ciphertext
    s = context.ctxtPrimes;
    return;
  }

  assert(verifyPrimeSet());
  bool halfSize = context.containsSmallPrime();
  double curNoise = log(getNoiseVar())/2;
  double firstNoise = context.logOfPrime(0);
  double noiseThreshold = log(modSwitchAddedNoiseVar())*0.55;
#ifdef VERBOSE
std::cout << "\tcurNoise: " << curNoise << " (" << curNoise/log(2) << ")" << std::endl;
std::cout << "\tnoiseThreshold: " << noiseThreshold << curNoise << " (" << noiseThreshold/log(2) << ")" << std::endl;
#endif
  // FIXME: The above should have been 0.5. Making it a bit more means
  // that we will mod-switch a little less frequently, whether this is
  // a good thing needs to be tested.

  // remove special primes, if they are included in this->primeSet
  s = getPrimeSet();
  if (!s.disjointFrom(context.specialPrimes)) { 
    // scale down noise
    curNoise -= context.logOfProduct(context.specialPrimes);
#ifdef VERBOSE
    std::cout << "\tcurNoise -= context.logOfProduct(context.specialPrimes);" << std::endl
    		  << "\tcontext.logOfProduct(context.specialPrimes): " << context.logOfProduct(context.specialPrimes) << std::endl
    		  << "\tcurNoise: " << curNoise << std::endl;
#endif
    s.remove(context.specialPrimes);
  }

  /* We compare below to noiseThreshold+1 rather than to noiseThreshold
   * to make sure that if you mod-switch down to c.findBaseSet() and
   * then immediately call c.findBaseSet() again, it will not tell you
   * to mod-switch further down. Note that mod-switching adds close to
   * noiseThreshold to the scaled noise, so if the scaled noise was
   * equal to noiseThreshold then after mod-switchign you would have
   * roughly twice as much noise. Since we're mesuring the log, it means
   * that you may have as much as noiseThreshold+log(2), which we round
   * up to noiseThreshold+1 in the test below.
   */
  if (curNoise<=noiseThreshold+1) return; // no need to mod down

  // if the first prime in half size, begin by removing it
  if (halfSize && s.contains(0)) {
    curNoise -= firstNoise;
#ifdef VERBOSE
    std::cout << "\tcurNoise -= firstNoise;" << std::endl
    		  << "\tfirstNoise: " << firstNoise << std::endl
    		  << "\tcurNoise: " << curNoise << std::endl;
#endif
    s.remove(0);
  }

  // while noise is larger than threshold, scale down by the next prime
  while (curNoise>noiseThreshold && !empty(s)) {
    curNoise -= context.logOfPrime(s.last());
#ifdef VERBOSE
    std::cout << "\tcurNoise -= context.logOfPrime(s.last());" << std::endl
    		  << "\tcontext.logOfPrime(s.last()): " << context.logOfPrime(s.last()) << std::endl
    		  << "\tcurNoise: " << curNoise << std::endl;
#endif
    s.remove(s.last());
  }

  // Add 1st prime if s is empty or if this does not increase noise too much
  if (empty(s) || (!s.contains(0) && curNoise+firstNoise<=noiseThreshold)) {
    s.insert(0);
    curNoise += firstNoise;
#ifdef VERBOSE
    std::cout << "\tcurNoise += firstNoise;" << std::endl
    		  << "\tfirstNoise: " << firstNoise << std::endl
    		  << "\tcurNoise: " << curNoise << std::endl;
#endif
  }

  if (curNoise>noiseThreshold && log_of_ratio()>-0.5)
    cerr << "Ctxt::findBaseSet warning: already at lowest level\n";
}


/********************************************************************/
// Ciphertext arithmetic

// Add/subtract a ciphertext part to a ciphertext.
// With negative=true we subtract, otherwise we add.
void Ctxt::addPart(const DoubleCRT& part, const SKHandle& handle, 
		   bool matchPrimeSet, bool negative)
{
  assert (&part.getContext() == &context);

  // add to the the prime-set of *this, if needed (this is expensive)
  if (matchPrimeSet && !(part.getIndexSet() <= primeSet)) {
    IndexSet setDiff = part.getIndexSet() / primeSet; // set minus
    for (size_t i=0; i<parts.size(); i++) parts[i].addPrimes(setDiff);
    primeSet.insert(setDiff);
  }

  if (parts.size()==0) { // inserting 1st part 
    primeSet = part.getIndexSet();
    parts.push_back(CtxtPart(part,handle));
    if (negative) parts.back().Negate(); // not thread-safe??
  } else {               // adding to a ciphertext with existing parts
    assert(part.getIndexSet() <= primeSet);  // Sanity check

    DoubleCRT tmp(context, IndexSet::emptySet());
    const DoubleCRT* ptr = &part;

    // mod-UP the part if needed
    IndexSet s = primeSet / part.getIndexSet();
    if (!empty(s)) { // if need to mod-UP, do it on a temporary copy
      tmp = part;
      tmp.addPrimesAndScale(s);
      ptr = &tmp;
    }
    long j = getPartIndexByHandle(handle);
    if (j>=0) { // found a matching part, add them up
      if (negative) parts[j] -= *ptr;
      else          parts[j] += *ptr;
    } else {    // no mathing part found, just append this part
      parts.push_back(CtxtPart(*ptr,handle));
      if (negative) parts.back().Negate(); // not thread-safe??
    }
  }
}

// Add/subtract a ciphertext part to a ciphertext.
// With negative=true we subtract, otherwise we add.
void Ctxt::addPart(const DoubleCRT& part, const SKHandle& handle,
		   bool matchPrimeSet, bool negative, const FHESecKey& secretKey, ZZX& m)
{
  assert (&part.getContext() == &context);

  // add to the the prime-set of *this, if needed (this is expensive)
  if (matchPrimeSet && !(part.getIndexSet() <= primeSet)) {
    IndexSet setDiff = part.getIndexSet() / primeSet; // set minus
    for (size_t i=0; i<parts.size(); i++) parts[i].addPrimes(setDiff);
    primeSet.insert(setDiff);
#ifdef VERBOSE
    std::cout << " \033[31m/!\\ Noise after add to the the prime-set of *this: " << secretKey.Noise(m, *this) << "\033[0m" << std::endl;
#endif
  }

  if (parts.size()==0) { // inserting 1st part
    primeSet = part.getIndexSet();
    parts.push_back(CtxtPart(part,handle));
#ifdef VERBOSE
    std::cout << " \033[31m/!\\ Noise after push back if parts.size()==0 : " << secretKey.Noise(m, *this) << "\033[0m" << std::endl;
#endif
    if (negative) parts.back().Negate(); // not thread-safe??
  } else {               // adding to a ciphertext with existing parts
    assert(part.getIndexSet() <= primeSet);  // Sanity check

    DoubleCRT tmp(context, IndexSet::emptySet());
    const DoubleCRT* ptr = &part;
#ifdef VERBOSE
    ZZX ns;
    ptr->toPoly(ns);
    std::cout << " \033[31m/!\\ ptr size: " << MaxBits(ns) << "\033[0m" << std::endl;
#endif
    // mod-UP the part if needed
    IndexSet s = primeSet / part.getIndexSet();
    if (!empty(s)) { // if need to mod-UP, do it on a temporary copy
      tmp = part;
#ifdef VERBOSE
      tmp.toPoly(ns);
      std::cout << " \033[31m/!\\ tmp size before addPrimesAndScale: " << MaxBits(ns) << "\033[0m" << std::endl;
#endif
      tmp.addPrimesAndScale(s);
#ifdef VERBOSE
      tmp.toPoly(ns);
      std::cout << " \033[31m/!\\ tmp size after addPrimesAndScale: " << MaxBits(ns) << "\033[0m" << std::endl;
#endif
      ptr = &tmp;
    }
    long j = getPartIndexByHandle(handle);
    if (j>=0) { // found a matching part, add them up
      if (negative) parts[j] -= *ptr;
      else          {parts[j] += *ptr;
#ifdef VERBOSE
      std::cout << " \033[31m/!\\ Noise after add ptr to parts[" << j << "]: " << secretKey.Noise(m, *this) << "\033[0m" << std::endl;
#endif
      }
    } else {    // no mathing part found, just append this part
      parts.push_back(CtxtPart(*ptr,handle));
      if (negative) parts.back().Negate(); // not thread-safe??
    }
  }
}

#ifndef BIG_P
// Add a constant polynomial
void Ctxt::addConstant(const DoubleCRT& dcrt, double size)
{
  // If the size is not given, we use the default value phi(m)*(ptxtSpace/2)^2
  if (size < 0.0) {
    // WARNING: the following line is written to prevent integer overflow
    size = ((double) context.zMStar.getPhiM()) * ptxtSpace*ptxtSpace /4.0;
  }

  // Scale the constant, then add it to the part that points to one
  long f = (ptxtSpace>2)? rem(context.productOfPrimes(primeSet),ptxtSpace): 1;
  noiseVar += (size*f)*f;

  IndexSet delta = dcrt.getIndexSet() / primeSet; // set minus
  if (f==1 && empty(delta)) { // just add it
    addPart(dcrt, SKHandle(0,1,0));
    return;
  }

  // work with a local copy
  DoubleCRT tmp = dcrt;
  if (!empty(delta)) tmp.removePrimes(delta);
  if (f!=1)          tmp *= f;
  addPart(tmp, SKHandle(0,1,0));
}
#else
// Add a constant polynomial
void Ctxt::addConstant(const DoubleCRT& dcrt, double size)
{
#ifdef VERBOSE
	  std::cout << "ctxt.addConstant:" << std::endl
			    << "\tinitSize: " << size << std::endl;
#endif
  xdouble xsize = to_xdouble(size);
  // If the size is not given, we use the default value phi(m)*(ptxtSpace/2)^2
  if (size < 0.0) {
	// WARNING: the following line is written to prevent integer overflow
	xsize = to_xdouble(context.zMStar.getPhiM()) * to_xdouble(ptxtSpace) * to_xdouble(ptxtSpace) /4.0;
  }

  // Scale the constant, then add it to the part that points to one
  ZZ f;
  rem(f, context.productOfPrimes(primeSet),ptxtSpace);
#ifdef VERBOSE
  std::cout << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif
  noiseVar += (xsize*to_xdouble(f))*to_xdouble(f);
#ifdef VERBOSE
  std::cout << "\tnoiseVar += (size*to_xdouble(f))*to_xdouble(f);" << std::endl
  		    << "\trem(f, context.productOfPrimes(primeSet),ptxtSpace);" << std::endl
  		    << "\tsize: " << size << std::endl
  		    << "\tto_xdouble(f): " << to_xdouble(f) << std::endl
            << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif

  IndexSet delta = dcrt.getIndexSet() / primeSet; // set minus
  if (f==1 && empty(delta)) { // just add it
    addPart(dcrt, SKHandle(0,1,0));
    return;
  }

  // work with a local copy
  DoubleCRT tmp = dcrt;
  if (!empty(delta)) tmp.removePrimes(delta);
  if (f!=1)          tmp *= f;
  addPart(tmp, SKHandle(0,1,0));
}
#endif

void Ctxt::addConstant(const DoubleCRT& dcrt, double size, const FHESecKey& secretKey, ZZX& m)
{
#ifdef VERBOSE
	  std::cout << "ctxt.addConstant:" << std::endl
			    << "\tinitSize: " << size << std::endl;
#endif
  // If the size is not given, we use the default value phi(m)*(ptxtSpace/2)^2
  if (size < 0.0) {
	// WARNING: the following line is written to prevent integer overflow
	size = ((double) context.zMStar.getPhiM()) * to_double(ptxtSpace)*to_double(ptxtSpace) /4.0;
  }

  // Scale the constant, then add it to the part that points to one
  ZZ f;
  rem(f, context.productOfPrimes(primeSet),ptxtSpace);
#ifdef VERBOSE
  std::cout << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif
  noiseVar += (size*to_xdouble(f))*to_xdouble(f);
#ifdef VERBOSE
  std::cout << "\tnoiseVar += (size*to_xdouble(f))*to_xdouble(f);" << std::endl
  		    << "\trem(f, context.productOfPrimes(primeSet),ptxtSpace);" << std::endl
  		    << "\tsize: " << size << std::endl
  		    << "\tto_xdouble(f): " << to_xdouble(f) << std::endl
            << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif

  IndexSet delta = dcrt.getIndexSet() / primeSet; // set minus
  if (f==1 && empty(delta)) { // just add it
    addPart(dcrt, SKHandle(0,1,0));
#ifdef VERBOSE
    std::cout << " \033[31m/!\\ Noise after addPart if f==1 && empty(delta): " << secretKey.Noise(m, *this) << "\033[0m" << std::endl;
#endif
    return;
  }
ZZX zero(0);
#ifdef VERBOSE
  std::cout << " \033[31m/!\\ Noise before addPart: " << secretKey.Noise(zero, *this) << "\033[0m" << std::endl;
#endif
  // work with a local copy
  DoubleCRT tmp = dcrt;
#ifdef VERBOSE
  ZZX ns;
  tmp.toPoly(ns);
  std::cout << " \033[31m/!\\ tmp noise before removePrimes(delta): " << MaxBits(ns) << "\033[0m" << std::endl;
#endif
  if (!empty(delta)) tmp.removePrimes(delta);
#ifdef VERBOSE
  tmp.toPoly(ns);
  std::cout << " \033[31m/!\\ tmp noise before tmp *= f: " << MaxBits(ns) << "\033[0m" << std::endl;
#endif
  if (f!=1)          tmp *= f;
#ifdef VERBOSE
  tmp.toPoly(ns);
  std::cout << " \033[31m/!\\ tmp noise after tmp *= f: " << MaxBits(ns) << "\033[0m" << std::endl;
#endif

  addPart(tmp, SKHandle(0,1,0), secretKey, m);
#ifdef VERBOSE
  std::cout << " \033[31m/!\\ Noise after addPart: " << secretKey.Noise(m, *this) << "\033[0m" << std::endl;
#endif
}

void Ctxt::addConstant(const ZZX& poly, double size, const FHESecKey& secretKey, ZZX& m)
{
#ifdef VERBOSE
	  std::cout << "ctxt.addConstant:" << std::endl
			    << "\tinitSize: " << size << std::endl;
#endif
// If the size is not given, we use the default value phi(m)*(ptxtSpace/2)^2
if (size < 0.0) {
	// WARNING: the following line is written to prevent integer overflow
	size = ((double) context.zMStar.getPhiM()) * to_double(ptxtSpace)*to_double(ptxtSpace) /4.0;
}

// Scale the constant, then add it to the part that points to one
ZZ f;
rem(f, context.productOfPrimes(primeSet),ptxtSpace);
#ifdef VERBOSE
std::cout << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif
noiseVar += (size*to_xdouble(f))*to_xdouble(f);
#ifdef VERBOSE
std::cout << "\tnoiseVar += (size*to_xdouble(f))*to_xdouble(f);" << std::endl
		    << "\trem(f, context.productOfPrimes(primeSet),ptxtSpace);" << std::endl
		    << "\tsize: " << size << std::endl
		    << "\tto_xdouble(f): " << to_xdouble(f) << std::endl
          << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif

//tmp.toPoly(ns);
ZZX polyTmp = poly;

if (f!=1) {ZZ tp = poly[0]; tp *= f; tp %= ptxtSpace; polyTmp = to_ZZX(tp);}
DoubleCRT tmp = DoubleCRT(polyTmp,context,primeSet);
IndexSet delta = tmp.getIndexSet() / primeSet; // set minus
if (f==1 && empty(delta)) { // just add it
  addPart(tmp, SKHandle(0,1,0));
#ifdef VERBOSE
  std::cout << " \033[31m/!\\ Noise after addPart if f==1 && empty(delta): " << secretKey.Noise(m, *this) << "\033[0m" << std::endl;
#endif
  return;
}
if (!empty(delta)) tmp.removePrimes(delta);
#ifdef VERBOSE
ZZX ns;
tmp.toPoly(ns);
std::cout << " \033[31m/!\\ tmp noise after tmp *= f: " << MaxBits(ns) << "\033[0m" << std::endl;
#endif
addPart(tmp, SKHandle(0,1,0), secretKey, m);
#ifdef VERBOSE
std::cout << " \033[31m/!\\ Noise after addPart: " << secretKey.Noise(m, *this) << "\033[0m" << std::endl;
#endif
}

#ifndef BIG_P
// Add a constant polynomial
void Ctxt::addConstant(const ZZ& c)
{
  DoubleCRT dcrt(getContext(), getPrimeSet());
  long cc = rem(c, ptxtSpace); // reduce modulo plaintext space
  dcrt = cc;

  if (cc > ptxtSpace/2) cc -= ptxtSpace;
  double size = to_double(cc);

  addConstant(dcrt, size*size);
}
#else
// Add a constant polynomial
void Ctxt::addConstant(const ZZ& c)
{
	  DoubleCRT dcrt(getContext(), getPrimeSet());
	  ZZ cc;
	  rem(cc, c, ptxtSpace); // reduce modulo plaintext space
	  dcrt = cc;

	  if (cc > ptxtSpace/2) cc -= ptxtSpace;
	  double size = to_double(cc);

	  addConstant(dcrt, size*size);
}
#endif

void Ctxt::negate()
{
  for (size_t i=0; i<parts.size(); i++) parts[i].Negate();
}

#ifndef BIG_P
// Add/subtract another ciphertxt (depending on the negative flag)
void Ctxt::addCtxt(const Ctxt& other, bool negative)
{
  // Sanity check: same context and public key
  assert (&context==&other.context && &pubKey==&other.pubKey);

  // Special case: if *this is empty then just copy other
  if (this->isEmpty()) {
    *this = other;
    if (negative) negate();
    return;
  }

  // Sanity check: verify that the plaintext spaces are compatible
  long g = GCD(this->ptxtSpace, other.ptxtSpace);
  assert (g>1);
  this->ptxtSpace = g;

  // Match the prime-sets, mod-UP the arguments if needed
  IndexSet s = other.primeSet / primeSet; // set-minus
  if (!empty(s)) modUpToSet(s);

  const Ctxt* other_pt = &other;
  Ctxt tmp(pubKey, other.ptxtSpace); // a temporaty empty ciphertext

  s = primeSet / other.primeSet; // set-minus
  if (!empty(s)) { // need to mod-UP the other, use a temporary copy
    tmp = other;
    tmp.modUpToSet(s);
    other_pt = &tmp;
  }

  // Go over the parts of other, for each one check if
  // there is a matching part in *this
  for (size_t i=0; i<other_pt->parts.size(); i++) {
    const CtxtPart& part = other_pt->parts[i];
    long j = getPartIndexByHandle(part.skHandle);
    if (j>=0) { // found a matching part, add them up
      if (negative) parts[j] -= part;
      else          parts[j] += part;
    } else {    // no mathing part found, just append this part
      parts.push_back(part);
      if (negative) parts.back().Negate(); // not thread safe??
    }
  }
  noiseVar += other_pt->noiseVar;
}
#else
// Add/subtract another ciphertxt (depending on the negative flag)
void Ctxt::addCtxt(const Ctxt& other, bool negative)
{
  // Sanity check: same context and public key
  assert (&context==&other.context && &pubKey==&other.pubKey);

  // Special case: if *this is empty then just copy other
  if (this->isEmpty()) {
    *this = other;
    if (negative) negate();
    return;
  }

  // Sanity check: verify that the plaintext spaces are compatible
  ZZ g;
  GCD(g, this->ptxtSpace, other.ptxtSpace);
  assert (g>1);
  this->ptxtSpace = g;

  // Match the prime-sets, mod-UP the arguments if needed
  IndexSet s = other.primeSet / primeSet; // set-minus
  if (!empty(s)) modUpToSet(s);

  const Ctxt* other_pt = &other;
  Ctxt tmp(pubKey, other.ptxtSpace); // a temporaty empty ciphertext

  s = primeSet / other.primeSet; // set-minus
  if (!empty(s)) { // need to mod-UP the other, use a temporary copy
    tmp = other;
    tmp.modUpToSet(s);
    other_pt = &tmp;
  }

  // Go over the parts of other, for each one check if
  // there is a matching part in *this
  for (size_t i=0; i<other_pt->parts.size(); i++) {
    const CtxtPart& part = other_pt->parts[i];
    long j = getPartIndexByHandle(part.skHandle);
    if (j>=0) { // found a matching part, add them up
      if (negative) parts[j] -= part;
      else          parts[j] += part;
    } else {    // no mathing part found, just append this part
      parts.push_back(part);
      if (negative) parts.back().Negate(); // not thread safe??
    }
  }
  noiseVar += other_pt->noiseVar;
}
#endif

#ifndef BIG_P
// Create a tensor product of c1,c2. It is assumed that *this,c1,c2
// are defined relative to the same set of primes and plaintext space,
// and that *this DOES NOT point to the same object as c1,c2
void Ctxt::tensorProduct(const Ctxt& c1, const Ctxt& c2)
{
  // c1,c2 may be scaled, so multiply by the inverse scalar if needed
  long f = 1;
  if (c1.ptxtSpace>2) 
    f = rem(context.productOfPrimes(c1.primeSet),c1.ptxtSpace);
  if (f!=1) f = InvMod(f,c1.ptxtSpace);

  clear();                // clear *this, before we start adding things to it
  primeSet = c1.primeSet; // set the correct prime-set before we begin

  // The actual tensoring
  CtxtPart tmpPart(context, IndexSet::emptySet()); // a scratch CtxtPart
  for (size_t i=0; i<c1.parts.size(); i++) {
    CtxtPart thisPart = c1.parts[i];
    if (f!=1) thisPart *= f;
    for (size_t j=0; j<c2.parts.size(); j++) {
      tmpPart = c2.parts[j];
      // What secret key will the product point to?
      if (!tmpPart.skHandle.mul(thisPart.skHandle, tmpPart.skHandle))
	Error("Ctxt::tensorProduct: cannot multiply secret-key handles");

      tmpPart *= thisPart; // The element of the tensor product

      // Check if we already have a part relative to this secret-key handle
      long k = getPartIndexByHandle(tmpPart.skHandle);
      if (k >= 0) // found a matching part
	parts[k] += tmpPart;
      else
	parts.push_back(tmpPart);
    }
  }

  /* Compute the noise estimate as c1.noiseVar * c2.noiseVar * factor
   * where the factor depends on the handles of c1,c2. Specifically,
   * if the largest powerOfS in c1,c2 are n1,n2, respectively, then we
   * have factor = ((n1+n2) choose n2).
   */
  long n1=0,  n2=0;
  for (size_t i=0; i<c1.parts.size(); i++) // get largest powerOfS in c1
    if (c1.parts[i].skHandle.getPowerOfS() > n1)
      n1 = c1.parts[i].skHandle.getPowerOfS();
  for (size_t i=0; i<c2.parts.size(); i++) // get largest powerOfS in c2
    if (c2.parts[i].skHandle.getPowerOfS() > n2)
      n2 = c2.parts[i].skHandle.getPowerOfS();

  // compute ((n1+n2) choose n2)
  long factor = 1;
  for (long i=n1+1; i<=n1+n2; i++) factor *= i;
  for (long i=n2  ; i>1     ; i--) factor /= i;

  noiseVar = c1.noiseVar * c2.noiseVar * factor * context.zMStar.get_cM();
  if (f!=1) {
    // WARNING: the following line is written just so to prevent overflow
    noiseVar = (noiseVar*f)*f; // because every product was scaled by f
  }
}
#else
// Create a tensor product of c1,c2. It is assumed that *this,c1,c2
// are defined relative to the same set of primes and plaintext space,
// and that *this DOES NOT point to the same object as c1,c2
void Ctxt::tensorProduct(const Ctxt& c1, const Ctxt& c2)
{
#ifdef VERBOSE
   std::cout << "ctxt.tensorProduct" << std::endl;
#endif

  // c1,c2 may be scaled, so multiply by the inverse scalar if needed
  ZZ f(1);
  if (c1.ptxtSpace>2)
	rem(f, context.productOfPrimes(c1.primeSet),c1.ptxtSpace);
  if (f!=1) InvMod(f,f,c1.ptxtSpace);

  clear();                // clear *this, before we start adding things to it
  primeSet = c1.primeSet; // set the correct prime-set before we begin

  // The actual tensoring
  CtxtPart tmpPart(context, IndexSet::emptySet()); // a scratch CtxtPart
  for (size_t i=0; i<c1.parts.size(); i++) {
	CtxtPart thisPart = c1.parts[i];
	if (f!=1) thisPart *= f;
	for (size_t j=0; j<c2.parts.size(); j++) {
	  tmpPart = c2.parts[j];
	  // What secret key will the product point to?
	  if (!tmpPart.skHandle.mul(thisPart.skHandle, tmpPart.skHandle))
	Error("Ctxt::tensorProduct: cannot multiply secret-key handles");

	  tmpPart *= thisPart; // The element of the tensor product

	  // Check if we already have a part relative to this secret-key handle
	  long k = getPartIndexByHandle(tmpPart.skHandle);
	  if (k >= 0) // found a matching part
	parts[k] += tmpPart;
	  else
	parts.push_back(tmpPart);
	}
  }

  /* Compute the noise estimate as c1.noiseVar * c2.noiseVar * factor
   * where the factor depends on the handles of c1,c2. Specifically,
   * if the largest powerOfS in c1,c2 are n1,n2, respectively, then we
   * have factor = ((n1+n2) choose n2).
   */
  long n1=0,  n2=0;
  for (size_t i=0; i<c1.parts.size(); i++) // get largest powerOfS in c1
	if (c1.parts[i].skHandle.getPowerOfS() > n1)
	  n1 = c1.parts[i].skHandle.getPowerOfS();
  for (size_t i=0; i<c2.parts.size(); i++) // get largest powerOfS in c2
	if (c2.parts[i].skHandle.getPowerOfS() > n2)
	  n2 = c2.parts[i].skHandle.getPowerOfS();

  // compute ((n1+n2) choose n2)
  long factor = 1;
  for (long i=n1+1; i<=n1+n2; i++) factor *= i;
  for (long i=n2  ; i>1     ; i--) factor /= i;

  noiseVar = c1.noiseVar * c2.noiseVar * factor * context.zMStar.get_cM();
#ifdef VERBOSE
  std::cout << "\tnoiseVar = c1.noiseVar * c2.noiseVar * factor * context.zMStar.get_cM();" << std::endl
            << "\tc1.noiseVar: " << c1.noiseVar << " (" << log(c1.noiseVar)/log(2)/2 << ")" << std::endl
            << "\tc2.noiseVar: " << c2.noiseVar << " (" << log(c2.noiseVar)/log(2)/2 << ")" << std::endl
            << "\tfactor: " << factor << std::endl
            << "\tcontext.zMStar.get_cM(): " << context.zMStar.get_cM() << std::endl
            << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif

  if (f!=1) {
#ifdef VERBOSE
	  std::cout << "f!=1" << std::endl;
#endif
	// WARNING: the following line is written just so to prevent overflow
	noiseVar = (noiseVar*to_xdouble(f))*to_xdouble(f); // because every product was scaled by f
#ifdef VERBOSE
	std::cout << "\tnoiseVar = (noiseVar*to_xdouble(f))*to_xdouble(f);" << std::endl
			  << "\tto_xdouble(f): " << to_xdouble(f) << std::endl
			  << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif
  }
}
#endif

long Ctxt::uselessPrimes(FHESecKey& secretKey, ZZX& m)
{
	IndexSet s = getPrimeSet();;
  if (getNoiseVar()<=0.0) { // an empty ciphertext
    s = context.ctxtPrimes;
    return card(s);
  }
  assert(verifyPrimeSet());
  bool halfSize = context.containsSmallPrime();
  double curNoise = secretKey.Noise(m, *this);
  double firstNoise = context.logOfPrime(0);
  double noiseThreshold = 0;

  if (!s.disjointFrom(context.specialPrimes)) {
    curNoise -= context.logOfProduct(context.specialPrimes);
    s.remove(context.specialPrimes);
  }

  if (curNoise<=noiseThreshold+1) return card(s); // no need to mod down

  if (halfSize && s.contains(0)) {
    curNoise -= firstNoise;
    s.remove(0);
  }

  while (curNoise>noiseThreshold && !empty(s)) {
    curNoise -= context.logOfPrime(s.last());
    s.remove(s.last());
  }

  return card(s);
//	IndexSet s = primeSet;
//	if (!s.disjointFrom(context.specialPrimes)) {
//		modDownToSet(s / context.specialPrimes);
//	}
//	double curNoise = log(noiseVar)/2;
//	long nbPrimes = card(s);
//	while (curNoise > 0 && !empty(s))
//	{
//		curNoise -= context.logOfPrime(s.last());
//		s.remove(s.last());
//	}
//	return card(s);
}

#ifndef BIG_P
Ctxt& Ctxt::operator*=(const Ctxt& other)
{
  FHE_TIMER_START;
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return  *this;

  // Sanity check: plaintext spaces are compatible
  long g = GCD(ptxtSpace, other.ptxtSpace);
  assert (g>1);
  this->ptxtSpace = g;
  Ctxt tmpCtxt(this->pubKey, this->ptxtSpace); // a scratch ciphertext

  if (this == &other) { // a squaring operation
    modDownToLevel(findBaseLevel());      // mod-down if needed
    tmpCtxt.tensorProduct(*this, other);  // compute the actual product
    tmpCtxt.noiseVar *= 2;     // a correction factor due to dependency
  }
  else {                // standard multiplication between two ciphertexts
    // Sanity check: same context and public key
    assert (&context==&other.context && &pubKey==&other.pubKey);

    // Match the levels, mod-DOWN the arguments if needed
    long lvl = findBaseLevel();
    long otherLvl = other.findBaseLevel();
    if (lvl > otherLvl) lvl = otherLvl; // the smallest of the two

    // mod-DOWN *this, if needed (also removes special primes, if any)
    modDownToLevel(lvl);

    // mod-DOWN other, if needed
    if (primeSet!=other.primeSet){ // use temporary copy to mod-DOWN other
      Ctxt tmpCtxt1 = other;
      tmpCtxt1.modDownToLevel(lvl);
      tmpCtxt.tensorProduct(*this,tmpCtxt1); // compute the actual product
    }
    else 
      tmpCtxt.tensorProduct(*this, other);   // compute the actual product
  }
  *this = tmpCtxt; // copy the result into *this

  FHE_TIMER_STOP;
  return *this;
}
#else
Ctxt& Ctxt::operator*=(const Ctxt& other)
{
  FHE_TIMER_START;
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return  *this;

  // Sanity check: plaintext spaces are compatible
  ZZ g;
  GCD(g, ptxtSpace, other.ptxtSpace);
  assert (g>1);
  this->ptxtSpace = g;
  Ctxt tmpCtxt(this->pubKey, this->ptxtSpace); // a scratch ciphertext

  if (this == &other) { // a squaring operation
    modDownToLevel(findBaseLevel());      // mod-down if needed
    tmpCtxt.tensorProduct(*this, other);  // compute the actual product
    tmpCtxt.noiseVar *= 2;     // a correction factor due to dependency
  }
  else {                // standard multiplication between two ciphertexts
    // Sanity check: same context and public key
    assert (&context==&other.context && &pubKey==&other.pubKey);

    // Match the levels, mod-DOWN the arguments if needed
    long lvl = findBaseLevel();
    long otherLvl = other.findBaseLevel();
    if (lvl > otherLvl) lvl = otherLvl; // the smallest of the two

    // mod-DOWN *this, if needed (also removes special primes, if any)
    modDownToLevel(lvl);

    // mod-DOWN other, if needed
    if (primeSet!=other.primeSet){ // use temporary copy to mod-DOWN other
      Ctxt tmpCtxt1 = other;
      tmpCtxt1.modDownToLevel(lvl);
      tmpCtxt.tensorProduct(*this,tmpCtxt1); // compute the actual product
    }
    else
      tmpCtxt.tensorProduct(*this, other);   // compute the actual product
  }
  *this = tmpCtxt; // copy the result into *this

  FHE_TIMER_STOP;
  return *this;
}
#endif

// Higher-level multiply routines that include also modulus-switching
// and re-linearization

void Ctxt::multiplyBy(const Ctxt& other, FHESecKey secretKey, ZZX m)
{
#ifdef VERBOSE
	std::cout << "Ctxt::multiplyBy" << std::endl
			  << "\tnoiseVar: " << noiseVar << " (" << log(noiseVar)/log(2)/2 << ")" << std::endl;
#endif
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;
  FHE_TIMER_START;

  // Sanity check: plaintext spaces are compatible
  ZZ g;
  GCD(g, ptxtSpace, other.ptxtSpace);
  assert (g>1);
  this->ptxtSpace = g;
  Ctxt tmpCtxt(this->pubKey, this->ptxtSpace); // a scratch ciphertext

  if (this == &other) { // a squaring operation
    modDownToLevel(findBaseLevel());      // mod-down if needed
#ifdef VERBOSE
	  std::cout << "e mod: " << log(this->getNoiseVar())/log(2)/2 << " (" << this->getNoiseVar() << ")" << std::endl;
	  std::cout << "r mod: " << secretKey.Noise(m, *this) << " (" << xexp(2*secretKey.Noise(m, *this)*log(2)) << ")" << std::endl;
#endif
	  //this->noiseVar = xexp(2*secretKey.Noise(m, *this));
    tmpCtxt.tensorProduct(*this, other);  // compute the actual product
    tmpCtxt.noiseVar *= 2;     // a correction factor due to dependency
#ifdef VERBOSE
    std::cout << "\ttmpCtxt.noiseVar *= 2;" << std::endl
    		  << "tmpCtxt.noiseVar: " << tmpCtxt.noiseVar << " (" << log(tmpCtxt.noiseVar)/log(2)/2 << ")" << std::endl;
#endif
#ifdef VERBOSE
	  std::cout << "e *: " << log(tmpCtxt.getNoiseVar())/log(2)/2 << " (" << tmpCtxt.getNoiseVar() << ")" << std::endl;
	  std::cout << "r *: " << secretKey.Noise(m, tmpCtxt) << " (" << xexp(2*secretKey.Noise(m, tmpCtxt)*log(2)) << ")" << std::endl;
#endif
  }
  else {                // standard multiplication between two ciphertexts
    // Sanity check: same context and public key
    assert (&context==&other.context && &pubKey==&other.pubKey);

    // Match the levels, mod-DOWN the arguments if needed
    long lvl = findBaseLevel();
    long otherLvl = other.findBaseLevel();
    if (lvl > otherLvl) lvl = otherLvl; // the smallest of the two

    // mod-DOWN *this, if needed (also removes special primes, if any)
    modDownToLevel(lvl);

    // mod-DOWN other, if needed
    if (primeSet!=other.primeSet){ // use temporary copy to mod-DOWN other
      Ctxt tmpCtxt1 = other;
      tmpCtxt1.modDownToLevel(lvl);
      tmpCtxt.tensorProduct(*this,tmpCtxt1); // compute the actual product
    }
    else
      tmpCtxt.tensorProduct(*this, other);   // compute the actual product
  }
  *this = tmpCtxt; // copy the result into *this

  FHE_TIMER_STOP;

  reLinearize();   // re-linearize
#ifdef VERBOSE
  std::cout << "e relin: " << log(this->getNoiseVar())/log(2)/2 << " (" << this->getNoiseVar() << ")" << std::endl;
  std::cout << "r relin: " << secretKey.Noise(m, *this) << " (" << xexp(2*secretKey.Noise(m, *this)*log(2)) << ")" << std::endl;
#endif

//FIXME: magic number
//  noiseVar *= xexp(400);
}

void Ctxt::multiplyBy(const Ctxt& other)
{
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;
  *this *= other;  // perform the multiplication

  reLinearize();   // re-linearize

//FIXME: The must magic number is 42
//  noiseVar *= xexp(2*700);
}
void Ctxt::multiplyBy2(const Ctxt& other1, const Ctxt& other2)
{
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  long lvl = findBaseLevel();
  long lvl1 = other1.findBaseLevel();
  long lvl2 = other2.findBaseLevel();

  if (lvl<lvl1 && lvl<lvl2){ // if both others at higher levels than this,
    Ctxt tmp = other1;       // multiply others by each other, then by this
    if (&other1 == &other2) tmp *= tmp; // squaring rather than multiplication
    else                    tmp *= other2;

    *this *= tmp;
    reLinearize(); // re-linearize after all the multiplications
    return;
  }

  const Ctxt *first, *second;
  if (lvl<lvl2) { // lvl1<=lvl<lvl2, multiply by other2, then by other1
    first = &other2;
    second = &other1;
  }
  else { // multiply first by other1, then by other2
    first = &other1;
    second = &other2;
  }

  if (this == second) { // handle pointer collision
    Ctxt tmp = *second;
    *this *= *first;
    *this *= tmp;
    if (this == first) // cubing operation
      noiseVar *= 3;   // a correction factor due to dependency
    else
      noiseVar *= 2;   // a correction factor due to dependency
  } else {
    *this *= *first;
    *this *= *second;
  }
  reLinearize(); // re-linearize after all the multiplications
}

#ifndef BIG_P
// Multiply-by-constant
void Ctxt::multByConstant(const ZZ& c)
{
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;
  FHE_TIMER_START;

  long cc = rem(c, ptxtSpace); // reduce modulo plaintext space
  ZZ tmp = to_ZZ(cc);

  // multiply all the parts by this constant
  for (size_t i=0; i<parts.size(); i++) parts[i] *= tmp;

  if (cc > ptxtSpace/2) cc -= ptxtSpace;
  double size = to_double(cc);
  noiseVar *= size*size * context.zMStar.get_cM();
}
#else
// Multiply-by-constant
void Ctxt::multByConstant(const ZZ& c)
{
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;
  FHE_TIMER_START;

  ZZ tmp;
  rem(tmp, c, ptxtSpace); // reduce modulo plaintext space

  // multiply all the parts by this constant
  for (size_t i=0; i<parts.size(); i++) parts[i] *= tmp;

  if (tmp > ptxtSpace/2) tmp -= ptxtSpace;
  double size = to_double(tmp);
  noiseVar *= size*size * context.zMStar.get_cM();

//  noiseVar *= xexp(2*350);
}
#endif

#ifndef BIG_P
// Multiply-by-constant
void Ctxt::multByConstant(const DoubleCRT& dcrt, double size)
{
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;
  FHE_TIMER_START;

   // If the size is not given, we use the default value phi(m)*ptxtSpace^2/2
  if (size < 0.0) {
    // WARNING: the following line is written just so to prevent overflow
    size = ((double) context.zMStar.getPhiM()) * ptxtSpace * (ptxtSpace /4.0);
  }

  // multiply all the parts by this constant
  for (size_t i=0; i<parts.size(); i++) 
    parts[i].Mul(dcrt,/*matchIndexSets=*/false);

  noiseVar *= size * context.zMStar.get_cM();
}
#else
// Multiply-by-constant
void Ctxt::multByConstant(const DoubleCRT& dcrt, double size)
{
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;
  FHE_TIMER_START;

   // If the size is not given, we use the default value phi(m)*ptxtSpace^2/2
  if (size < 0.0) {
	// WARNING: the following line is written just so to prevent overflow
	size = ((double) context.zMStar.getPhiM()) * to_double(ptxtSpace) * (to_double(ptxtSpace) /4.0);
  }

  // multiply all the parts by this constant
  for (size_t i=0; i<parts.size(); i++)
	parts[i].Mul(dcrt,/*matchIndexSets=*/false);

  noiseVar *= size * context.zMStar.get_cM();
}
#endif

void Ctxt::multByConstant(const ZZX& poly, double size)
{
  if (this->isEmpty()) return;
  FHE_TIMER_START;
  DoubleCRT dcrt(poly,context,primeSet);
  multByConstant(dcrt,size);
}

// Divide a cipehrtext by 2. It is assumed that the ciphertext
// encrypts an even polynomial and has plaintext space 2^r for r>1.
// As a side-effect, the plaintext space is halved from 2^r to 2^{r-1}
// If these assumptions are not met then the result will not be a
// valid ciphertext anymore.
void Ctxt::divideBy2()
{
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;
  assert (ptxtSpace % 2 == 0 && ptxtSpace>2);

  // multiply all the parts by (productOfPrimes+1)/2
  ZZ twoInverse; // set to (Q+1)/2
  getContext().productOfPrimes(twoInverse, getPrimeSet());
  twoInverse += 1;
  twoInverse /= 2;
  for (size_t i=0; i<parts.size(); i++)
    parts[i] *= twoInverse;

  noiseVar /= 4;  // noise is halved by this operation
  ptxtSpace /= 2; // and so is the plaintext space
}

#ifndef BIG_P
// Divide a cipehrtext by p, for plaintext space p^r, r>1. It is assumed
// that the ciphertext encrypts a polynomial which is zero mod p. If this
// is not the case then the result will not be a valid ciphertext anymore.
// As a side-effect, the plaintext space is reduced from p^r to p^{r-1}.
void Ctxt::divideByP()
{
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  long p = getContext().zMStar.getP();
  assert (ptxtSpace>p);

  // multiply all the parts by p^{-1} mod Q (Q=productOfPrimes)
  ZZ pInverse, Q;
  getContext().productOfPrimes(Q, getPrimeSet());
  InvMod(pInverse, conv<ZZ>(p), Q);
  for (size_t i=0; i<parts.size(); i++)
    parts[i] *= pInverse;

  noiseVar /= (p * (double)p);  // noise is reduced by a p factor
  ptxtSpace /= p;               // and so is the plaintext space
}
#else
// Divide a cipehrtext by p, for plaintext space p^r, r>1. It is assumed
// that the ciphertext encrypts a polynomial which is zero mod p. If this
// is not the case then the result will not be a valid ciphertext anymore.
// As a side-effect, the plaintext space is reduced from p^r to p^{r-1}.
void Ctxt::divideByP()
{
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  ZZ p;
  p = getContext().ModulusP();
  assert (ptxtSpace>p);

  // multiply all the parts by p^{-1} mod Q (Q=productOfPrimes)
  ZZ pInverse, Q;
  getContext().productOfPrimes(Q, getPrimeSet());
  InvMod(pInverse, p, Q);
  for (size_t i=0; i<parts.size(); i++)
    parts[i] *= pInverse;

  noiseVar /= (to_xdouble(p) * to_xdouble(p));  // noise is reduced by a p factor
  ptxtSpace /= p;               // and so is the plaintext space
}
#endif

void Ctxt::automorph(long k) // Apply automorphism F(X)->F(X^k) (gcd(k,m)=1)
{
  FHE_TIMER_START;
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  // Sanity check: verify that k \in Zm*
  assert (context.zMStar.inZmStar(k));
  long m = context.zMStar.getM();

  // Apply this automorphism to all the parts
  for (size_t i=0; i<parts.size(); i++) { 
    parts[i].automorph(k);
    if (!parts[i].skHandle.isOne()) {
      parts[i].skHandle.powerOfX = MulMod(parts[i].skHandle.powerOfX,k,m);
    }
  }
  // no change in noise variance
  FHE_TIMER_STOP;
}


// Apply F(X)->F(X^k) followed by re-liearization. The automorphism is possibly
// evaluated via a sequence of steps, to ensure that we can re-linearize the
// result of every step.
void Ctxt::smartAutomorph(long k) 
{
  FHE_TIMER_START;
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  long m = context.zMStar.getM();

  k = mcMod(k, m);

  // Sanity check: verify that k \in Zm*
  assert (context.zMStar.inZmStar(k));

  long keyID=getKeyID();
  if (!inCanonicalForm(keyID)) {     // Re-linearize the input, if needed
    reLinearize(keyID);
    assert (inCanonicalForm(keyID)); // ensure that re-linearization succeeded
  }
  assert (pubKey.isReachable(k,keyID)); // reachable from 1

  while (k != 1) {
    const KeySwitch& matrix = pubKey.getNextKSWmatrix(k,keyID);
    long amt = matrix.fromKey.getPowerOfX();

    automorph(amt);
    reLinearize(keyID);

    k = MulMod(k, InvMod(amt,m), m);
  }
  FHE_TIMER_STOP;
}


#ifndef BIG_P
// applies the Frobenius automorphism p^j
void Ctxt::frobeniusAutomorph(long j)
{
  FHE_TIMER_START;
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  long m = context.zMStar.getM();
  long p = context.zMStar.getP();
  long d = context.zMStar.getOrdP();

  j = mcMod(j, d);
  long val = PowerMod(p, j, m);
  smartAutomorph(val);
  FHE_TIMER_STOP;
}
#endif


/********************************************************************/
// Utility methods

const long Ctxt::getKeyID() const
{
  for (size_t i=0; i<parts.size(); i++)
    if (!parts[i].skHandle.isOne()) return parts[i].skHandle.getSecretKeyID();

  return -1;
}

#ifndef BIG_P
// Estimates the added noise variance from mod-switching down
xdouble Ctxt::modSwitchAddedNoiseVar() const
{
  xdouble addedNoise = to_xdouble(0.0);

  // incorporate the secret keys' Hamming-weight
  for (size_t i=0; i<parts.size(); i++) {
    if (parts[i].skHandle.isOne()) {
      addedNoise += 1.0;
    }
    else {
      long keyId = parts[i].skHandle.getSecretKeyID();
      long d = parts[i].skHandle.getPowerOfS();
      xdouble h, t;
      h = pubKey.getSKeyWeight(keyId);


      // added noise is d! h^d
      t = h;
      for (long j = 2; j <= d; j++)
        t = t * h * j;

      addedNoise += t;
    }
  }

  // WARNING: the following line is written just so to prevent overflow
  addedNoise = (addedNoise * context.zMStar.getPhiM()) 
               * ptxtSpace * (ptxtSpace/ 12.0);

  return addedNoise;
}
#else
// Estimates the added noise variance from mod-switching down
xdouble Ctxt::modSwitchAddedNoiseVar() const
{
  xdouble addedNoise = to_xdouble(0.0);

  // incorporate the secret keys' Hamming-weight
  for (size_t i=0; i<parts.size(); i++) {
	if (parts[i].skHandle.isOne()) {
	  addedNoise += 1.0;
	}
	else {
	  long keyId = parts[i].skHandle.getSecretKeyID();
	  long d = parts[i].skHandle.getPowerOfS();
	  xdouble h, t;
	  h = pubKey.getSKeyWeight(keyId);


	  // added noise is d! h^d
	  t = h;
	  for (long j = 2; j <= d; j++)
		t = t * h * j;

	  addedNoise += t;
	}
  }


  // WARNING: the following line is written just so to prevent overflow
  addedNoise = (addedNoise * context.zMStar.getPhiM())
			   * to_xdouble(ptxtSpace) * (to_xdouble(ptxtSpace)/ 12.0);
  return addedNoise;
}
#endif


void Ctxt::reduce() const
{
  long n = parts.size();
  for (long i = 0; i < n; i++) parts[i].reduce();
}

istream& operator>>(istream& str, SKHandle& handle)
{
  //  cerr << "SKHandle[";
  seekPastChar(str,'['); // defined in NumbTh.cpp
  str >> handle.powerOfS;
  str >> handle.powerOfX;
  str >> handle.secretKeyID;
  seekPastChar(str,']');
  //  cerr << "]";  
  return str;
}

ostream& operator<<(ostream& str, const CtxtPart& p)
{
  return str << "[" << ((const DoubleCRT&)p) << endl 
	     << p.skHandle << "]"; 
}

istream& operator>>(istream& str, CtxtPart& p)
{
  //  cerr << "CtxtPart[";
  seekPastChar(str,'['); // defined in NumbTh.cpp
  str >> (DoubleCRT&) p;
  str >> p.skHandle;
  seekPastChar(str,']');
  //  cerr << "]";
  return str;
}

ostream& operator<<(ostream& str, const Ctxt& ctxt)
{
  str << "["<<ctxt.ptxtSpace<<" "<<ctxt.noiseVar<<" "<<ctxt.primeSet
      << " "<<ctxt.parts.size() << endl;
  for (size_t i=0; i<ctxt.parts.size(); i++)
    str << ctxt.parts[i] << endl;
  return str << "]";
}

istream& operator>>(istream& str, Ctxt& ctxt)
{
  //  cerr << "Ctxt[";
  seekPastChar(str,'['); // defined in NumbTh.cpp
  str >> ctxt.ptxtSpace >> ctxt.noiseVar >> ctxt.primeSet;
  long nParts;
  str >> nParts;
  ctxt.parts.resize(nParts, CtxtPart(ctxt.context,IndexSet::emptySet()));
  for (long i=0; i<nParts; i++) {
    str >> ctxt.parts[i];
    assert (ctxt.parts[i].getIndexSet()==ctxt.primeSet); // sanity-check
  }
  seekPastChar(str,']');
  //  cerr << "]";
  return str;
}


void CheckCtxt(const Ctxt& c, const char* label)
{
  cerr << "  "<<label << ", level=" << c.findBaseLevel() << ", log(noise/modulus)~" << c.log_of_ratio() << endl;
}

// The recursive incremental-product function that does the actual work
static void recursiveIncrementalProduct(Ctxt array[], long n)
{
  if (n <= 1) return; // nothing to do

  // split the array in two, first part is the highest power of two smaller
  // than n and second part is the rest

  long ell = NTL::NumBits(n-1); // 2^{l-1} <= n-1
  long n1 = 1UL << (ell-1);     // n/2 <= n1 = 2^l < n

  // Call the recursive procedure separately on the first and second parts
  recursiveIncrementalProduct(array, n1);
  recursiveIncrementalProduct(&array[n1], n-n1);

  // Multiply the last product in the 1st part into every product in the 2nd
  for (long i=n1; i<n; i++) array[i].multiplyBy(array[n1-1]);
}

// For i=n-1...0, set v[i]=prod_{j<=i} v[j]
// This implementation uses depth log n and (nlog n)/2 products
void incrementalProduct(vector<Ctxt>& v)
{
  long n = v.size();  // how many ciphertexts do we have
  if (n > 0) recursiveIncrementalProduct(&v[0], n); // do the actual work
}

// Compute the inner product of two vectors of ciphertexts, this routine uses
// the lower-level *= operator and does only one re-linearization at the end.
void innerProduct(Ctxt& result, const vector<Ctxt>& v1, const vector<Ctxt>& v2)
{
  long n = min(v1.size(), v2.size());
  if (n<=0) {
    result.clear();
    return;
  }
  result = v1[0]; result *= v2[0];
  for (long i=1; i<n; i++) {
    Ctxt tmp = v1[i];
    tmp *= v2[i];
    result += tmp;
  }
  result.reLinearize();
}

// Compute the inner product of a ciphertext vector and a constant vector
void innerProduct(Ctxt& result,
		  const vector<Ctxt>& v1, const vector<DoubleCRT>& v2)
{
  long n = min(v1.size(), v2.size());
  if (n<=0) {
    result.clear();
    return;
  }
  result = v1[0]; result.multByConstant(v2[0]);
  for (long i=1; i<n; i++) {
    Ctxt tmp = v1[i];
    tmp.multByConstant(v2[i]);
    result += tmp;
  }
}

void innerProduct(Ctxt& result,
		  const vector<Ctxt>& v1, const vector<ZZX>& v2)
{
  long n = min(v1.size(), v2.size());
  if (n<=0) {
    result.clear();
    return;
  }
  result = v1[0]; result.multByConstant(v2[0]);
  for (long i=1; i<n; i++) {
    Ctxt tmp = v1[i];
    tmp.multByConstant(v2[i]);
    result += tmp;
  }
}


#ifndef BIG_P
// Special-purpose modulus-switching for bootstrapping.
// Mod-switch to an externally-supplied modulus. The modulus need not be in
// the moduli-chain in the context, and does not even need to be a prime.
// The ciphertext *this is not affected, instead the result is returned in
// the zzParts vector, as a vector of ZZX'es. Returns an extimate for the
// noise variance after mod-switching.
#include "powerful.h"
double Ctxt::rawModSwitch(vector<ZZX>& zzParts, long toModulus) const
{
  // Ensure that new modulus is co-prime with plaintetx space
  assert(toModulus>1 && GCD(toModulus, getPtxtSpace())==1);
  const long p2r = getPtxtSpace();

  // Compute the ratio between the current modulus and the new one.
  // NOTE: toModulus is a long int, so a double for the logarithms and
  //       xdouble for the ratio itself is sufficient
  xdouble ratio = xexp(log((double)toModulus)
		       - context.logOfProduct(getPrimeSet()));

  // Compute also the ratio modulo ptxtSpace
  const ZZ fromModulus = context.productOfPrimes(getPrimeSet());
  long ratioModP = MulMod(toModulus % p2r, 
			  InvMod(rem(fromModulus,p2r),p2r), p2r);

  mulmod_precon_t precon = PrepMulModPrecon(ratioModP, p2r);

  // cerr << "## converting from mod-"<<context.productOfPrimes(getPrimeSet())
  //      << " to mod-"<<toModulus<<" (ratio="<<ratio
  //      << "), ptxtSpace="<<p2r<<endl;

  // Scale and round all the integers in all the parts
  zzParts.resize(parts.size());
  const PowerfulDCRT& p2d_conv = *context.rcData.p2dConv;
  for (size_t i=0; i<parts.size(); i++) {

    Vec<ZZ> powerful;
    p2d_conv.dcrtToPowerful(powerful, parts[i]); // conver to powerful rep

    for (long j=0; j<powerful.length(); j++) {
      const ZZ& coef = powerful[j];
      long c_mod_p = MulModPrecon(rem(coef,p2r), ratioModP, p2r, precon);
      xdouble xcoef = ratio*to_xdouble(coef); // the scaled coefficient

      // round xcoef to an integer which is equal to c_mod_p modulo ptxtSpace
      long rounded = conv<long>(floor(xcoef));
      long r_mod_p = rounded % p2r;
      if (r_mod_p < 0) r_mod_p += p2r; // r_mod_p in [0,p-1]

      if (r_mod_p != c_mod_p) {
        long delta = SubMod(c_mod_p, r_mod_p, p2r);
	// either add delta or subtract toModulus-delta
	rounded += delta;
	if (delta > toModulus-delta) rounded -= p2r;
      }
      // SetCoeff(zzParts[i],j,rounded);
      conv(powerful[j], rounded);  // store back in the powerful vector
    }
    //  if (deg(zzParts[i])<20) cerr << "]\n   scaled poly="<<zzParts[i]<<endl;
    p2d_conv.powerfulToZZX(zzParts[i],powerful); // conver to ZZX
  }

  // Return an estimate for the noise
  return to_double(noiseVar*ratio*ratio + modSwitchAddedNoiseVar());
}
#endif
