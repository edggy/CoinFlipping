#include <Python.h>

#include "gf2/gf2.h"

static void printPoly(GF2* a, size_t len) {
	int j;
	for(j = 0; j < len; ++j) {
		printf("0x%llx, ", a[j]);
	}
}

GF2* polyCopy(GF2* p, GF2* copy, size_t len) {
	// Copies the polynomial P into memory `copy`
	size_t i;
	for(i = 0; i < len; ++i) {
		copy[i] = p[i];
	}
	return copy;
}

size_t polyDegree(GF2* p, size_t guess) {
	while(p[guess] == 0 && guess > 0) --guess;
	return guess;
}

GF2* polyAdd(GF2* p, GF2* q, size_t len, unsigned int lgsize) {
	// P = P + Q, P, Q are polynomials
	size_t i;
	for(i = 0; i < len; ++i) {
		p[i] = gf2add(p[i], q[i], lgsize);
	}
	return p;
}

GF2* polySub(GF2* p, GF2* q, size_t len, unsigned int lgsize) {
	// P = P - Q, P, Q are polynomials
	size_t i;
	for(i = 0; i < len; ++i) {
		p[i] = gf2sub(p[i], q[i], lgsize);
	}
	return p;
}

GF2* polyMulSubShift(GF2* p, GF2* q, GF2 s, size_t len, int shift, GF2 mod, unsigned int lgsize) {
	// P = P - s*Q * x^shift, P, Q are polynomials, s is a constant
	size_t i;
	for(i = 0; i < len; ++i) {
		p[i+shift] = gf2sub(p[i+shift], gf2mulmod(q[i], s, mod), lgsize);
	}
	return p;
}

GF2* polyMultiplyC(GF2* p, GF2 c, size_t len, GF2 mod) {
	// P = c*P, P is a polynomial, c is a constant

	size_t i;
	for(i = 0; i < len; ++i) {
		p[i] = gf2mulmod(p[i], c, mod);
	}
	return p;
}

GF2* polyMultiply2(GF2* p, GF2* q, size_t len, 
                                   GF2 mod, unsigned int lgsize) {
	// Calculates P = P*Q when degree of Q <= 1
	// Assumptions: P[-1] is allocated
	//              Q = q[1] x + q[0] i.e. degree 1
	// q = a + bx
	// pq = p(a + bx) = a*p + b*px

	GF2* px = p - 1;
	px[0] = gf2mulmod(q[0], p[0], mod);

	size_t i;
	for(i = 1; i <= len; ++i) {
		// px[i] = q[1]*px[i] + q[0]*p[i]
		px[i] = gf2add(gf2mulmod(q[1], px[i], mod), gf2mulmod(q[0], p[i], mod), lgsize);
	}

	return px;
}

GF2* polyMultiply(GF2* p, GF2* q, size_t lenP, size_t lenQ,
                             GF2 mod, unsigned int lgsize) {
	// Calculates P = P*Q
	// Assumptions: Memory from P to P+lenP+lenQ-2 must be allocated

	int i, k;
	for(k = lenP + lenQ - 2; k >= 0; --k) {
		if(k-(lenP-1) > 0) {
			p[k] = gf2mulmod(p[lenP-1], q[k-(lenP-1)], mod);
			i = lenP-2;
		}
		else {
			p[k] = gf2mulmod(p[k], q[0], mod);
			i = k-1;
		}
		for(; k-i < (int)lenQ && i >= 0; --i) {
			p[k] = gf2add(p[k], gf2mulmod(p[i], q[k-i], mod), lgsize);
		}
	}

	return p;
}


GF2* polyDivideC(GF2* n, GF2 d, size_t len, GF2 mod, bool* err) {
	// N = N/d, N is a polynomial, d is a constant

	size_t i;
	for(i = 0; i < len; ++i) {
		n[i] = gf2divmod(n[i], d, mod, err);
	}

	return n;
}

GF2** polyDivRem(GF2* n, GF2* d, size_t lenN, size_t lenD, GF2** res, 
			GF2 mod, unsigned int lgsize, bool* err) {
	// res = Q,R such that N = Q*D + R, N and D are polynomials
	// Assumptions: d != 0
	//              d[-
	
	printf("Entering polyDivRem\n");
	printf("n = ");
	printPoly(n, lenN);
	printf("\n");
	printf("d = ");
	printPoly(d, lenD);
	printf("\n");

	size_t i;
	size_t degN = lenN-1;
	size_t degD = lenD-1;
	GF2* q = res[0];
	GF2* r = res[1];
	polyCopy(n, r, lenN);
	for(i = 0; i <= degN; ++i) {
		q[i] = 0;
	}
	size_t degR = degN;
	size_t degS;
	
	GF2 c = d[degD];
	GF2 s;
	//GF2* s = (GF2*) malloc((lenN+lenD) * sizeof(GF2));

	while(degR >= degD) {
		/*printf("q = ");
		printPoly(q, lenN);
		printf("\n");

		printf("r = ");
		printPoly(r, degR + 1);
		printf("\n");*/

		degS = degR-degD;
		// s = lc(r) / lc(d) x^degS
		
		//printf("r[degR] = 0x%llx\n", r[degR]);
		//printf("c = 0x%llx\n", c);
		
		s = gf2divmod(r[degR], c, mod, err);
		
		//printf("s = 0x%llx\n", s);


		// q = q + s
		q[degS] = gf2add(q[degS], s, lgsize);

		// r = r - s*d
		//   = r - lc(r) / lc(d) x^degS * (d0 + d1 x + d2 x^2 + ...)
		//   = r - lc(r) / lc(d) * (d0 x^degS + d1 x^(degS+1) + d2 x^(degS+2) + ...)
		r = polyMulSubShift(r, d, s, degD+1, degS, mod, lgsize);

		degR = polyDegree(r, degR);
	}

	/*printf("q = ");
	printPoly(q, lenN);
	printf("\n");

	printf("r = ");
	printPoly(r, degR + 1);
	printf("\n");
	printf("Leving polyDivRem\n");
	*/
	return res;
}

GF2** polyExtEuc(GF2* a, GF2* b, size_t lenA, size_t lenB, size_t stop, GF2** res, 
			GF2 mod, unsigned int lgsize, bool* err) {
	printf("Entering polyExtEuc\n");
	// Returns polynomials g, u, v such that:
	// g is GCD(a, b)
	// a*u + b*v = g
	
	// arow = [r0, s0, t0]
	// brow = [r1, s1, t1]
	// res = [g,u,v]
	GF2* arowData[3];
	GF2** arow = arowData;
	GF2** brow = res;
	GF2* rrowData[2];
	GF2** rrow = rrowData;
	GF2** tmprow;

	
	size_t arowDegData[3];
	size_t browDegData[3];
	size_t* arowDeg = arowDegData;
	size_t* browDeg = browDegData;
	size_t rrowDeg[2];
	size_t* tmpDeg;

	GF2* scrap = (GF2*) malloc(6 * 2*lenA * sizeof(GF2));
	size_t scrapDeg;

	arow[0] = scrap + 1*2*lenA;
	arow[1] = scrap + 2*2*lenA;
	arow[2] = scrap + 3*2*lenA;
	rrow[0] = scrap + 4*2*lenA;
	rrow[1] = scrap + 5*2*lenA;

	printf("&brow[0] = %p\n", brow[0]);
	size_t i;
	for(i = 0; i < 6*2*lenA; ++i) {
		scrap[i] = 0;
	}

	printf("init arow and brow\n");
	polyCopy(a, arow[0], lenA);
	arow[1][0] = 1; arow[2][0] = 0;
	arowDeg[0] = lenA-1; arowDeg[1] = 0; arowDeg[2] = 0; 
	polyCopy(b, brow[0], lenB);
	brow[1][0] = 0; brow[2][0] = 1;
	browDeg[0] = lenB-1; browDeg[1] = 0; browDeg[2] = 0;

	printf("Entering Main Loop\n");
	printf("arow[0] = ");
	printPoly(arow[0], arowDeg[0] + 1);
	printf("\n");

	printf("brow[0] = ");
	printPoly(brow[0], browDeg[0] + 1);
	printf("\n");
	while(1) {
		rrow = polyDivRem(arow[0], brow[0], arowDeg[0]+1, browDeg[0]+1, rrow, mod, lgsize, err);
		
		rrowDeg[0] = polyDegree(rrow[0], arowDeg[0]);
		rrowDeg[1] = polyDegree(rrow[1], arowDeg[0]);


		printf("arow[0] / brow[0] = ");
		printPoly(rrow[0], rrowDeg[0] + 1);
		printf("\n");

		printf("arow[0] % brow[0] = ");
		printPoly(rrow[1], rrowDeg[1] + 1);
		printf("\n");
		
		for(i = 0; i < 3; ++i) {
			//printf("Copy\n");
			//for(i = 0; i < lenA; ++i) {
			//	scrap[i] = 0;
			//}
			polyCopy(rrow[0], scrap, rrowDeg[0]+1);
			scrapDeg = polyDegree(scrap, rrowDeg[0]);
			printf("scrapDeg = %d\n", scrapDeg);
			printf("scrap = ");
			printPoly(scrap, scrapDeg+1);
			printf("\n");
			printf("brow[%d] = ", i);
			printPoly(brow[i], browDeg[i] + 1);
			printf("\n");

			polyMultiply(scrap, brow[i], scrapDeg+1, browDeg[i]+1, mod, lgsize);
			printf("scrapDeg+browDeg[i]+1 = %d\n", scrapDeg+browDeg[i]+1);
			scrapDeg = polyDegree(scrap, scrapDeg+browDeg[i]);
			printf("scrapDeg = %d\n", scrapDeg);

			printf("scrap = ");
			printPoly(scrap, scrapDeg +1);
			printf("\n");

			printf("arow[%d] = ", i);
			printPoly(arow[i], arowDeg[i] + 1);
			printf("\n");
			arow[i] = polySub(arow[i], scrap, scrapDeg+1, lgsize);
			arowDeg[i] = polyDegree(arow[i], scrapDeg);
			printf("arow[%d] = ", i);
			printPoly(arow[i], arowDeg[i] + 1);
			printf("\n");
		}
		tmprow = arow;
		tmpDeg = arowDeg;
		arow = brow;
		arowDeg = browDeg;
		brow = tmprow;
		browDeg = tmpDeg;

		for(i = 0; i < 3; ++i) {
			printf("arow[%d] = ", i);
			printPoly(arow[i], arowDeg[i] + 1);
			printf("\n");
		}
		for(i = 0; i < 3; ++i) {
			printf("brow[%d] = ", i);
			printPoly(brow[i], browDeg[i] + 1);
			printf("\n");
		}

		printf("browDeg[0] = %d stop = %d\n", browDeg[0], stop);
		if(browDeg[0] < stop) break;
	}
	if(brow != res) {
		for(i = 0; i < 3; ++i) {
			polyCopy(brow[i], res[i], browDeg[i]+1);
		}
	}
	free(scrap);

	printf("Leaving polyExtEuc\n");
	return res;
}

GF2* lagrangeBasisGF2(GF2* res, size_t j, GF2* xs, size_t len, 
                                      GF2 mod, unsigned int lgsize, bool* err) {

	// res is an array such that the memory res[-len] though res[len-1] is allocated (or more)
	res[0] = 1;
	size_t curLen = 1;

	size_t i;	
	for(i = 1; i < len; ++i) {
		res[i] = 0;
	}
	GF2 xj = xs[j];
	GF2 xi;
	GF2 tmpPoly[2] = {0, 1};
	for(i = 0; i < len; ++i) {
		if(i == j) continue;
		xi = xs[i];
		// res *= x - xi
		tmpPoly[0] = xi;
		res = polyMultiply2(res, tmpPoly, curLen, mod, lgsize);
		++curLen;

		res = polyDivideC(res, gf2sub(xj, xi, lgsize), curLen, mod, err);
	}
	return res;
}

GF2* interpolateGF2(GF2* xs, GF2* ys, size_t len, GF2 mod, bool* err) {
	unsigned int lgsize = gf2bitlength(mod) - 1;

	GF2* res = (GF2*) malloc(len * sizeof(GF2));

	size_t i;	
	for(i = 0; i < len; ++i) {
		res[i] = 0;
	}

	GF2* polyData = (GF2*) malloc(2 * len * sizeof(GF2));
	GF2* tmpPoly;

	for(i = 0; i < len; ++i) {
		tmpPoly = polyData + len;
		tmpPoly = lagrangeBasisGF2(tmpPoly, i, xs, len, mod, lgsize, err);

		// res += y_i*l_i(x)
		res = polyAdd(res, polyMultiplyC(tmpPoly, ys[i], len, mod), len, lgsize);
	}
	free(polyData);
	return res;
}

// TODO: Bulk interpolation when the x values are the same.

GF2* decodeReedSolomonGF2(GF2* xs, GF2* ys, size_t len, size_t k, GF2 mod, bool* err) {
	printf("Entering decodeReedSolomonGF2\n");
	unsigned int lgsize = gf2bitlength(mod) - 1;

	printf("Interpolating\n");
	GF2* g1 = interpolateGF2(xs, ys, len, mod, err);
	printf("Done Interpolating\n");
	size_t g1Len = polyDegree(g1, len-1) + 1;

	printf("n = %d, k = %d, g1Len = %d, mod = %llu\n", len, k, g1Len, mod);

	if(g1Len-1 <= k) {
		return g1;
	}

	size_t i;
	GF2* memory = (GF2*) malloc(6*len * sizeof(GF2));
	for(i = 0; i < 6*len; ++i) {
		memory[i] = 0;
	}
	printf("memory = %p to %p\n", memory, memory+6*len);
	//GF2* gStart[6];
	//for(i = 0; i < 6; ++i) {
	//	gStart[i] = (GF2*) malloc(len * sizeof(GF2));
	//}

	GF2* g0 = memory + len - 1;
	g0[0] = 1;
	size_t g0Len = 1;

		
	GF2 tmpPoly[2] = {0, 1};
	for(i = 0; i < len; ++i) {
		tmpPoly[0] = xs[i];
		g0 = polyMultiply2(g0, tmpPoly, g0Len, mod, lgsize);
		++g0Len;
	}

	printf("g0 = ");
	printPoly(g0, g0Len);
	printf("\n");

	printf("g1 = ");
	printPoly(g1, g1Len);
	printf("\n");

	GF2* res[3];
	res[0] = memory + len;
	res[1] = memory + 2*len;
	res[2] = memory + 3*len;
	printf("res[0] = %p\n", res[0]);
	//printf("gStart[0] = %p\n", gStart[0]);
	
	polyExtEuc(g0, g1, g0Len, g1Len, (len+k)/2, &res, mod, lgsize, err);

	printf("res[0] = %p\n", res[0]);
	//free(gStart[0]);
	//gStart[0] = NULL;
	//g0 = NULL;

	GF2* g = res[0];
	GF2* u = res[1];
	GF2* v = res[2];

	size_t degG = polyDegree(g, (len+k)/2);
	size_t degU = polyDegree(u, (len+k)/2);
	size_t degV = polyDegree(v, (len+k)/2);

	printf("g = ");
	printPoly(g, degG + 1);
	printf("\n");
	printf("u = ");
	printPoly(u, degU + 1);
	printf("\n");
	printf("v = ");
	printPoly(v, degV + 1);
	printf("\n");

	bool cont = false;
	for(i = 0; i <= degV; ++i) {
		if(v[i] != 0) {
			cont = true;
		}
	}
	if(!cont) {
		printf("Leaving decodeReedSolomonGF2: No errors\n");
		printf("memory = %p\n", memory);
		printf("memory = ");
		printPoly(memory, 6*len);
		printf("\n");
		//for(i = 0; i < 6*len; ++i) {
		//	printf("Erasing memory[%d]\n", i);
		//	memory[i] = 0;
		//}
		printf("Memory erased");
		//gStart[0] = NULL;gStart[1] = NULL;gStart[2] = NULL;gStart[3] = NULL;
		//g0 = NULL;
		free(memory);
		printf("Memory freed\n");
		return g1;
	}
	
	GF2* f[2];
	f[0] = memory + 4*len;
	f[1] = memory + 5*len;
	polyDivRem(g, v, degG+1, degV+1, &f, mod, lgsize, err);
	size_t degF0 = polyDegree(f[0], degG);
	size_t degF1 = polyDegree(f[1], degG);
	printf("f[0] = ");
	printPoly(f[0], degF0 + 1);
	printf("\n");
	printf("f[1] = ");
	printPoly(f[1], degF1 + 1);
	printf("\n");
	
	polyCopy(f[0], g1, degF0+1);
	free(memory);
	printf("Leaving decodeReedSolomonGF2\n");
	return g1;
	
}

//-------------------------------------------------------------

typedef struct {
	size_t degree;
	size_t len;
	size_t start;
	GF2* poly;
} Poly;

typedef struct {
	GF2 mod;
	unsigned int lgsize;
} GFdata;

Poly* _polyCopy(Poly* p, Poly* copy) {
	// Copies the polynomial P into memory `copy`
	// Assumptions: copy has enough room

	size_t i;
	for(i = 0; i < p->degree+1; ++i) {
		copy->poly[i] = p->poly[i];
	}
	copy->degree = p->degree;
	return copy;
}

Poly* _polyAdd(Poly* p, Poly* q, GFdata* gf) {
	// P = P + Q, P, Q are polynomials
	size_t i;
	size_t lenMax = p->degree+1;
	size_t lenMin = q->degree+1;
	Poly* longer = p;
	if(lenMin > lenMax) {
		lenMax = q->degree+1;
		lenMin = p->degree+1;
		longer = q;
	}
	for(i = 0; i < lenMin; ++i) {
		p->poly[i] = gf2add(p->poly[i], q->poly[i], gf->lgsize);
	}
	for(i = lenMin; i < lenMax; ++i) {
		p->poly[i] = longer->poly[i];
	}
	p->degree = lenMax;
	return p;
}

Poly* _polySub(Poly* p, Poly* q, GFdata* gf) {
	// P = P - Q, P, Q are polynomials
	size_t i;
	size_t lenMax = p->degree+1;
	size_t lenMin = q->degree+1;
	Poly* longer = p;
	if(lenMin > lenMax) {
		lenMax = q->degree+1;
		lenMin = p->degree+1;
		longer = q;
	}
	for(i = 0; i < lenMin; ++i) {
		p->poly[i] = gf2sub(p->poly[i], q->poly[i], gf->lgsize);
	}
	for(i = lenMin; i < lenMax; ++i) {
		p->poly[i] = longer->poly[i];
	}
	p->degree = lenMax;
	return p;
}


Poly* _polyMultiplyC(Poly* p, GF2 c, GFdata* gf) {
	// P = c*P, P is a polynomial, c is a constant

	size_t i;
	for(i = 0; i < p->degree+1; ++i) {
		p->poly[i] = gf2mulmod(p->poly[i], c, gf->mod);
	}
	return p;
}

Poly* _polyMultiply2(Poly* p, Poly* q, GFdata* gf) {
	// Calculates P = P*Q when degree of Q <= 1
	// Assumptions: P[-1] is allocated
	//              Q = q[1] x + q[0] i.e. degree 1
	// q = a + bx
	// pq = p(a + bx) = a*p + b*px
	// p[-1] = a*p[0], p[0] = a*p[1] + b*p[0]
	// p[i] = a*p[i+1] + b*p[i]

	Poly* px = p;
	px->poly -= 1;

	px->poly[0] = gf2mulmod(q->poly[0], px->poly[1], gf->mod);

	size_t i;
	for(i = 1; i <= px->degree; ++i) {
		// px[i] = q[1]*px[i] + q[0]*p[i]
		px->poly[i] = gf2add(gf2mulmod(q->poly[1], px->poly[i], gf->mod), gf2mulmod(q->poly[0], px->poly[i+1], gf->mod), gf->lgsize);
	}
	px->degree = px->degree+1;
	px->start = px->start-1;
	return px;
}

Poly* _polyMultiply(Poly* p, Poly* q, GFdata* gf) {
	// Calculates P = P*Q
	// Assumptions: Memory from P to P+p->degree+q->degree must be allocated

	int i, k;
	for(k = p->degree + q->degree; k >= 0; --k) {
		if(k-(p->degree) > 0) {
			p->poly[k] = gf2mulmod(p->poly[p->degree], q->poly[k-(p->degree)], gf->mod);
			i = p->degree - 1;
		}
		else {
			p->poly[k] = gf2mulmod(p->poly[k], q->poly[0], gf->mod);
			i = k-1;
		}
		for(; k-i <= (int)q->degree && i >= 0; --i) {
			p->poly[k] = gf2add(p->poly[k], gf2mulmod(p->poly[i], q->poly[k-i], gf->mod), gf->lgsize);
		}
	}
	p->degree = p->degree + q->degree;
	return p;
}


Poly* _polyDivideC(Poly* n, GF2 d, GFdata* gf, bool* err) {
	// N = N/d, N is a polynomial, d is a constant

	size_t i;
	for(i = 0; i < n->degree+1; ++i) {
		n->poly[i] = gf2divmod(n->poly[i], d, gf->mod, err);
	}

	return n;
}

Poly* _lagrangeBasisGF2(Poly* res, size_t j, GF2* xs, size_t len, 
                                      GFdata* gf, bool* err) {

	// res is an array such that the memory res[-len] though res[len-1] is allocated (or more)
	res->poly[0] = 1;
	res->len = 1;

	size_t i;	
	for(i = 1; i < len; ++i) {
		res->poly[i] = 0;
	}
	GF2 xj = xs[j];
	GF2 xi;

	GF2 tmp[2] = {0, 1};
	Poly tmpPoly;
	tmpPoly.poly = tmp;
	tmpPoly.len = 2;
	tmpPoly.degree = 1;
	tmpPoly.start = 0;

	for(i = 0; i < len; ++i) {
		if(i == j) continue;
		xi = xs[i];
		// res *= x - xi
		tmpPoly.poly[0] = xi;
		res = _polyMultiply2(res, &tmpPoly, gf);

		// res /= xj - xi
		res = _polyDivideC(res, gf2sub(xj, xi, gf->lgsize), gf, err);
	}
	return res;
}

Poly* interpolateGF2Poly(GF2* xs, GF2* ys, size_t len, GF2 mod, bool* err) {
	GFdata gf;
	gf.mod = mod;
	gf.lgsize = gf2bitlength(mod) - 1;

	Poly* res;
	
	if(((res = (Poly*) malloc(sizeof(Poly))) == NULL) || ((res->poly = (GF2*) malloc(len * sizeof(GF2))) == NULL)) {
		*err = true;
		return NULL;
	}
	
	res->len = len;
	res->start = 0;
	res->degree = 0;

	size_t i;	
	for(i = 0; i < res->len; ++i) {
		res->poly[i] = (GF2)0;
	}

	Poly polyData;
	if((polyData.poly = (GF2*) malloc(2 * len * sizeof(GF2))) == NULL) {
		*err = true;
		return NULL;
	}
	polyData.len = 2*len;
	polyData.start = 0;
	polyData.degree = 0;

	Poly tmpPoly;

	for(i = 0; i < len; ++i) {
		tmpPoly.poly = polyData.poly + len;
		tmpPoly.len = 2*len;
		tmpPoly.start = len;
		tmpPoly.degree = 0;

		_lagrangeBasisGF2(&tmpPoly, i, xs, len, &gf, err);

		// res += y_i*l_i(x)
		res = _polyAdd(res, _polyMultiplyC(&tmpPoly, ys[i], &gf), &gf);
	}
	free(polyData.poly);
	return res;
	
}

static PyObject* interpolatePolynomial( PyObject *self, PyObject *args ) {
	GF2 mod;
	size_t len;

	PyObject* pList;
	PyObject* pTuple;
	PyObject* pItem;
	size_t i;

	if (!PyArg_ParseTuple(args, "O!K", &PyList_Type, &pList, &mod)) {
		PyErr_SetString(PyExc_TypeError, "parameter must be a list and an int.");
		return NULL;
	}

	len = (size_t) PyList_Size(pList);
	GF2* xs = (GF2*) malloc(len * sizeof(GF2));
	GF2* ys = (GF2*) malloc(len * sizeof(GF2));
	for (i = 0; i < len; ++i) {
		pTuple = PyList_GetItem(pList, i);
		if(!PyTuple_Check(pTuple)) {
			PyErr_SetString(PyExc_TypeError, "List must contain tuples");
			return NULL;
		}
		pItem = PyTuple_GET_ITEM(pTuple, 0);
		xs[i] = PyLong_AsUnsignedLongLongMask(pItem);

		pItem = PyTuple_GET_ITEM(pTuple, 1);
		ys[i] = PyLong_AsUnsignedLongLongMask(pItem);
	}
	
	bool err = false;
	GF2* res = interpolateGF2(xs, ys, len, mod, &err);
	free(xs);
	free(ys);
	
	if(err) {
		PyErr_SetString(PyExc_ValueError, "Interpolation Error");
		free(res);
		return NULL;
	}

	PyObject* resList = PyList_New((Py_ssize_t) len);
	PyObject* pyLong;
	for (i = 0; i < len; ++i) {
		pyLong = PyLong_FromUnsignedLongLong(res[i]);
		PyList_SET_ITEM(resList, i, pyLong);
	}
	free(res);
	return resList;
}

static PyObject* decodeReedSolomon( PyObject *self, PyObject *args ) {
	GF2 mod;
	size_t len, k;

	PyObject* pList;
	PyObject* pTuple;
	PyObject* pItem;
	size_t i;

	if (!PyArg_ParseTuple(args, "O!nK", &PyList_Type, &pList, &k, &mod)) {
		PyErr_SetString(PyExc_TypeError, "parameter must be a list and two ints.");
		return NULL;
	}

	len = (size_t) PyList_Size(pList);
	GF2* xs = (GF2*) malloc(len * sizeof(GF2));
	GF2* ys = (GF2*) malloc(len * sizeof(GF2));
	for (i = 0; i < len; ++i) {
		pTuple = PyList_GetItem(pList, i);
		if(!PyTuple_Check(pTuple)) {
			PyErr_SetString(PyExc_TypeError, "List must contain tuples");
			return NULL;
		}
		pItem = PyTuple_GET_ITEM(pTuple, 0);
		xs[i] = PyLong_AsUnsignedLongLongMask(pItem);

		pItem = PyTuple_GET_ITEM(pTuple, 1);
		ys[i] = PyLong_AsUnsignedLongLongMask(pItem);
	}
	
	bool err = false;
	printf("Decoding\n");
	GF2* res = decodeReedSolomonGF2(xs, ys, len, k, mod, &err);
	printf("Done Decoding\n");
	free(xs);
	free(ys);
	
	if(err) {
		PyErr_SetString(PyExc_ValueError, "Interpolation Error");
		free(res);
		return NULL;
	}

	PyObject* resList = PyList_New((Py_ssize_t) len);
	PyObject* pyLong;
	for (i = 0; i < len; ++i) {
		pyLong = PyLong_FromUnsignedLongLong(res[i]);
		PyList_SET_ITEM(resList, i, pyLong);
	}
	free(res);
	return resList;
}

static PyMethodDef interpolateGF2_funcs[] = {
	{"interpolatePolynomial", interpolatePolynomial, METH_VARARGS, "Interpolates a polynomial."},
	{"decodeReedSolomon", decodeReedSolomon, METH_VARARGS, "Decodes and corrects a Reed Solomon encoding."},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef interpolateGF2_definition = { 
	PyModuleDef_HEAD_INIT,
	"interpolatePolynomial",
	"",
	-1, 
	interpolateGF2_funcs
};

PyMODINIT_FUNC PyInit_interpolateGF2(void)
{
	Py_Initialize();

	return PyModule_Create(&interpolateGF2_definition);
}





