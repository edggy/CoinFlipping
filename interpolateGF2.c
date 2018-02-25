#include <Python.h>

#include "gf2/gf2.h"

/*static void printArray(GF2* a, size_t len) {
	int j;
	for(j = 0; j < len; ++j) {
		printf("%llu, ", a[j]);
	}
}*/

static GF2* polyAdd(GF2* p, GF2* q, size_t len, unsigned int lgsize) {
	// P = P + Q, P, Q are polynomials
	size_t i;
	for(i = 0; i < len; ++i) {
		p[i] = gf2add(p[i], q[i], lgsize);
	}
	return p;
}

static GF2* polyMultiply(GF2* p, GF2 c, size_t len, GF2 mod) {
	// N = c*N, N is a polynomial, c is a constant

	size_t i;
	for(i = 0; i < len; ++i) {
		p[i] = gf2mulmod(p[i], c, mod);
	}
	return p;
}

static GF2* polyMultiply2(GF2* p, GF2* q, size_t len, 
                                   GF2 mod, GF2 lgsize) {
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

static GF2* polyDivide(GF2* n, GF2 d, size_t len, GF2 mod, bool* err) {
	// N = N/q, N is a polynomial, q is a constant

	size_t i;
	for(i = 0; i < len; ++i) {
		n[i] = gf2divmod(n[i], d, mod, err);
	}

	return n;
}

static GF2* lagrangeBasisGF2(GF2* res, size_t j, GF2* xs, size_t len, 
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

		res = polyDivide(res, gf2sub(xj, xi, lgsize), curLen, mod, err);
	}
	return res;
}

static GF2* interpolateGF2(GF2* xs, GF2* ys, size_t len, GF2 mod, bool* err) {
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
		res = polyAdd(res, polyMultiply(tmpPoly, ys[i], len, mod), len, lgsize);
	}
	free(polyData);
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

static PyMethodDef interpolateGF2_funcs[] = {
	{"interpolatePolynomial", interpolatePolynomial, METH_VARARGS, "Interpolates a polynomial."},
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





