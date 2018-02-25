#include <Python.h>
#include <stdbool.h>

#include "gf2/gf2.h"

GF2* gf2generateKey(GF2* result, GF2 gen, GF2 mod, GF2 sk) {
	result[1] = sk;
	result[0] = gf2powmod(gen, sk, mod);
	return result;
}

GF2* gf2encrypt(GF2* msg, GF2 pk, GF2 gen, GF2 mod, GF2 esk) {
	msg[1] = gf2mulmod(msg[0], gf2powmod(pk, esk, mod), mod);
	msg[0] = gf2powmod(gen, esk, mod);
	return msg;
}

GF2* gf2decrypt(GF2* ciphertext, GF2 sk, GF2 mod, bool* err) {
	ciphertext[0] = gf2divmod(ciphertext[1], gf2powmod(ciphertext[0], sk, mod), mod, err);
	//ciphertext[0] = gf2mulmod(ciphertext[1], gf2modinv(gf2powmod(ciphertext[0], sk, mod), mod, err), mod);
	ciphertext[1] = 0;
	return ciphertext;
}

GF2 gf2commit(GF2 msg, GF2 gen1, GF2 gen2, GF2 mod, GF2 r) {
	return gf2mulmod(gf2powmod(gen1, msg, mod), gf2powmod(gen2, r, mod), mod);
}

bool gf2verify(GF2 msg, GF2 com, GF2 gen1, GF2 gen2, GF2 mod, GF2 r) {
	return com == gf2commit(msg, gen1, gen2, mod, r);
}

static PyObject* _gf2generateKey( PyObject *self, PyObject *args ) {
	GF2 gen, mod, sk;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "KKK", &gen, &mod, &sk))
		return NULL;
	
	GF2 result[2];
	gf2generateKey(result, gen, mod, sk);
	return Py_BuildValue("(KK)", result[0], result[1]);
}

static PyObject* _gf2encrypt( PyObject *self, PyObject *args ) {
	GF2 msg, pk, gen, mod, esk;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "KKKKK", &msg, &pk, &gen, &mod, &esk))
		return NULL;
	
	GF2 result[2];
	result[0] = msg;
	gf2encrypt(result, pk, gen, mod, esk);
	return Py_BuildValue("(KK)", result[0], result[1]);
}

static PyObject* _gf2decrypt( PyObject *self, PyObject *args ) {
	GF2 c1, c2, sk, mod;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "(KK)KK", &c1, &c2, &sk, &mod))
		return NULL;
	
	GF2 result[2];
	result[0] = c1;
	result[1] = c2;
	bool err = false;
	gf2decrypt(result, sk, mod, &err);
	if(err) {
		PyErr_SetString(PyExc_ValueError, "Decryption Error");
		return NULL;
	}
	return Py_BuildValue("K", result[0]);
}

static PyObject* _gf2commit( PyObject *self, PyObject *args ) {
	GF2 msg, gen1, gen2, mod, r;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "KKKKK", &msg, &gen1, &gen2, &mod, &r))
		return NULL;
	
	return Py_BuildValue("K", gf2commit(msg, gen1, gen2, mod, r));
}

static PyObject* _gf2verify( PyObject *self, PyObject *args ) {
	GF2 msg, com, gen1, gen2, mod, r;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "KKKKKK", &msg, &com, &gen1, &gen2, &mod, &r))
		return NULL;
	
	if(gf2verify(msg, com, gen1, gen2, mod, r)) Py_RETURN_TRUE;
	Py_RETURN_FALSE;
}

static PyMethodDef ElGamalGF2_funcs[] = {
	{"generateKey", _gf2generateKey, METH_VARARGS, "Generates an ElGamal key."},
	{"encrypt", _gf2encrypt, METH_VARARGS, "Encrypts a message."},
	{"decrypt", _gf2decrypt, METH_VARARGS, "Decrypts a message."},
	{"commit", _gf2commit, METH_VARARGS, "Creates a commitment."},
	{"verify", _gf2verify, METH_VARARGS, "Verifies a commitment."},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef ElGamalGF2_definition = { 
	PyModuleDef_HEAD_INIT,
	"ElGamalGF2",
	"",
	-1, 
	ElGamalGF2_funcs
};

PyMODINIT_FUNC PyInit_ElGamalGF2(void)
{
	Py_Initialize();

	return PyModule_Create(&ElGamalGF2_definition);
}
