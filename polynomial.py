from functools import reduce
from collections import defaultdict
import math

class InverseException(Exception):
    pass   

class PolynomialMeta(type):
    cache = {}
        
    def __call__(cls, *args, **kwargs):
        ret = super(PolynomialMeta, cls).__call__(*args, **kwargs)
        
        if ret in cls.cache:
            return cls.cache[ret]
        else:
            cls.cache[ret] = ret
            
        return ret   

#class Polynomial(metaclass=PolynomialMeta):
class Polynomial:
    def __init__(self, *, degree = 0, coefficients = None, mod = None):
        self.mod = mod
        if coefficients is None:
            self.coefficients =  (0,) * degree + (1,)
        else:
            self.coefficients = tuple(coefficients)
            
        if self.mod is not None:
            self.coefficients = tuple(map(lambda x: x % self.mod, self.coefficients))
            
        while self.coefficients[-1] == 0 and len(self.coefficients) > 1:
            self.coefficients = self.coefficients[:-1]        
        
    def __str__(self):
        formatStr = '{sign}{coefficient}{times}{var}{exp}{power}'
        
        fstr = ''
        for i in range(len(self.coefficients)-1, -1, -1):
            fstr += formatStr.format(sign = '{sign[%s]}' % i, 
                                     coefficient = '{coefficient[%s]}' % i,
                                     times = '{times[%s]}' % i,
                                     var = '{var[%s]}' % i,
                                     exp = '{exp[%s]}' % i,
                                     power = '{power[%s]}' % i)   
        
        sign = list(map(lambda x: ' + ' if x > 0 else (' - ' if x < 0 else ''), self.coefficients))
        
        rightmost = len(self.coefficients) - 1
        while rightmost >= 0 and self.coefficients[rightmost] == 0:
            rightmost -= 1
    
        sign[rightmost] = '' if self.coefficients[rightmost] >= 0 else '-'
        
        coefficients = list(map(lambda x: repr(abs(x[1])) if (x[1] not in (0, 1, -1) or (x[0] == 0 and x[1] != 0)) else '', enumerate(self.coefficients)))
        times = list(map(lambda x: ' * ' if (x[1] not in (0, 1, -1) and x[0] != 0) else '', enumerate(self.coefficients)))
        var = list(map(lambda x: 'x' if (x[1] != 0 and x[0] != 0) else '', enumerate(self.coefficients)))
        exp = list(map(lambda x: '**' if (x[1] != 0 and x[0] not in (1, 0)) else '', enumerate(self.coefficients)))
        power = list(map(lambda x: x[0] if (x[1] != 0 and x[0] not in (1, 0)) else '', enumerate(self.coefficients)))
        
        mod = '' if self.mod is None else ' mod %d' % self.mod
        
        result = fstr.format(sign = sign, 
                           coefficient = coefficients,
                           times = times,
                           var = var,
                           exp = exp,
                           power = power)
        
        if result == '':
            result = str(self.coefficients[0])
        return result + mod
    def __repr__(self):
        return str(self)
    
    def __call__(self, x):
        val = self.coefficients[-1]
        for i in range(-2, -len(self.coefficients)-1, -1):
            val = self.coefficients[i] + x*val
            if self.mod is not None:
                val %= self.mod
        return val
        
    def __hash__(self):
        return hash((self.coefficients, self.mod))
    
    def __getitem__(self, index):
        return self.coefficients[index]
    
    def __eq__(self, other):
        if isinstance(other, Polynomial):
            return self.coefficients == other.coefficients and self.mod == other.mod
        return False
    
    def __op__(self, other, op):
        maxLen = max(len(self.coefficients), len(other.coefficients))
        coefficients = [0] * maxLen
        for i in range(maxLen):
            try:
                coefficients[i] = op(self.coefficients[i], other.coefficients[i])
            except IndexError:
                try:
                    coefficients[i] = op(self.coefficients[i], 0)
                except IndexError:
                    coefficients[i] = op(0, other.coefficients[i])
                
            
        return Polynomial(coefficients = coefficients, mod=self.mod)
    
    def __add__(self, other):
        return self.__op__(other, lambda x, y: x+y)
    
    def __radd__(self, other):
        return self.__op__(other, lambda x, y: y+x)    
    
    def __sub__(self, other):
        return self.__op__(other, lambda x, y: x-y)   
    
    def __rsub__(self, other):
        return self.__op__(other, lambda x, y: y-x)      
    
    def __mul__(self, other):
        # (a x^2 + b x + c)(d x^2 + e x + f) = ad x^4 + (ae + db) x^3 + (af + dc + eb) x^2 + (bf + ec) x + cf
        # (c, b, a) * (f, e ,d) = (a*d, a*e + b*d, a*f + b*e + c*d, b*f + c*e, c*f)
        
        try:
            other = Polynomial(coefficients=[int(other)])
        except TypeError:
            other = Polynomial(coefficients=other)
        
        coefficients = [0] * ((len(self.coefficients)-1) + (len(other.coefficients)-1) + 1)
        for i, ci in enumerate(self.coefficients):
            for j, cj in enumerate(other.coefficients):
                coefficients[i+j] += ci*cj
                if self.mod is not None:
                    coefficients[i+j] %= self.mod
        
        return Polynomial(coefficients = coefficients, mod=self.mod)
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __imul__(self, other):
        return self.__mul__(other)
    
    def __truediv__(self, other):
        # TODO: General polynomial division
        if self.mod is None:
            return Polynomial(coefficients = map(lambda x: x / other, self.coefficients), mod=self.mod)  
        else:
            # TODO: account for division not being defined
            return Polynomial(coefficients = map(lambda x: x / other, self.coefficients), mod=self.mod) 
            #return Polynomial(coefficients = map(lambda x: moddiv(x, other, self.mod), self.coefficients), mod=self.mod)  
    
    def __itruediv__(self, other):
        return self.__truediv__(other)
    
    def __floordiv__(self, other):
        if self.mod is None:
            return Polynomial(coefficients = map(lambda x: x // other, self.coefficients), mod=self.mod)
        else:
            return Polynomial(coefficients = map(lambda x: moddiv(x, other, self.mod), self.coefficients), mod=self.mod) 
        
    def __ifloordiv__(self, other):
        return self // other    
    
    def __mod__(self, other):
        if isinstance(other, int):
            return Polynomial(coefficients = self.coefficients, mod=other)
        else:
            return Polynomial(coefficients = map(lambda x: x % other, self.coefficients), mod=self.mod)
        
    def degree(self):
        rightmost = len(self.coefficients) - 1
        while rightmost >= 0 and self.coefficients[rightmost] == 0:
            rightmost -= 1
        return rightmost
    
    def extract(self, *args):
        data = b''
        for i in range(*args):
            output = self(i)
            data += output.to_bytes(t.bit_length()//8+1, 'big')
        return data

def lagrangeBasisPolynomial(j, points, mod = None):
    numerator = Polynomial(degree=0, mod=mod)
    denominator = 1
    xj, yj = points[j]
    for i, p in enumerate(points):
        if i != j:
            xi, yi = p
            numerator *= Polynomial(coefficients=[-xi, 1])
            denominator *= xj - xi
            #assert int((numerator/denominator)(xj)) == 1
            #assert int((numerator/denominator)(xi)) == 0
            #print([int(i) for i in (numerator/ denominator)])
    return [numerator, denominator]


def interpolatePolynomial(points, mod = None):
    res = Polynomial(coefficients=[0], mod=mod)
    numDens = []
    for j, p in enumerate(points):
        xj, yj = p
        #numerator, denominator = lagrangeBasisPolynomial(j, points, mod=mod)
        #res += yj * numerator / denominator
        numDens.append(lagrangeBasisPolynomial(j, points, mod=mod))
        numDens[-1][0] *= yj
        if mod is not None:
            numDens[-1][1] %= mod
        try:
            numDens[-1][0] /= numDens[-1][1]
            numDens[-1][1] = 1
        except InverseException:
            pass
        
    #print([i for i in map(lambda x: [int(i) for i in (x[0]/x[1])], numDens)])    
    if mod is None:
        for i, numDeni in enumerate(numDens):
            numerator, denominator = numDeni
            res += numerator / denominator
        return (res, 0)
    
    newDenominator = 1
    for i, numDeni in enumerate(numDens):
        numerator, denominator = numDeni
        if denominator == 1:
            continue
        newDenominator *= denominator
        for j, numDenj in enumerate(numDens):
            if i != j:
                numDens[j][0] *= denominator
                
    for i, numDeni in enumerate(numDens):
        numerator, denominator = numDeni
        res += numerator
        
    try:
        res = (res / (newDenominator % mod), 0)
        for xj, yj in points:
            assert res[0](xj) == yj
    except InverseException:
        res = (res, (newDenominator % mod))
    
    return res


if __name__ == '__main__':
   
    import random
    p = Polynomial(degree = 5)
    q = Polynomial(coefficients = [1,2,3], mod=10)
    r = Polynomial(coefficients = [90,37,63])
    
    x = 3
    print('%s @ %d = %d' % (p, x, p(x)))
    print('%s @ %d = %d' % (q, x, q(x)))
    print('%s @ %d = %d' % (r, x, r(x)))
    print('%s @ %d = %d' % (r+p, x, (r+p)(x)))
    print('%s @ %d = %d' % (r*p, x, (r*p)(x)))
    
    intPol, remainder = interpolatePolynomial([(1,r(1)),(2,r(2)),(3,r(3)),(4,r(4))])
    print('%s @ %d = %d' % (r, x, r(x)))
    print('%s @ %d = %d' % (intPol, x, intPol(x)))
    
    #inverse = modinv(Polynomial(degree = 2, mod=2), Polynomial(degree = 6, mod=2))
    
    '''print(p, p(5), eval(str(p).replace('x', '5')))
    print(q, q(5), eval(str(q).replace('x', '5')))
    print(r, r(5), eval(str(r).replace('x', '5')))'''
    
    '''gf2a = GF2(value=0b10100110)
    gf2b = GF2(value=0b10000111)
    gf2c = GF2(value=0b00000101)
    print(gf2a*gf2b)
    print(gf2a*gf2b/gf2b)
    print(~gf2c)
    print(gf2c + -gf2c)
    
    gfpoly = Polynomial(coefficients = [GF2(value=0b10, size=8), GF2(value=0b110, size=8), GF2(value=0b1110, size=8)], mod=256)
    inputValue = 0b1110011
    print('%s @ %s = %r' % (gfpoly, bin(inputValue), gfpoly(GF2(value=inputValue, size=8))))
    
    gfpoly = Polynomial(coefficients = [GF2(value=0b10, size=16), GF2(value=0b110, size=16), GF2(value=0b1110, size=16)], mod=256**2)
    print('%s @ %s = %r' % (gfpoly, bin(inputValue), gfpoly(GF2(value=inputValue, size=16))))'''