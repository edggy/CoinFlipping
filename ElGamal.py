try:
    from gf2c import findRandomIrreduciblePolynomial
    from gf2c import findRandomGeneratorPolynomial
except ImportError:
    from gf2.gf2 import findRandomIrreduciblePolynomial
    from gf2.gf2 import findRandomGeneratorPolynomial    

try:
    import ElGamalGF2
except ImportError:
    pass

def generateKey(generator, groupSize, random):
    secretKey = random.randrange(groupSize)
    try:
        return ElGamalGF2.generateKey(int(generator), generator.mod, secretKey)
    except NameError:
        publicKey = generator**secretKey
        return (publicKey, secretKey)

def encrypt(message, generator, groupSize, publicKey, random):
    ephemeralSecretKey = random.randrange(groupSize)
    try:
        return ElGamalGF2.encrypt(int(message), int(publicKey), int(generator), generator.mod, ephemeralSecretKey)
    except NameError:    
        ephemeralPublicKey = generator**ephemeralSecretKey
        sharedSecret = publicKey**ephemeralSecretKey
        c2 = int(message) * sharedSecret
        return (ephemeralPublicKey, c2)

def decrypt(ciphertext, secretKey, *, modulus = None):
    try:
        if modulus is None:
            modulus = secretKey.mod
        return ElGamalGF2.decrypt((int(ciphertext[0]), int(ciphertext[1])), int(secretKey), modulus)
    except NameError:    
        ephemeralPublicKey, c2 = ciphertext
        sharedSecret = ephemeralPublicKey**secretKey
        message = c2 / sharedSecret
        return message

class DecryptionError(Exception):
    pass

class ElGamal:
    def __init__(self, *, lgGroupSize = None , generator = None, random = None, secretKey = None, publicKey = None, newElement = None):
        self.lgGroupSize = lgGroupSize
        if self.lgGroupSize is None:
            self.lgGroupSize = 8
        
        self.random = random
        if self.random is None:
            import random
            self.random = random.SystemRandom()
            
        self.generator = generator
        if self.generator is None:
            self.mod = findRandomIrreduciblePolynomial(self.lgGroupSize, self.random)
            self.generator = findRandomGeneratorPolynomial(self.lgGroupSize, self.mod, self.random)   
        else:
            self.mod = self.generator.mod
        
        if secretKey is None and publicKey is None:
            self.publicKey, self.secretKey = generateKey(self.generator, 2**self.lgGroupSize, self.random)
            
        elif secretKey is not None:
            self.publicKey, self.secretKey = generator**secretKey, secretKey
            
        elif publicKey is not None:
            self.publicKey, self.secretKey = publicKey, None
        
    def encrypt(self, message):
        return encrypt(message, self.generator, 2**self.lgGroupSize, self.publicKey, self.random)
    
    def decrypt(self, ciphertext):
        if self.secretKey is not None:
            return decrypt(ciphertext, self.secretKey, modulus = self.mod)
        raise DecryptionError('No secret key')
    
    def __serialize__(self, buffer):
        # TODO
        pass

    @classmethod
    def __deserialize__(self, buffer):
        # TODO
        pass     
    
if __name__ == '__main__':
    import secrets
    
    from gf2 import GF2
    
    size = 32
    
    groupSize = 2**size
    random = secrets.SystemRandom()
    
    if size == 32:
        mod, generator = (0x199740c05,0xdd9345ba)
        generator = GF2(value=generator, size=size, mod=mod)
    else:
        mod = findRandomIrreduciblePolynomial(size, random)
        generator = findRandomGeneratorPolynomial(size, mod, random)        
    
    publicKey, secretKey = generateKey(generator, groupSize, random)
    
    message = random.randrange(groupSize)
    
    ciphertext = encrypt(message, generator, groupSize, publicKey, random)
    
    messagePrime = decrypt(ciphertext, secretKey, modulus = generator.mod)
    
    print('Using generator %x mod %x' % (generator, mod))
    print('Enc(%x) = (%x,%x)' % (message, ciphertext[0], ciphertext[1]))
    print('Dec((%x,%x)) = %x' % (ciphertext[0], ciphertext[1], messagePrime))
    print()
    
    lgGroupSize=32
    elgamal = ElGamal(lgGroupSize=lgGroupSize)
    generator, mod = elgamal.generator, elgamal.generator.mod
    message = random.randrange(2**lgGroupSize)
    ciphertext = elgamal.encrypt(message)
    
    keys = (int(elgamal.publicKey), int(elgamal.secretKey))
    elgamal2 = ElGamal(lgGroupSize=lgGroupSize, generator=generator, secretKey=keys[1])
    messagePrime = elgamal2.decrypt(ciphertext)    
    
    print('Using generator %x mod %x' % (generator, mod))
    print('Public key: %x, Secret key: %s' % keys)
    print('Enc(%x) = (%x,%x)' % (message, ciphertext[0], ciphertext[1]))
    print('Dec((%x,%x)) = %x' % (ciphertext[0], ciphertext[1], messagePrime))    
    
    '''
    n = 10
    t = n // 2

    
    coefficients = [polynomial.GF2(value=random.randrange(0, 2**size), size=size, mod=mod) for i in range(t+1)]
    gfpoly = polynomial.Polynomial(coefficients = coefficients, mod=mod)
    
    print('%r' % gfpoly)
    
    deal = [gfpoly(i + t) for i in range(n)]
    
    print(deal)
    
    mods = [polynomial.findRandomIrreduciblePolynomial(size, random) for i in range(n)]
    generators = [polynomial.findRandomGeneratorPolynomial(size, mods[i], random) for i in range(n)]
    keys = [generateKey(generators[i], groupSize, random) for i in range(n)]    
    
    print('%s\n%s\n%s' % (mods, generators, keys))
    
    encDeal = [encrypt(deal[i], generators[i], groupSize, keys[i][0], random) for i in range(n)]
    
    print('%r' % encDeal)
    '''