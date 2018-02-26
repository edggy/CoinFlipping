import time
import math
from functools import reduce

from gf2.gf2 import GF2
from polynomial import Polynomial

try:
    from gf2c import findRandomIrreduciblePolynomial
    from gf2c import findRandomGeneratorPolynomial
except ImportError:
    from gf2.gf2 import findRandomIrreduciblePolynomial
    from gf2.gf2 import findRandomGeneratorPolynomial 
    
import ElGamal

#import concurrent.futures

try:
    from interpolateGF2 import interpolatePolynomial as interpolate
    interpolatePolynomial = lambda points, polyMod, size: Polynomial(coefficients=[GF2(value=i, size=size, mod=polyMod) for i in interpolate(points, polyMod)])
except ImportError:
    from polynomial import interpolatePolynomial as interpolate
    interpolatePolynomial = lambda points, polyMod, size: Polynomial(coefficients=interpolate(points))[0]

getMod = lambda size, random: findRandomIrreduciblePolynomial(size, random)
getGen = lambda mod, size, random: findRandomGeneratorPolynomial(size, mod, random)
getKey = lambda gen, size, random: ElGamal.ElGamal(generator=gen, lgGroupSize=size, random=random)

hardcodedKeys = {
    32: [(0x1020609b3, 0xf2dff8e8),
         (0x10ae3a5b5,0xb3a61eea),
         (0x10c44a745,0xb46b3c1a),
         (0x112db649d,0x9e10e810),
         (0x112e52541,0xc7742a13),
         (0x118651375,0x510b540a),
         (0x11950d75b,0x9341edbe),
         (0x122d26c33,0x6c3e5484),
         (0x125f3797b,0x8a9b0b93),
         (0x12828e52f,0xf4fd377e),
         (0x12d232ce1,0xf16c0589),
         (0x12d36bf61,0xe1f65ad3),
         (0x1306a6fa5,0x41e23c41),
         (0x131ca7e5b,0xab469edb),
         (0x1322f2dd7,0x57e61cc6),
         (0x132cc9d0f,0xb3ecc0b2),
         (0x136941fdd,0xe0b7a98f),
         (0x138dfd1bb,0x2632ae57),
         (0x13e0d9af5,0x785fdb2f),
         (0x1494df651,0xf87b1bf0),
         (0x14b1aa05d,0x2b2278dc),
         (0x14e077749,0x3b714eac),
         (0x15c31b23b,0xbb248f8d),
         (0x16a66447d,0xee283233),
         (0x16df76d31,0x738763c9),
         (0x16e16d1ff,0xe3cc25f5),
         (0x170e2a4f9,0x7ad302d1),
         (0x17113ff8f,0x501d7d2),
         (0x1749d92ff,0xe9a7cc7d),
         (0x1785697b3,0x2b59f8fd),
         (0x188e3491b,0x8212121d),
         (0x18cb94dbd,0xc467da7d),
         (0x18f165e83,0x67feb58b),
         (0x1957dcd91,0x48199453),
         (0x195e486cd,0x31479fe5),
         (0x199740c05,0xdd9345ba),
         (0x19b98d4e1,0x947b658c),
         (0x1a688aecd,0xc776e3a0),
         (0x1a85c5fd5,0xf13e16f2),
         (0x1aeb6425f,0xdcf8fe94),
         (0x1c161ab4f,0x2a9b8784),
         (0x1c61943eb,0x78501824),
         (0x1c811fbed,0x893dd8d4),
         (0x1c93bbbf1,0xaefa4d2c),
         (0x1cfde74f5,0x856eec7a),
         (0x1d31abc89,0xe3c0f1bb),
         (0x1e2dbe967,0x844071a2),
         (0x1e1ee55cf,0xbf47c9d0),
         (0x1e41f4bbf,0xcc98d5eb),
         (0x1e684059b,0xa6d827c3),
         (0x1ea9d620d,0xba5cba95),
         (0x1ee2d6291,0xb9740837),
         (0x1f039fedf,0xfa3a5a3d),
         (0x1ffd5e933,0x320c9e72),         
         (0x1ffe77383,0x8b4caadb)]
}

def genKey(size, random, *, hardcode = False):
    m, g = None, None
    if hardcode and size in hardcodedKeys:
        m, g = hardcodedKeys[size][random.randrange(len(hardcodedKeys[size]))]
        g = GF2(value=g, size=size, mod=m)
        
        from gf2.gf2_math import _exteuc
        
        # Find s such that s and 2**size-1 are coprime
        s = random.randrange(1, 2**size)
        while _exteuc(s, 2**size-1)[0] != 1:
            s = random.randrange(1, 2**size)
        
        g = g**s
    
    if m is None or g is None:
        m = getMod(size, random)
        g = getGen(m, size, random)
        
    k = getKey(g, size, random)
    return (m, g, k)

class SecretSharing:
    def __init__(self, n, lgSize, random):
        self.n = n
        self.t = self.n // 2
        self.size = lgSize
        self.random = random
        self.keys = None
        self.publicKeys = None
        self.privateKeys = None
        self.gfpoly = None
        self.deal = None
        self.encDeal = None
        self.summedPoly = None
        self.userWarnings = [None] * self.n
        
    def __repr__(self):
        return '%s' % ((self.n, self.size, self.publicKeys, self.privateKeys, self.gfpoly, self.deal, self.encDeal, self.summedPoly),) 
        
    def generateKeys(self, *, hardcode = False):
        if self.keys is None:
            self.keys = [i for i in map(lambda x: genKey(self.size, self.random, hardcode = hardcode), range(self.n))]
            self.publicKeys = [(int(mod), int(gen), int(key.publicKey)) for mod, gen, key in self.keys]
            self.privateKeys = [int(key.secretKey) for mod, gen, key in self.keys]
        
        return self.publicKeys
        
    def share(self, sharedPublicKeys, polyMod = None, *, _testing = None):
        
        if polyMod is None:
            self.polyMod = findRandomIrreduciblePolynomial(self.size, self.random)
        else:
            self.polyMod = polyMod
        
        coefficients = [GF2(value=self.random.randrange(0, 2**self.size), size=self.size, mod=self.polyMod) for i in range(self.t+1)]
        
        if _testing is not None and 'degree' in _testing:
            coefficients = [GF2(value=self.random.randrange(0, 2**self.size), size=self.size, mod=self.polyMod) for i in range(_testing['degree']+1)]
            
        self.gfpoly = Polynomial(coefficients = coefficients)        
        
        self.deal = [self.gfpoly(i + self.t + 1) for i in range(self.n)]        
        
        sharedPublicKeys = [ElGamal.ElGamal(generator=GF2(value=generator, size=self.size, mod=mod), lgGroupSize=self.size, random=self.random, publicKey=GF2(value=publicKey, size=self.size, mod=mod)) for mod, generator, publicKey in sharedPublicKeys]

        self.encDeal = [tuple(map(int, sharedPublicKeys[i].encrypt(self.deal[i]))) for i in range(self.n)]
        
        return self.encDeal
        
    def reconstruct(self, encShares, sharedPublicKeys, sharedSecretKeys, polyMod):
        '''
        param list<list<int>> encShares: An array of encrypted shares to be reconstructed
        '''
        # Transpose the encrypted shares array so that each row (instead of each column) can be decrypted by a single user
        encShares = list(zip(*encShares))
        
        # Decrypt all of the shares
        shares = []
        shareIndex = 0
        GF2GenPoly = lambda x: GF2(value=x, size=self.size, mod=polyMod)
        for pkRow, secretKeyRow, encSharesRow in zip(sharedPublicKeys, sharedSecretKeys, encShares):
            
            if pkRow is None or secretKeyRow is None or encSharesRow is None:
                shares.append([None] * self.n)
                self.userWarnings[shareIndex] = 'Aborted'
                shareIndex += 1
                continue
            
            sharesRow = []
            
            for pk, secretKey, encShare in zip(pkRow, secretKeyRow, encSharesRow):
                if pk is None or secretKey is None or encShare is None:
                    sharesRow.append(None)
                    self.userWarnings[shareIndex] = 'Aborted'
                    continue
                
                mod, generator, publicKey = pk
                
                GF2Gen = lambda x: GF2(value=x, size=self.size, mod=mod)
                
                generator = GF2Gen(generator)
                secretKey = GF2Gen(secretKey)
                encShare = tuple(map(GF2Gen, encShare))
                
                key = ElGamal.ElGamal(generator=generator, lgGroupSize=self.size, random=self.random, secretKey=secretKey)
                if int(key.publicKey) != publicKey:
                    sharesRow.append(None)
                    self.userWarnings[shareIndex] = 'Malicious'
                    continue
                
                share = int(key.decrypt(encShare))
                share = GF2GenPoly(share)
                sharesRow.append((GF2GenPoly(shareIndex + self.t + 1), share))
                
            shareIndex += 1
            shares.append(sharesRow)
            
        # Transpose the shares so each row corrosponds to a polynomial
        shares = list(zip(*shares))
        
        pointList = [list(filter(lambda x: x is not None, points)) for points in shares]

        polynomials = [interpolatePolynomial(points[:self.t+2], polyMod, self.size) for points in pointList]
        
        self.summedPoly = Polynomial(coefficients = [GF2GenPoly(0)]) 
        for i, poly in enumerate(polynomials):
            if poly.degree() > self.t:
                polynomials[i] = None
                self.userWarnings[i] = 'Malicious'
            else:
                self.summedPoly += poly
        
        return reduce(lambda x, y: x + int(y).to_bytes(math.ceil(self.size/8), 'big'), (self.summedPoly(i) for i in range(self.t)), b'')
        
                
if __name__ == '__main__':       
    def keygen(partyData, *, hardcode = False):
        def keygenUser(name, data):
            if not isinstance(name, int):
                return    
            print('Generating keys for %d' % name)
            data['keys'] = {(i, name):key for i, key in enumerate(data['ss'].generateKeys(hardcode=hardcode))}
            for oName, oData in partyData.items():
                oData['sharedKeys'][name] = data['keys'][(oName, name)]        
        
        [i for i in map(lambda x: keygenUser(x[0], x[1]), partyData.items())]
        return partyData
            
    def shareData(partyData, polyMod, *, _testing = None):
        for name, data in partyData.items():
            print('Sharing data for %d' % name)
            # Simulate users generating a polynomial that's too big
            if _testing is not None and name in _testing:
                data['shares'] = data['ss'].share(data['sharedKeys'], polyMod=polyMod, _testing = _testing[name])
            else:
                data['shares'] = data['ss'].share(data['sharedKeys'], polyMod=polyMod)
        return partyData
    
    def recon(partyData, n, lgSize, polyMod, *, badParties = []):
        encShares = [partyData[i]['shares'] for i in range(n)]
        sharedPublicKeys = [partyData[i]['ss'].publicKeys for i in range(n)]
        sharedSecretKeys = [partyData[i]['ss'].privateKeys for i in range(n)]
        
        # Simulate parties aborting (or having invalid signatures)
        for badParty in badParties:
            sharedSecretKeys[badParty] = None
        
        publicSS = SecretSharing(n, lgSize, random)
        publicRandomness = publicSS.reconstruct(encShares, sharedPublicKeys, sharedSecretKeys, polyMod)    
    
        return publicSS, publicRandomness
    
    def testing(n = 20, lgSize = 8):
        partyData = {i:{'ss':SecretSharing(n, lgSize, random), 'keys':{}, 'sharedKeys':[None]*n} for i in range(n)}
        
        partyData = keygen(partyData)
        polyMod = findRandomIrreduciblePolynomial(lgSize, random)
        
        partyData = shareData(partyData)
        
        publicSS, publicRandomness = recon(partyData, [])
        
        print(publicSS)
        print(publicRandomness)
        print(publicSS.userWarnings)
        
        #print('Reconstruct time: %.3f, per user: %.3f' % (reconsrtuctTime, reconsrtuctTime/n))
    
if __name__ == '__main__':
    import random
    #random.seed(0)
    random = random.SystemRandom()
    from collections import defaultdict    
    import cProfile
    
    makeKeys = False
    
    if makeKeys:
        print('Generating Polynomials')
        size = 32
        with open('polys-%s' % size, 'a') as f:
            while True:
                m = getMod(size, random)
                print(hex(m), end=' ', flush=True)
                g = getGen(m, size, random)
                print(hex(g), flush=True)
                f.write('%x\t%x\n' % (m, g))
                f.flush()
        
    n = 100
    lgSize = 32
    hardcode = True
    
    partyData = {i:{'ss':SecretSharing(n, lgSize, random), 'keys':{}, 'sharedKeys':[None]*n} for i in range(n)}
    
    cProfile.run('partyData = keygen(partyData, hardcode=hardcode)', sort='cumulative')
    
    # TODO: move hardcodedKeys to seperate file
    if hardcode and lgSize in hardcodedKeys:
        polyMod, g = hardcodedKeys[lgSize][random.randrange(len(hardcodedKeys[lgSize]))]
    else:
        polyMod = findRandomIrreduciblePolynomial(lgSize, random)
    
    cProfile.run('partyData = shareData(partyData, polyMod)', sort='cumulative')
    
    cProfile.run('publicSS, publicRandomness = recon(partyData, n = n, lgSize = lgSize, polyMod = polyMod, badParties = [])', sort='cumulative')

    print(publicSS)
    print(publicRandomness)
    print(publicSS.userWarnings)    