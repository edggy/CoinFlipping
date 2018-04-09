import time
import math
from functools import reduce
from itertools import combinations

from gf2 import GF2
from polynomial import Polynomial

from gf2 import findRandomIrreduciblePolynomial
from gf2 import findRandomGeneratorPolynomial 
    
import ElGamal

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
    8: [(0x15f, 0x1e),
        (0x1a3, 0x98),
        (0x165, 0x36),
        (0x14d, 0xd2),
        (0x1d7, 0x48),
        (0x1f3, 0xad),
        (0x177, 0x6c),
        (0x1d7, 0x31),
        (0x1dd, 0x70),
        (0x1cf, 0x73),
        (0x14d, 0x4b),
        (0x13f, 0xa9),
        (0x12d, 0xd1),
        (0x171, 0x29),
        (0x11b, 0x2c),
        (0x1dd, 0x39),
        (0x169, 0x91),
        (0x139, 0x60),
        (0x177, 0xa5),
        (0x1b1, 0xcf),
        (0x139, 0x5f),
        (0x1c3, 0x54),
        (0x12b, 0x6),
        (0x1cf, 0x54),
        (0x1f5, 0xeb),
        (0x1d7, 0x38),
        (0x1bd, 0xd4),
        (0x165, 0xc2),
        (0x13f, 0x8d),
        (0x177, 0x75),
        (0x165, 0x8a),
        (0x18b, 0x39),
        (0x171, 0x24),
        (0x18d, 0x80),
        (0x177, 0x5b),
        (0x165, 0x87),
        (0x1d7, 0x38),
        (0x1a3, 0x47),
        (0x12b, 0x10),
        (0x169, 0xdb),
        (0x18d, 0x6e),
        (0x18b, 0x75),
        (0x139, 0x3f),
        (0x171, 0x38),
        (0x1f9, 0x82),
        (0x1bd, 0x56),
        (0x187, 0x99),
        (0x1a3, 0xa4),
        (0x12d, 0x43),
        (0x13f, 0x7),
        (0x1e7, 0xd7),
        (0x1f3, 0xfc),
        (0x19f, 0xb5),
        (0x1bd, 0xb9),
        (0x14d, 0x4b),
        (0x19f, 0xa5),
        (0x139, 0x74),
        (0x169, 0x21),
        (0x12b, 0x86),
        (0x171, 0x68),
        (0x14d, 0xaf),
        (0x1c3, 0xf1),
        (0x15f, 0x87)],
    16: [(0x1a2fd, 0x9ae5),
         (0x1ed5f, 0xb46),
         (0x151d3, 0xce65),
         (0x1640d, 0xfe3c),
         (0x1d441, 0xec58),
         (0x17447, 0xca95),
         (0x17515, 0xf3db),
         (0x15851, 0x3955),
         (0x1315d, 0xc8bd),
         (0x183d5, 0x669a),
         (0x1683b, 0xa4d2),
         (0x1bd65, 0x4a5c),
         (0x1ecad, 0x379a),
         (0x1e233, 0x3223),
         (0x1cf09, 0x8c3c),
         (0x12295, 0xa3c6),
         (0x1ece9, 0x2811),
         (0x15289, 0xd23c),
         (0x167a1, 0x625d),
         (0x148f5, 0x6deb),
         (0x1ecd9, 0x2b4f),
         (0x14a6d, 0xa014),
         (0x1d72d, 0x5715),
         (0x19b5d, 0xfff),
         (0x179ab, 0xd40a),
         (0x11127, 0xc5d0),
         (0x19517, 0xb040),
         (0x19055, 0x159c),
         (0x1dd99, 0x59c9),
         (0x1af93, 0x6dc9),
         (0x17a61, 0x5d9d),
         (0x141e1, 0xdbcf),
         (0x1728d, 0x65fa),
         (0x1ef85, 0xfba4),
         (0x158d9, 0xdaa4),
         (0x1df29, 0xa3),
         (0x1e87d, 0xa18),
         (0x18315, 0xbb10),
         (0x1ddf3, 0x30a2),
         (0x10d43, 0xee1b),
         (0x1554b, 0xcc91),
         (0x1e7a5, 0xd89d),
         (0x1db47, 0xb916),
         (0x1cb23, 0xf5f7),
         (0x189ad, 0x4132),
         (0x17e41, 0x49c3),
         (0x1aadd, 0xabc4),
         (0x14c67, 0x3c66),
         (0x17e65, 0x6044),
         (0x1ec2f, 0x49e8),
         (0x19335, 0x781c),
         (0x1e6ad, 0xd69b),
         (0x19d49, 0x1f52),
         (0x1c1df, 0xc189),
         (0x14a91, 0x174),
         (0x16749, 0x52dd),
         (0x1ab3d, 0xcb11),
         (0x14dd1, 0x331e),
         (0x1844d, 0x165a),
         (0x17711, 0x35d7),
         (0x11a6b, 0x1911),
         (0x1b05b, 0x4f79),
         (0x1e439, 0xc646),
         (0x1ba85, 0x1771),
         (0x1706b, 0xf1ba),
         (0x12ec9, 0x7d26),
         (0x161ef, 0x9751),
         (0x182bb, 0xc451),
         (0x1562d, 0x72a6),
         (0x159ff, 0x501e),
         (0x1d8c3, 0x5fd3),
         (0x1ff05, 0xa6c7),
         (0x1e233, 0xdf8e),
         (0x1c527, 0x2b57),
         (0x1c8ef, 0xdc7e),
         (0x183b3, 0x24de),
         (0x1f60f, 0x6150),
         (0x16539, 0xada9),
         (0x1e0d5, 0xfcc0),
         (0x1708f, 0xfd37),
         (0x1a4cd, 0xd030),
         (0x175a7, 0xc1ac),
         (0x1cf33, 0x5dca),
         (0x11bbd, 0x89f9),
         (0x1bc45, 0x5b7e),
         (0x14d39, 0x719b),
         (0x18a0d, 0x51ba),
         (0x1abd9, 0x539b),
         (0x12cad, 0x586a),
         (0x1aa21, 0x6506),
         (0x1733d, 0x3747),
         (0x12051, 0xfd31),
         (0x1c1cd, 0xdb),
         (0x1ea51, 0x49cc),
         (0x1a6b7, 0x41e5),
         (0x1916b, 0xbe1e),
         (0x10175, 0xb0bd),
         (0x15ad7, 0xef65),
         (0x1915d, 0xad36),
         (0x1f93f, 0xe2e8),
         (0x164b9, 0x1047),
         (0x11ae5, 0x4b93),
         (0x15b45, 0x8603),
         (0x18985, 0x5eb2),
         (0x1dd39, 0x9b0a),
         (0x19ba1, 0xeb7),
         (0x115df, 0xe7ae),
         (0x19681, 0xe9fd),
         (0x11bbd, 0x25cf),
         (0x191a1, 0xc093),
         (0x16a3f, 0x1790),
         (0x18637, 0xbc40),
         (0x11c57, 0xb065)],
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
    '''
    Generate keys for ElGamal
    
    @param size - The number of bits of each key
    @param random - The randomness to use to choose keys (requires randrange method)
    @param hardcode - Should we use the hardcoded moduli or generate new moduli
    
    @return - a tuple (m, g, k) where (m, g) is a modulus, generator pair and k is an ElGamal key
    '''
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

class CoinFlipping:
    '''
    An algorithm to generate randomness as long as more than half of the parties are honest
    '''
    def __init__(self, n, lgSize, random):
        
        # Number of parties
        self.n = n
        
        # Max number of corruptions
        self.t = self.n // 2
        
        # The number of bits for each party to generate (for a total of self.size * self.t)
        self.size = lgSize
        
        # The randomness to use
        self.random = random
        
        # The local keys generated
        self.keys = None
        self.publicKeys = None
        self.privateKeys = None
        
        # The randomly generated polynomial in GF2
        self.gfpoly = None
        
        # Points on self.gfpoly
        self.deal = None
        
        # Encrypted points on self.gfpoly
        self.encDeal = None
        
        
        self.summedPoly = None
        self.userWarnings = [None] * self.n
        
    def __repr__(self):
        return '%s' % ((self.n, self.size, self.publicKeys, self.privateKeys, self.gfpoly, self.deal, self.encDeal, self.summedPoly),) 
        
    def generateKeys(self, *, hardcode = False):
        '''
        Generate keys to share with other parties
        '''
        if self.keys is None:
            self.keys = [i for i in map(lambda x: genKey(self.size, self.random, hardcode = hardcode), range(self.n))]
            self.publicKeys = [(int(mod), int(gen), int(key.publicKey)) for mod, gen, key in self.keys]
            self.privateKeys = [int(key.secretKey) for mod, gen, key in self.keys]
        
        return self.publicKeys
        
    def share(self, sharedPublicKeys, polyMod = None, *, _testing = None):
        '''
        Given public keys from other parties generate and share a random polynomial
        '''
        
        if polyMod is None:
            self.polyMod = findRandomIrreduciblePolynomial(self.size, self.random)
        else:
            self.polyMod = polyMod
        
        # Generate t+1 random coefficients
        coefficients = [GF2(value=self.random.randrange(0, 2**self.size), size=self.size, mod=self.polyMod) for i in range(self.t+1)]
        
        # Used for testing to allow a party to change the degree of the polynomial
        if _testing is not None:
            if 'degree' in _testing:
                coefficients = [GF2(value=self.random.randrange(0, 2**self.size), size=self.size, mod=self.polyMod) for i in range(_testing['degree']+1)]
            elif 'coefficients' in _testing:
                coefficients = [GF2(value=i, size=self.size, mod=self.polyMod) for i in _testing['coefficients']]
            elif 'polynomial' in _testing:
                self.gfpoly = _testing['polynomial']
            
        if self.gfpoly is None:
            # Create the polynomial
            self.gfpoly = Polynomial(coefficients = coefficients)        
        
        # Deal out the polynomial
        self.deal = [self.gfpoly(i + self.t + 1) for i in range(self.n)]        
        
        # Determine which keys to use
        sharedPublicKeys = [ElGamal.ElGamal(generator=GF2(value=generator, size=self.size, mod=mod), lgGroupSize=self.size, random=self.random, publicKey=GF2(value=publicKey, size=self.size, mod=mod)) for mod, generator, publicKey in sharedPublicKeys]
        
        # Encrypt each share with the apropriate public key
        self.encDeal = [tuple(map(int, sharedPublicKeys[i].encrypt(self.deal[i]))) for i in range(self.n)]
        
        return self.encDeal
        
    def reconstruct(self, encShares, sharedPublicKeys, sharedSecretKeys, polyMod):
        '''
        @param list<list<int>> encShares - An array of encrypted shares to be reconstructed
        '''
        # Transpose the encrypted shares array so that each row (instead of each column) can be decrypted by a single user
        encShares = list(zip(*encShares))
        
        # Decrypt all of the shares
        shares = []
        shareIndex = 0
        GF2GenPoly = lambda x: GF2(value=x, size=self.size, mod=polyMod)
        for publicKeyRow, secretKeyRow, encSharesRow in zip(sharedPublicKeys, sharedSecretKeys, encShares):
            
            # Check that all data is available
            if publicKeyRow is None or secretKeyRow is None or encSharesRow is None:
                shares.append([None] * self.n)
                self.userWarnings[shareIndex] = 'Aborted'
                shareIndex += 1
                continue
            
            # A row of decrypted shares
            sharesRow = []
            
            for publicKey, secretKey, encShare in zip(publicKeyRow, secretKeyRow, encSharesRow):
                # Check that all data is available
                if publicKey is None or secretKey is None or encShare is None:
                    sharesRow.append(None)
                    self.userWarnings[shareIndex] = 'Aborted'
                    continue
                
                # Seperate the public key
                mod, generator, publicKey = publicKey
                
                GF2Gen = lambda x: GF2(value=x, size=self.size, mod=mod)
                
                # Cast generator, secretKey and encShare into GF2 elements
                generator = GF2Gen(generator)
                secretKey = GF2Gen(secretKey)
                encShare = tuple(map(GF2Gen, encShare))
                
                key = ElGamal.ElGamal(generator=generator, lgGroupSize=self.size, random=self.random, secretKey=secretKey)
                
                # Unique witness detection (check that the public key generated from the secret key is the same as the original public key)
                if int(key.publicKey) != publicKey:
                    sharesRow.append(None)
                    self.userWarnings[shareIndex] = 'Malicious'
                    continue
                
                # Decrypt the share
                share = int(key.decrypt(encShare))
                
                # Cast share into a GF2 element
                share = GF2GenPoly(share)
                
                # Add share to row
                sharesRow.append((GF2GenPoly(shareIndex + self.t + 1), share))
                
            shareIndex += 1
            shares.append(sharesRow)
            
        # Transpose the shares so each row corrosponds to a polynomial
        shares = list(zip(*shares))
        
        # Remove null values
        pointList = [list(filter(lambda x: x is not None, points)) for points in shares]
        
        # Use the first t+1 points to interpolate a unique polynomial
        polynomials = [interpolatePolynomial(points, polyMod, self.size) for points in pointList]
        
        # Sum all of the valid polynomials together
        self.summedPoly = Polynomial(coefficients = [GF2GenPoly(0)]) 
        
        # Check that polynomials are of the correct degree
        for i, poly in enumerate(polynomials):
            if poly.degree() > self.t:
                polynomials[i] = None
                self.userWarnings[i] = 'Malicious'
                combos = combinations(reversed(pointList[i]), self.t+2)
                for combo in combos:
                    newPoly = interpolatePolynomial(list(combo), polyMod, self.size)
                    if newPoly.degree() <= self.t:
                        print(combo)
                        print(newPoly)
            else:
                self.summedPoly += poly
        
        # Evaluate and concatinate the sum of all of the valid polynomials
        return reduce(lambda x, y: x + int(y).to_bytes(math.ceil(self.size/8), 'big'), (self.summedPoly(i) for i in range(self.t)), b'')
        
                
if __name__ == '__main__':       
    def keygen(partyData, *, hardcode = False, verbose = False):
        def keygenUser(name, data):
            if not isinstance(name, int):
                return  
            if verbose: 
                print('Generating keys for %d' % name)
            data['keys'] = {(i, name):key for i, key in enumerate(data['ss'].generateKeys(hardcode=hardcode))}
            for oName, oData in partyData.items():
                oData['sharedKeys'][name] = data['keys'][(oName, name)]        
        
        [i for i in map(lambda x: keygenUser(x[0], x[1]), partyData.items())]
        return partyData
            
    def shareData(partyData, polyMod, *, _testing = None, verbose = False ):
        for name, data in partyData.items():
            if verbose: 
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
        
        publicSS = CoinFlipping(n, lgSize, None)
        publicRandomness = publicSS.reconstruct(encShares, sharedPublicKeys, sharedSecretKeys, polyMod)    
    
        return publicSS, publicRandomness
    
    def testing(n = 20, lgSize = 8):
        partyData = {i:{'ss':CoinFlipping(n, lgSize, random), 'keys':{}, 'sharedKeys':[None]*n} for i in range(n)}
        
        partyData = keygen(partyData)
        polyMod = findRandomIrreduciblePolynomial(lgSize, random)
        
        partyData = shareData(partyData)
        
        publicSS, publicRandomness = recon(partyData, [])
        
        print(publicSS)
        print(publicRandomness)
        print(publicSS.userWarnings)
        

def statistics(n, lgSize = 8, file = None, sortby = 'cumulative', 
               statsRegex = '(\\(generateKeys\\)|\\(share\\)|\\(reconstruct\\))', 
               hardcode = True, seed = None, verbose = False):
    import random
    import cProfile, pstats
    from io import StringIO
    #from collections import defaultdict
    
    if seed is None:
        seed = random.randrange(0, 2**32)
        random.seed(seed)
    elif seed == 'SystemRandom':
        seed = -1
        random = random.SystemRandom()
    else:
        random.seed(seed)
        
    buffer = StringIO()
    
    pr = cProfile.Profile()
    
    #data = defaultdict(lambda: defaultdict(int))
    data = {}
    data['n'] = n
    data['lgSize'] = lgSize
    data['seed'] = seed
    data['hardcode'] = hardcode
    
    buffer.write('n = %s, lgSize = %s, seed = 0x%x, hardcode = %s\n' % (n, lgSize, seed, hardcode))
    
    partyData = {i:{'ss':CoinFlipping(n, lgSize, random), 'keys':{}, 'sharedKeys':[None]*n} for i in range(n)}
    data['partyData'] = partyData
    
    pr.enable()
    partyData = keygen(partyData, hardcode=hardcode, verbose=verbose) 
    pr.disable()
    
    # TODO: move hardcodedKeys to seperate file
    if hardcode and lgSize in hardcodedKeys:
        polyMod, g = hardcodedKeys[lgSize][random.randrange(len(hardcodedKeys[lgSize]))]
    else:
        polyMod = findRandomIrreduciblePolynomial(lgSize, random)
    
    pr.enable()
    partyData = shareData(partyData, polyMod, verbose=verbose)
    pr.disable()
    
    pr.enable()
    publicSS, publicRandomness = recon(partyData, n = n, lgSize = lgSize, polyMod = polyMod, badParties = [])
    pr.disable()
    
    data['publicSS'] = publicSS
    data['publicRandomness'] = publicRandomness
    data['publicPoly'] = publicSS.summedPoly
    data['userWarnings'] = publicSS.userWarnings
    
    ps = pstats.Stats(pr, stream=buffer).sort_stats(sortby)
    ps.print_stats(statsRegex)    
    
    for line in buffer.getvalue().split('\n'):
        functs = ('generateKeys', 'share', 'reconstruct')
        toks = line.split()
        if 'function calls' in line:
            data['functionCalls'] = int(toks[0])
        
        for funName in functs:
            if ('(%s)' % funName) in line:
                # ncalls  tottime  percall  cumtime  percall
                data['%s_ncalls' % funName] = int(toks[0])
                data['%s_tottime'  % funName] = float(toks[1])
                data['%s_tottime_percall'  % funName] = float(toks[2])
                data['%s_cumtime'  % funName] = float(toks[3])
                data['%s_cumtime_percall'  % funName] = float(toks[4])
    
    buffer.write('%r\n' % publicSS)
    buffer.write('%r\n' % publicRandomness)
    buffer.write('%r\n' % publicSS.userWarnings)
    
    totalBits = 0
    for i in range(n):
        encShares = partyData[i]['shares']
        for encShare in encShares:
            totalBits += reduce(lambda x, y: x + y.bit_length(), encShare, 0)
        sharedPublicKeys = partyData[i]['ss'].publicKeys
        for publicKey in sharedPublicKeys:
            totalBits += reduce(lambda x, y: x + y.bit_length(), publicKey, 0)
        sharedSecretKeys = partyData[i]['ss'].privateKeys
        for secretKey in sharedSecretKeys:
            totalBits += secretKey.bit_length()
    
    buffer.write('Total storage needed = %s bits\n' % totalBits)
    buffer.write('Bits generated = %s bits\n' % (len(publicRandomness)*8,))
    buffer.write('Storage per bit = %s\n' % (totalBits/(len(publicRandomness)*8),))    
    buffer.write('\n')
    file.write(buffer.getvalue())
    file.flush()
    
    data['storage'] = totalBits
    data['totalBits'] = len(publicRandomness)*8
    data['storagePerBit'] = totalBits/(len(publicRandomness)*8)
                             
    return file, data
    
def doExperiment(ns = (4,8,16,32,64,128), lgSizes = (8, 16, 32), file = None, csvFile = None, seeds = (0,), repeats = 10, verbose = False):
    import csv
    import os
    header = next(csv.reader(csvFile))
    csvFile.seek(0, os.SEEK_END)
    for n in ns:
        for lgSize in lgSizes:
            for seed in seeds:
                for rep in range(repeats):
                    print('n = %s, lgSize = %s, seed = %s, rep = %s' % (n, lgSize, seed, rep))
                    output, data = statistics(n, lgSize=lgSize, file=file, seed=seed, verbose=verbose)
                    if file is None:
                        print(output.getvalue())
                    if csvFile is not None:
                        
                        w = csv.DictWriter(csvFile, header, extrasaction='ignore')
                        #w.writeheader()
                        w.writerow(data)                        
    

def makeKeys(size=32, amount = None):
    print('Generating Polynomials')
    count = 0
    with open('polys-%s' % size, 'a') as f:
        while amount is None or count < amount:
            m = getMod(size, random)
            print(hex(m), end=' ', flush=True)
            g = getGen(m, size, random)
            print(hex(g), flush=True)
            f.write('%x\t%x\n' % (m, g))
            f.flush()   
            count += 1

if __name__ == '__main__':
    #import random
    #random.seed(0)
    #random = random.SystemRandom()
    #from collections import defaultdict    
    #import cProfile
    
    #makeKeys(size=32, amount=None)
    with open('tmpStats.txt' , 'a') as f, open('tmpStats.csv' , 'r+', newline='') as g:
        #statistics(32, lgSize = 8, file = f, sortby = 'cumulative', hardcode = True, seed = 0)
        doExperiment(file=f, csvFile=g, ns=(4,8,16,32,64,128,256,512), lgSizes=(8,16,32), seeds=(0,1,2,None,'SystemRandom'), repeats=10, verbose=True)
        
    '''n = 32
    lgSize = 16
    hardcode = True
    
    print('n = %s, lgSize = %s, hardcode = %s' % (n, lgSize, hardcode))
    
    partyData = {i:{'ss':CoinFlipping(n, lgSize, random), 'keys':{}, 'sharedKeys':[None]*n} for i in range(n)}
    
    cProfile.run('partyData = keygen(partyData, hardcode=hardcode)', sort='cumulative')
    
    # TODO: move hardcodedKeys to seperate file
    if hardcode and lgSize in hardcodedKeys:
        polyMod, g = hardcodedKeys[lgSize][random.randrange(len(hardcodedKeys[lgSize]))]
    else:
        polyMod = findRandomIrreduciblePolynomial(lgSize, random)
    
    GF2GenPoly = lambda x: GF2(value=x, size=lgSize, mod=polyMod)
    makePoly = lambda x, start = 1: interpolatePolynomial([(GF2GenPoly(i + n//2 + start), GF2GenPoly(v)) for i,v in enumerate(x)] , polyMod, lgSize)
    badPoints0 = [n//2 + 1] + [i + n//2 + 1 for i in range(1,n//2)]
    badPoly0 = makePoly(badPoints0)
    """badPoints0A = [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    badPoly0A = makePoly(badPoints0A)
    badPoints0B = [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    badPoly0B = makePoly(badPoints0B, 2)    
    badPoly0 = badPoly0A + badPoly0B"""
    
    badPoints1 = [2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    badPoly1 = makePoly(badPoints1)
    
    shareStats = cProfile.run('partyData = shareData(partyData, polyMod, _testing = {0:{"polynomial":badPoly0},1:{"polynomial":badPoly1}})', sort='cumulative')
    
    print(badPoly0)
    print(badPoly1)    
    cProfile.run('publicSS, publicRandomness = recon(partyData, n = n, lgSize = lgSize, polyMod = polyMod, badParties = [2,3,4])', sort='cumulative')

    print(publicSS)
    print(publicRandomness)
    print(publicSS.userWarnings)
    
    totalBits = 0
    for i in range(n):
        encShares = partyData[i]['shares']
        for encShare in encShares:
            totalBits += reduce(lambda x, y: x + y.bit_length(), encShare, 0)
        sharedPublicKeys = partyData[i]['ss'].publicKeys
        for publicKey in sharedPublicKeys:
            totalBits += reduce(lambda x, y: x + y.bit_length(), publicKey, 0)
        sharedSecretKeys = partyData[i]['ss'].privateKeys
        for secretKey in sharedSecretKeys:
            totalBits += secretKey.bit_length()
    
    print('Total storage needed = %s bits' % totalBits)
    print('Bits generated = %s bits' % (len(publicRandomness)*8,))
    print('Storage per bit = %s' % (totalBits/len(publicRandomness)*8,))
    '''