#!/usr/bin/env python3
### Curve25519 stuff from RFC_7748
import binascii

def FP(p):
    # caveat: caller must ensure that p is prime
    class F:
        def __init__(self,x):
            self.int = x % p
        def __str__(self):
            return str(self.int)
        __repr__ = __str__
        def __eq__(a,b):
            return a.int == b.int
        def __ne__(a,b):
            return a.int != b.int
        def __add__(a,b):
            return F(a.int + b.int)
        def __sub__(a,b):
            return F(a.int - b.int)
        def __mul__(a,b):
            return F(a.int * b.int)
        def __truediv__(a,b):
            # caveat: caller must ensure that b is nonzero
            # Use Fermat little theorem
            return a*F(pow(b.int,p-2,p))
    return F

def fieldToInt(s, bits):
    bytes = s.int.to_bytes((bits+7)//8,byteorder='little')
    return int.from_bytes(bytes,byteorder='big')

def strToInt(s):
    bytes = binascii.unhexlify(s)
    return int.from_bytes(bytes,byteorder='little')

def prepScalar25519(k):
    maskAnd = ((1<<256)-1) ^ ((127 ^ 0xFF)<<(8*31)) ^ 248 ^ 0xFF
    maskOr = (64<<(8*31))
    k &= maskAnd
    k |= maskOr
    return k

def prepScalar448(k):
    maskAnd = ((1<<448)-1) ^ 252 ^ 0xFF
    maskOr = (128<<(8*55))
    k &= maskAnd
    k |= maskOr
    return k

def cswap(swap, x_2, x_3):
    """shall be implemented in constant time (so NOT like it is now)"""
    if swap:
        tmp = x_2
        x_2 = x_3
        x_3 = tmp
    return (x_2, x_3)

def MontgommeryScalarMul(field,k,u,bits,a24):
    x_1 = u
    x_2 = field(1)
    z_2 = field(0)
    x_3 = u
    z_3 = field(1)
    swap = 0

    for t in range(bits-1, -1,-1):
        k_t = (k.int >> t) & 1
        swap ^= k_t
        (x_2, x_3) = cswap(swap, x_2, x_3)
        (z_2, z_3) = cswap(swap, z_2, z_3)
        swap = k_t

        A = x_2 + z_2
        AA = A * A
        B = x_2 - z_2
        BB = B * B
        E = AA - BB
        C = x_3 + z_3
        D = x_3 - z_3
        DA = D * A
        CB = C * B
        x_3 = (DA + CB)
        x_3 = x_3 * x_3
        z_3 = (DA - CB)
        z_3 = z_3 * z_3
        z_3 = x_1 * z_3
        x_2 = AA * BB
        z_2 = E * (AA + a24 * E)

    (x_2, x_3) = cswap(swap, x_2, x_3)
    (z_2, z_3) = cswap(swap, z_2, z_3)
    return x_2 / z_2

def X25519(k,u):
    bits = 255
    p = (1<<255) - 19
    field = FP(p)
    a24 = field(121665)
    k=field(prepScalar25519(k))
    u=field(u)
    return MontgommeryScalarMul(field,k,u,bits,a24).int

def DH25519(private_key,other_public_key=None):
    if other_public_key is None:
        return X25519(private_key,9)
    K = X25519(private_key,other_public_key)
    assert(K != 0)
    return K

def X448(k,u):
    bits = 448
    p = (1<<448) - (1<<224) - 1
    field = FP(p)
    a24 = field(39081)
    k=field(prepScalar448(k))
    u=field(u)
    return MontgommeryScalarMul(field,k,u,bits,a24).int

def DH448(private_key,other_public_key=None):
    if other_public_key is None:
        return X448(private_key,5)
    K = X448(private_key,other_public_key)
    assert(K != 0)
    return K

def test25519():
    bits = 255
    p = (1<<255) - 19
    field = FP(p)
    a24 = field(121665)
    expected=field(strToInt('c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552'))
    k=field(prepScalar25519(31029842492115040904895560451863089656472772604678260265531221036453811406496))
    u=field(34426434033919594451155107781188821651316167215306631574996226621102155684838)
    s=MontgommeryScalarMul(field,k,u,bits,a24)
    print("e=%x"%expected.int)
    print("s=%x"%s.int)
    assert(s==expected)
    expected=field(strToInt('95cbde9476e8907d7aade45cb4b873f88b595a68799fa152e6f8f7647aac7957'))
    k=field(prepScalar25519(35156891815674817266734212754503633747128614016119564763269015315466259359304))
    u=field(8883857351183929894090759386610649319417338800022198945255395922347792736741)
    s=MontgommeryScalarMul(field,k,u,bits,a24)
    print("e=%x"%expected.int)
    print("s=%x"%s.int)
    assert(s==expected)
    expected1 = field(strToInt('422c8e7a6227d7bca1350b3e2bb7279f7897b87bb6854b783c60e80311ae3079'))
    expected1000 = field(strToInt('684cf59ba83309552800ef566f2f4d3c1c3887c49360e3875f2eb94d99532c51'))
    k=field(strToInt('0900000000000000000000000000000000000000000000000000000000000000'))
    u = k
    for i in range(0,1000):
        oldk = k
        k=field(prepScalar25519(k.int))
        k = MontgommeryScalarMul(field,k,u,bits,a24)
        u = oldk
    e = expected1000
    print("e=%x"%e.int)
    print("k=%x"%k.int)
    assert(k==e)

def test25519_DH():
    alice_private_key = strToInt('77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a')
    alice_public_key = DH25519(alice_private_key)
    print("alice_public_key=%x"%alice_public_key)
    bob_private_key = strToInt('5dab087e624a8a4b79e17f8b83800ee66f3bb1292618b6fd1c2f8b27ff88e0eb')
    bob_public_key = DH25519(bob_private_key)
    print("bob_public_key  =%x"%bob_public_key)
    alice_K = DH25519(alice_private_key,bob_public_key)
    bob_K = DH25519(bob_private_key,alice_public_key)
    assert(alice_K==bob_K)
    expected = strToInt('4a5d9d5ba4ce2de1728e3bf480350f25e07e21c947d19e3376f09b3c1e161742')
    assert(bob_K==expected)

def test448():
    bits = 448
    p = (1<<448) - (1<<224) - 1
    field = FP(p)
    a24 = field(39081)
    expected=field(strToInt('ce3e4ff95a60dc6697da1db1d85e6afbdf79b50a2412d7546d5f239fe14fbaadeb445fc66a01b0779d98223961111e21766282f73dd96b6f'))
    k=field(prepScalar448(599189175373896402783756016145213256157230856085026129926891459468622403380588640249457727683869421921443004045221642549886377526240828))
    u=field(382239910814107330116229961234899377031416365240571325148346555922438025162094455820962429142971339584360034337310079791515452463053830)
    s=MontgommeryScalarMul(field,k,u,bits,a24)
    print("e=%x"%expected.int)
    print("s=%x"%s.int)
    assert(s==expected)
    expected=field(strToInt('884a02576239ff7a2f2f63b2db6a9ff37047ac13568e1e30fe63c4a7ad1b3ee3a5700df34321d62077e63633c575c1c954514e99da7c179d'))
    k=field(prepScalar448(633254335906970592779259481534862372382525155252028961056404001332122152890562527156973881968934311400345568203929409663925541994577184))
    u=field(622761797758325444462922068431234180649590390024811299761625153767228042600197997696167956134770744996690267634159427999832340166786063)
    s=MontgommeryScalarMul(field,k,u,bits,a24)
    print("e=%x"%expected.int)
    print("s=%x"%s.int)
    assert(s==expected)
    expected1 = field(strToInt('3f482c8a9f19b01e6c46ee9711d9dc14fd4bf67af30765c2ae2b846a4d23a8cd0db897086239492caf350b51f833868b9bc2b3bca9cf4113'))
    expected1000 = field(strToInt('aa3b4749d55b9daf1e5b00288826c467274ce3ebbdd5c17b975e09d4af6c67cf10d087202db88286e2b79fceea3ec353ef54faa26e219f38'))
    k=field(strToInt('0500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000'))
    u = k
    for i in range(0,1000):
        oldk = k
        k=field(prepScalar448(k.int))
        k = MontgommeryScalarMul(field,k,u,bits,a24)
        u = oldk
    e = expected1000
    print("e=%x"%e.int)
    print("k=%x"%k.int)
    assert(k==e)

def test448_DH():
    alice_private_key = strToInt('9a8f4925d1519f5775cf46b04b5800d4ee9ee8bae8bc5565d498c28dd9c9baf574a9419744897391006382a6f127ab1d9ac2d8c0a598726b')
    alice_public_key = DH448(alice_private_key)
    print("alice_public_key=%x"%alice_public_key)
    bob_private_key = strToInt('1c306a7ac2a0e2e0990b294470cba339e6453772b075811d8fad0d1d6927c120bb5ee8972b0d3e21374c9c921b09d1b0366f10b65173992d')
    bob_public_key = DH448(bob_private_key)
    print("bob_public_key  =%x"%bob_public_key)
    alice_K = DH448(alice_private_key,bob_public_key)
    bob_K = DH448(bob_private_key,alice_public_key)
    assert(alice_K==bob_K)
    expected = strToInt('07fff4181ac6cc95ec1c16a94a0f74d12da232ce40a77552281d282bb60c0b56fd2464c335543936521c24403085d59a449a5037514a879d')
    assert(bob_K==expected)

if __name__ == "__main__":
    test448()
    test25519()
    test448_DH()
    test25519_DH()
