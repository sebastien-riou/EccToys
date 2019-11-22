#!/usr/bin/env python3

from hashlib import sha256
import hmac

from collections import namedtuple
Point = namedtuple("Point", "x y")

class ModArithmetic(object):
    def __init__(self, modulus):
        self.modulus = modulus

    def add(self, a, b):
        return (a+b)%self.modulus

    def sub(self, a, b):
        o = a-b
        if o<0:
            o+=self.modulus
        return o % self.modulus

    def mul(self, a, b):
        return (a*b)%self.modulus

    def div(self, a, b):
        return (a*self.inv(b))%self.modulus

    def inv(self, a):
        return ModArithmetic.modinv(a,self.modulus)

    @staticmethod
    def egcd(a, b):
        if a == 0:
            return (b, 0, 1)
        else:
            g, y, x = ModArithmetic.egcd(b % a, a)
            return (g, x - (b // a) * y, y)

    @staticmethod
    def modinv(a, m):
        g, x, y = ModArithmetic.egcd(a, m)
        if g != 1:
            raise Exception('modular inverse does not exist')
        else:
            return x % m



class GroupOperation(object):
    O = 'Origin'

    def __init__(self, a,b,modulus):
        self.a=a
        self.b=b
        self._ = ModArithmetic(modulus)

    def modulus(self):
        return self._.modulus

    def IsNeutralPoint(self,p):
        return p==self.O

    def __PointDouble(self,a):
        if a == self.O:
            return self.O
        x = a.x
        y = a.y
        x2 = self._.mul(x,x)
        s = self._.mul(3,x2)
        s = self._.add(s,self.a)
        d = self._.mul(2,y)
        try:
            s = self._.div(s,d)
        except:
            return self.O
        ox = self._.mul(s,s)
        ox = self._.sub(ox,x)
        ox = self._.sub(ox,x)
        oy = self._.sub(x,ox)
        oy = self._.mul(s,oy)
        oy = self._.sub(oy,y)
        return Point(ox,oy)

    def PointAdd(self,a,b):
        if a == self.O:
            return b
        if b == self.O:
            return a
        if a == b:
            return self.__PointDouble(a)
        x1 = a.x
        x2 = b.x
        y1 = a.y
        y2 = b.y
        n = self._.sub(y2,y1)
        d = self._.sub(x2,x1)
        try:
            s = self._.div(n,d)
        except:
            return self.O
        ox = self._.mul(s,s)
        ox = self._.sub(ox,x1)
        ox = self._.sub(ox,x2)
        oy = self._.sub(x1,ox)
        oy = self._.mul(s,oy)
        oy = self._.sub(oy,y1)
        return Point(ox,oy)

    def ScalarMul(self,scalar,p):
        #print(scalar,"*",p,": ",end="")
        t = p
        l = scalar.bit_length()
        for i in range(l-2,-1,-1):
            #print("D",end="")
            t = self.__PointDouble(t)
            if 0 != (scalar & (1<<i)):
                t = self.PointAdd(t,p)
                #print("A",end="")
        #print()
        return t

    def IsOnCurve(self,p):
        if p == self.O:
            return True
        x3 = self._.mul(p.x,p.x)
        x3 = self._.mul(p.x,x3)
        tmp = self._.mul(self.a,p.x)
        tmp = self._.add(self.b,tmp)
        tmp = self._.add(tmp,x3)
        y2 = self._.mul(p.y,p.y)
        return y2==tmp

class GroupOperationSec(object):
    O = 'Origin'

    def __init__(self, a,b,modulus):
        self.a=a
        self.b=b
        self._ = ModArithmetic(modulus)

    def ec_inv(self,a):
        """
        Inverse of the point a on the elliptic curve y^2 = x^3 + ax + b.
        """
        if a == self.O:
            return a
        return Point(a.x, (-a.y)%self._.modulus)

    def __PointAdd_BrierJoye(self,a,b):
        #from Brier and Joye 2002
        #λ = (x1^2 + x1 x2 + x2^2 + a2 x1 + a2 x2 + a4 − a1 y1) / (y1 + y2 + a1 x2 + a3)
        #a4 = self.a
        #a6 = self.b
        #others a = 0
        x1 = a.x
        x2 = b.x
        y1 = a.y
        y2 = b.y

        x1x2 = self._.mul(x1,x2)
        #compute (a+b)^2 - ab rather than a^2+ab+b^2 to save a squaring
        n = self._.add(x1,x2)
        n = self._.mul(n,n)
        n = self._.sub(n,x1x2)

        n = self._.add(n,self.a)
        d = self._.add(y1,y2)

        n_fix = self._.sub(y1,y2)
        d_fix = self._.sub(x1,x2)

        #delta and fix formulaes from "2005 - Stebila and Theriault - Unified Point Addition Formulae and Side-Channel Attacks"
        delta = self._.add(y1,y2)
        delta = self._.add(delta,x1)
        delta = self._.sub(delta,x2)
        if delta != 0:
            #add
            n = self._.add(n,n_fix)
            d = self._.add(d,d_fix)
        else:
            #sub
            n = self._.sub(n,n_fix)
            d = self._.sub(d,d_fix)

        l = self._.div(n,d)

        #x3 = λ^2 + a1 λ − a2 − x1 − x2
        x3 = self._.mul(l,l)
        x3 = self._.sub(x3,x1)
        x3 = self._.sub(x3,x2)
        #μ = y1 − λx1
        #compute -u instead, -u = λx1 - y1
        mu = self._.mul(l,x1)
        mu = self._.sub(mu,y1)
        #y3 = −(λ + a1 )x3 − μ − a3
        #y3 = -u - (λ + a1 )x3 - a3
        y3 = self._.mul(l,x3)
        y3 = self._.sub(mu,y3)
        return Point(x3,y3)

    def __PointAdd_BrierJoye2(self,a,b):
        #from "2005 - Stebila and Theriault - Unified Point Addition Formulae and Side-Channel Attacks"
        #choose m=m0=0 as described in the appendix A.3
        #can't get it to work :-(
        raise("does not work, don't use")
        #a4 = self.a
        #a6 = self.b
        #others a = 0
        x1 = a.x
        x2 = b.x
        y1 = a.y
        y2 = b.y

        x1x2 = self._.mul(x1,x2)
        #n = (x1+x2)^2 + x1x2 + a(x1+x2)
        n = self._.add(x1,x2)
        n = self._.mul(n,n)
        n = self._.sub(n,x1x2)

        tmp = self._.add(x1,x2)
        tmp = self._.mul(tmp,self.a)

        n = self._.add(n,tmp)
        d = self._.add(y1,y2)

        delta = self._.add(y1,y2)
        delta = self._.add(delta,x1)
        if delta != 0:
            #add
            n = self._.add(n,y2)
            d = self._.add(d,x1)
        else:
            delta = self._.add(y1,y2)
            delta = self._.add(delta,x2)
            assert delta != 0
            #sub
            n = self._.sub(n,y1)
            d = self._.sub(d,x2)

        l = self._.div(n,d)

        #x3 = λ^2 + a1 λ − a2 − x1 − x2
        x3 = self._.mul(l,l)
        x3 = self._.sub(x3,x1)
        x3 = self._.sub(x3,x2)
        #μ = y1 − λx1
        #compute -u instead, -u = λx1 - y1
        mu = self._.mul(l,x1)
        mu = self._.sub(mu,y1)
        #y3 = −(λ + a1 )x3 − μ − a3
        #y3 = -u - (λ + a1 )x3 - a3
        y3 = self._.mul(l,x3)
        y3 = self._.sub(mu,y3)
        return Point(x3,y3)

    def PointAdd(self,a,b):
        if a == self.O:
            return b
        if b == self.O:
            return a
        if a == self.ec_inv(b):
            return self.O
        return self.__PointAdd_BrierJoye(a,b)


class ECDSA_SignatureGenerator(object):
    def __init__(self, group_operation, base_point, bp_order):
        self._ = group_operation
        self.modulus = group_operation.modulus()
        self.bp = base_point
        self.bp_order = bp_order
        assert(self._.IsOnCurve(base_point))
        assert(self._.IsNeutralPoint(self._.ScalarMul(bp_order,base_point)))

    def ComputePublicKey(self,private_key):
        assert(private_key>0)
        assert(private_key<self.bp_order)
        public_key = self._.ScalarMul(private_key,self.bp)
        assert(not self._.IsNeutralPoint(public_key))
        assert(self._.IsOnCurve(public_key))
        assert(self._.IsNeutralPoint(self._.ScalarMul(self.bp_order,public_key)))
        return public_key

    @staticmethod
    def int_to_bytes(i):
        l = (i.bit_length() + 7)//8
        return i.to_bytes(l,byteorder='big')

    @staticmethod
    def bytes_to_int(b):
        return int.from_bytes(b,byteorder='big')

    @staticmethod
    def bytes_bitlen(b):
        i = ECDSA_SignatureGenerator.bytes_to_int(b)
        return i.bit_length()

    def generate_k(self,h1,private_key):
        #deterministic generation of k according to RFC 6979:
        #V = 0x010101...01, same length as h1
        #k = HMAC_K(V || 0x00 || int2octets(x) || bits2octets(h1))
        hlen = 32
        qlen = self.modulus.bit_length()
        V = b'\x01'*hlen
        x = self.int_to_bytes(private_key)
        K = b'\x00'*hlen
        mac = hmac.new(K,digestmod=sha256)
        mac.update(V)
        mac.update(b'\x00')
        mac.update(x)
        mac.update(h1)
        K = mac.digest()
        #e: V = HMAC_K(V)
        V = hmac.new(K,V,digestmod=sha256).digest()
        #f: K = HMAC_K(V || 0x01 || int2octets(x) || bits2octets(h1))
        mac = hmac.new(K,digestmod=sha256)
        mac.update(V)
        mac.update(b'\x01')
        mac.update(x)
        mac.update(h1)
        K = mac.digest()
        #g: V = HMAC_K(V)
        V = hmac.new(K,V,digestmod=sha256).digest()
        k=0
        bitmask = (1<<qlen)-1
        while (k==0) | (k>=self.modulus):
            T = b''
            while self.bytes_bitlen(T)<qlen:
                V = hmac.new(K,V,digestmod=sha256).digest()
                T = T+V

            #not specified but seems necessary to avoid infinte loop
            k = self.bytes_to_int(T)
            k = k & bitmask
            #print("k=",k)

            #K = HMAC_K(V || 0x00)
            K = hmac.new(K,V+b'\x00',digestmod=sha256).digest()
            #V = HMAC_K(V)
            V = hmac.new(K,V,digestmod=sha256).digest()
        return k

    def sign(self, message, private_key):
        calc = ModArithmetic(self.bp_order)
        ln = self.bp_order.bit_length()
        e = sha256(message).digest()
        z = self.bytes_to_int(e[0:ln])
        r=0
        s=0
        while s==0:
            while r==0:
                k = self.generate_k(e,private_key)
                print("k=%x"%k)
                kG = self._.ScalarMul(k,self.bp)
                r = kG.x % self.bp_order
            s = calc.mul(r,private_key)
            s = calc.add(s,z)
            s = calc.div(s,k)
        print("k=%x"%k)
        return (r,s)

class ECDSA_SignatureVerifier(object):
    def __init__(self, group_operation, base_point, bp_order):
        self._ = group_operation
        self.modulus = group_operation.modulus()
        self.bp = base_point
        self.bp_order = bp_order
        assert(self._.IsOnCurve(base_point))
        assert(self._.IsNeutralPoint(self._.ScalarMul(bp_order,base_point)))

    def verify(self,message, signature, public_key):
        calc = ModArithmetic(self.bp_order)
        (r,s) = signature
        assert(r>0)
        assert(r<self.bp_order)
        assert(s>0)
        assert(s<self.bp_order)
        assert(not self._.IsNeutralPoint(public_key))
        assert(self._.IsOnCurve(public_key))
        assert(self._.IsNeutralPoint(self._.ScalarMul(self.bp_order,public_key)))
        e = sha256(message).digest()
        ln = self.bp_order.bit_length()
        z = ECDSA_SignatureGenerator.bytes_to_int(e[0:ln])
        w = calc.inv(s)
        u1 = calc.mul(z,w)
        u2 = calc.mul(r,w)
        u1G = self._.ScalarMul(u1,self.bp)
        u2Qa = self._.ScalarMul(u2,public_key)
        p = self._.PointAdd(u1G,u2Qa)
        out = calc.sub(p.x,r)
        if out != 0:
            raise(Exception('signature verification failed'))



def TestGroupOperation(CalcClass):
    calc = CalcClass(2,2,17)
    p = Point(5,1)
    acc = calc.O
    for i in range(1,23):
        acc = calc.PointAdd(acc,p)
        print("%2d*p="%i,acc,calc.ScalarMul(i,p))
        assert(acc == calc.ScalarMul(i,p))
        assert(calc.IsOnCurve(acc))

def RFC6979_testvector(CalcClass):
    #q = 0x4000000000000000000020108A2E0CC0D99F8A5EF
    #x = 0x09A4D6792295A7F730FC3F2B49CBC0F62E862272F
    #Ux = 0x79AEE090DB05EC252D5CB4452F356BE198A4FF96F
    #Uy = 0x782E29634DDC9A31EF40386E896BAA18B53AFA5A3
    ##Curve K163
    #a=1
    #n=5846006549323611672814741753598448348329118574063

    #A.2.3.  ECDSA, 192 Bits (Prime Field)
    #Key pair:
    #curve: NIST P-192
    q = 0xFFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831
    x = 0x6FAB034934E4C0FC9AE67F5B5659A9D7D1FEFD187EE09FD4
    U = Point(
        0xAC2C77F529F91689FEA0EA5EFEC7F210D8EEA0B9E047ED56,
        0x3BC723E57670BD4887EBC732C523063D0A7C957BC97C1C43)
    k = 0x32B1B6D7D42A05CB449065727A84804FB1A3E34D8F261496
    r = 0x4B0B8CE98A92866A2820E20AA6B75B56382E0F9BFD5ECB55
    s = 0xCCDB006926EA9565CBADC840829D8C384E06DE1F1E381B85
    expected_sig = (r,s)

    n = 6277101735386680763835789423176059013767194773182842284081
    modulus = 6277101735386680763835789423207666416083908700390324961279
    a = -3
    b = 0x64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1
    G = Point(
        0x188da80eb03090f67cbf20eb43a18800f4ff0afd82ff1012,
        0x07192b95ffc8da78631011ed6b24cdd573f977a11e794811)

    signer = ECDSA_SignatureGenerator(CalcClass(a,b,modulus),G,n)
    private_key = x
    public_key = signer.ComputePublicKey(private_key)
    print("public_key=",public_key)
    assert(U==public_key)
    message = bytes("sample", 'utf-8')
    sig = signer.sign(message,private_key)
    print("signature=",sig)
    assert(sig==expected_sig)

def ECDSA_basictest(CalcClass):
    signer = ECDSA_SignatureGenerator(CalcClass(2,2,17),Point(5,1),19)
    private_key = 12
    public_key = signer.ComputePublicKey(private_key)
    print("public_key=",public_key)
    message = bytes("hello world!" , 'utf-8')
    sig = signer.sign(message,private_key)
    print("signature=",sig)
    verifier = ECDSA_SignatureVerifier(CalcClass(2,2,17),Point(5,1),19)
    try:
        message = bytes("hellp world!" , 'utf-8')
        verifier.verify(message,sig,public_key)
        raise(Exception("message corruption not detected"))
    except:
        print("message corruption detected")
    message = bytes("hello world!" , 'utf-8')
    try:
        sig2 = sig
        sig2[0] = sig[0] + 1
        verifier.verify(message,sig2,public_key)
        raise(Exception("signature corruption not detected"))
    except:
        print("signature corruption detected")
    verifier.verify(message,sig,public_key)
    print("correct message verified")

if __name__ == "__main__":
    calc = GroupOperation
    TestGroupOperation(calc)
    ECDSA_basictest(calc)
    RFC6979_testvector(calc)
