#!/usr/bin/env python3

from hashlib import sha256
import hmac
import binascii

from collections import namedtuple
Point = namedtuple("Point", "x y")
#E : y 2 + a 1 xy + a 3 y = x 3 + a 2 x 2 + a 4 x + a 6
Curve = namedtuple("Curve", "a1 a3 a2 a4 a6")

def CurveToString(curve):
    out = "y2"
    if curve.a1:
        if curve.a1>0:
            out+="+"
        out+="%dxy"%curve.a1
    if curve.a3:
        if curve.a3>0:
            out+="+"
        out+="%dy"%curve.a3
    out+=" = x3"
    if curve.a2:
        if curve.a2>0:
            out+="+"
        out+="%dx2"%curve.a2
    if curve.a4:
        if curve.a4>0:
            out+="+"
        out+="%dx"%curve.a4
    if curve.a6:
        if curve.a6>0:
            out+="+"
        out+="%d"%curve.a6
    return out


class ModArithmetic(object):
    def __init__(self, modulus):
        self.__modulus = modulus

    def modulus(self):
        return self.__modulus

    def add(self, a, b):
        return (a+b)%self.__modulus

    def sub(self, a, b):
        o = a-b
        if o<0:
            o+=self.__modulus
        return o % self.__modulus

    def square(self, a):
        return self.mul(a,a)

    def mul(self, a, b):
        return (a*b)%self.__modulus

    def div(self, a, b):
        return (a*self.inv(b))%self.__modulus

    def inv(self, a):
        return ModArithmetic.modinv(a,self.__modulus)

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

class Mod2Arithmetic(object):
    maxshift =0
    sumshift =0
    cntshift =0
    def __init__(self, irreducible_poly):
        self.irreducible_poly = irreducible_poly
        self.msb = 1<<(irreducible_poly.bit_length()-1)

    def modulus(self):
        return self.irreducible_poly

    def add(self, a, b):
        return a^b

    def sub(self, a, b):
        return a^b

    def div(self,a,b):
        return self.mul(a,self.inv(b))

    @staticmethod
    def plain_div(a, b):
        assert(b>0)
        if a<b:
            return (0,a)
        if b==1:
            return (a,0)
        divider = b
        divMsb = divider.bit_length()-1
        remainder=a
        remMsb = remainder.bit_length()-1
        quotient = 0
        maxshift=Mod2Arithmetic.maxshift
        while remMsb>=divMsb:
            shift = remMsb - divMsb
            maxshift = max(maxshift,shift)
            Mod2Arithmetic.cntshift+=1
            Mod2Arithmetic.sumshift+=shift
            quotient ^= 1 << shift
            remainder = remainder ^ (divider << shift)
            remMsb = remainder.bit_length()-1
        if maxshift>Mod2Arithmetic.maxshift:
            print("maxshift=",maxshift)
            Mod2Arithmetic.maxshift = maxshift
            print("cntshift=",Mod2Arithmetic.cntshift)
            print("average shift=",Mod2Arithmetic.sumshift/Mod2Arithmetic.cntshift)
        return (quotient, remainder)

    @staticmethod
    def mod(a, m):
        (quotient, remainder) = Mod2Arithmetic.plain_div(a, b)
        return remainder

    def plain_mul(self, a, b):
        o = 0
        for i in range(0,a.bit_length()):
            if(1 & (a>>i)):
                o ^= b<<i
        return o

    def mul1(self, a, b):
        o = self.plain_mul(a,b)
        return self.mod(o,self.irreducible_poly)

    def square(self, a):
        return self.mul(a,a)

    def mul(self, a, b):
        p = 0
        for i in range(a.bit_length()-1,-1,-1):
            p=p<<1
            if a & (1<<i):
                p ^= b
            if p & self.msb:
                p ^= self.irreducible_poly
        return p

    def inv(self, a):
        (d,s,t) = self.egcd(a,self.irreducible_poly)
        return s


    def egcd(self,g, h):
        #print("g=%x"%g)
        #print("h=%x"%h)
        #print(type(g))
        #print(type(h))
        if h == 0:
            return (g, 1, 0)
        s2=1
        s1=0
        t2=0
        t1=1
        while h:
            (q,r) = self.plain_div(g,h)
            s = self.sub(s2,self.plain_mul(q,s1))
            t = self.sub(t2,self.plain_mul(q,t1))
            g = h
            h = r
            s2 = s1
            s1 = s
            t2 = t1
            t1 = t
        d = g
        s = s2
        t = t2
        return (d,s,t)




class GroupOperation(object):
    O = 'Origin'

    def __init__(self, curve,modulus,arithmetic = ModArithmetic):
        #E : y 2 + a 1 xy + a 3 y = x 3 + a 2 x 2 + a 4 x + a 6
        self.curve=curve
        #print(CurveToString(curve))
        self.arithmetic = arithmetic
        self._ = arithmetic(modulus)

    def modulus(self):
        return self._.modulus()

    def IsNeutralPoint(self,p):
        return p==self.O

    def __PointDouble(self,a):
        if a == self.O:
            return self.O
        if self.arithmetic==Mod2Arithmetic:
            assert(self.curve.a1==1)
            assert(self.curve.a3==0)
            assert(self.curve.a4==0)
            s = self._.div(a.y,a.x)
            l = self._.add(a.x,s)
            ox = self._.mul(l,l)
            ox = self._.add(ox,l)
            ox = self._.add(ox,self.curve.a2)
            oy = self._.mul(a.x,a.x)
            tmp = self._.mul(l,ox)
            oy = self._.add(oy,tmp)
            oy = self._.add(oy,ox)
        else:
            x = a.x
            y = a.y
            x2 = self._.mul(x,x)
            s = self._.mul(3,x2)
            tmp = self._.mul(2,self.curve.a2)
            tmp = self._.mul(tmp,x)
            s = self._.add(s,tmp)
            s = self._.add(s,self.curve.a4)
            tmp = self._.mul(y,self.curve.a1)
            s = self._.sub(s,tmp)
            d = self._.mul(2,y)
            tmp = self._.mul(x,self.curve.a1)
            d = self._.add(d,tmp)
            d = self._.add(d,self.curve.a3)
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
        if self.arithmetic==Mod2Arithmetic:
            assert(self.curve.a1==1)
            assert(self.curve.a3==0)
            assert(self.curve.a4==0)
            n = self._.sub(a.y,b.y)
            d = self._.sub(a.x,b.x)
            l = self._.div(n,d)
            #print(a)
            #print(b)
            #print("n=%x"%n)
            #print("d=%x"%d)
            #print("l=%x"%l)
            #print("self.curve.a2=%x"%self.curve.a2)
            ox = self._.mul(l,l)
            ox = self._.add(ox,l)
            ox = self._.add(ox,self.curve.a2)
            ox = self._.add(ox,a.x)
            ox = self._.add(ox,b.x)
            tmp = self._.add(a.x,ox)
            oy = self._.mul(l,tmp)
            oy = self._.add(oy,ox)
            oy = self._.add(oy,a.y)
        else:
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

    #def ScalarMulSec(self,scalar,p):
    #    rand1 = 0x123456789
    #    temp = self._.mul(r,d)
    #    return self.ScalarMul(temp,p)

    def ShamirsTrick(self,a,p,b,q):
        """returns ap+bq"""
        pq = self.PointAdd(p,q)
        t = self.O
        l = max(a.bit_length(),b.bit_length())
        for i in range(l-1,-1,-1):
            t = self.__PointDouble(t)
            ai = (a & (1<<i))
            bi = (b & (1<<i))
            if ai & bi:
                t = self.PointAdd(t,pq)
            elif ai:
                t = self.PointAdd(t,p)
            elif bi:
                t = self.PointAdd(t,q)
        return t

    def IsOnCurve(self,p):
        if p == self.O:
            return True
        #print("p=",p)
        x2 = self._.mul(p.x,p.x)
        x3 = self._.mul(p.x,x2)
        a2x2 = self._.mul(self.curve.a2,x2)
        right = self._.mul(self.curve.a4,p.x)
        right = self._.add(right,self.curve.a6)
        right = self._.add(right,x3)
        right = self._.add(right,a2x2)
        y2 = self._.mul(p.y,p.y)
        a1xy = self._.mul(self.curve.a1,p.x)
        a1xy = self._.mul(a1xy,p.y)
        a3y = self._.mul(self.curve.a3,p.y)
        left = self._.add(y2,a1xy)
        left = self._.add(left,a3y)
        #print("right=",right)
        #print("left =",left)
        return left==right

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
        self.deterministicECDSA=True
        self._ = group_operation
        self.modulus = group_operation.modulus()
        self.bp = base_point
        self.bp_order = bp_order
        self.q = bp_order
        self.qlen = self.bp_order.bit_length()
        self.rlen = ((self.qlen+7)//8)*8
        assert(self._.IsOnCurve(base_point))
        if self._.arithmetic != Mod2Arithmetic:
            assert(self._.IsNeutralPoint(self._.ScalarMul(bp_order,base_point)))
        else:
            print("base_point:",self.bp)
            print("neutral ?",self._.ScalarMul(bp_order,base_point))

    def ComputePublicKey(self,private_key):
        assert(private_key>0)
        assert(private_key<self.bp_order)
        public_key = self._.ScalarMul(private_key,self.bp)
        assert(not self._.IsNeutralPoint(public_key))
        assert(self._.IsOnCurve(public_key))
        if self._.arithmetic != Mod2Arithmetic:
            assert(self._.IsNeutralPoint(self._.ScalarMul(self.bp_order,public_key)))
        else:
            print("base_point:",self.bp)
            print("neutral ?",self._.ScalarMul(self.bp_order,public_key))
        return public_key

    @staticmethod
    def int_to_bytes(i,bitlen=None):
        if bitlen is None:
            bitlen = (i.bit_length() + 7)//8
        i = i & ((1<<bitlen)-1)
        return i.to_bytes((bitlen+7)//8,byteorder='big')

    @staticmethod
    def bytes_to_int(b):
        return int.from_bytes(b,byteorder='big')

    def bits2int(self,b):
        v = int.from_bytes(b,byteorder='big')
        vlen = len(b)*8
        if(vlen > self.qlen):
            v = v>>(vlen-self.qlen)
        return v

    def bits2octets(self,b):
        z2 = self.bits2int(b)
        if z2 >= self.q:
            z2 = z2-self.q
        return ECDSA_SignatureGenerator.int_to_bytes(z2,self.rlen)

    @staticmethod
    def bytes_bitlen(b):
        i = ECDSA_SignatureGenerator.bytes_to_int(b)
        return i.bit_length()

    def generate_k(self,h1,private_key):
        #deterministic generation of k according to RFC 6979:
        #V = 0x010101...01, same length as h1
        #k = HMAC_K(V || 0x00 || int2octets(x) || bits2octets(h1))
        hlen = 32
        #print("h1=",binascii.hexlify(h1))
        h1 = self.bits2octets(h1)
        #print("h1=",binascii.hexlify(h1))
        V = b'\x01'*hlen
        #print("V=",binascii.hexlify(V))
        x = self.int_to_bytes(private_key,self.rlen)
        #print("x=",binascii.hexlify(x))
        K = b'\x00'*hlen
        mac = hmac.new(K,digestmod=sha256)
        mac.update(V)
        mac.update(b'\x00')
        mac.update(x)
        mac.update(h1)
        K = mac.digest()
        #print("K=",binascii.hexlify(K))
        #e: V = HMAC_K(V)
        V = hmac.new(K,V,digestmod=sha256).digest()
        #print("V=",binascii.hexlify(V))

        #f: K = HMAC_K(V || 0x01 || int2octets(x) || bits2octets(h1))
        mac = hmac.new(K,digestmod=sha256)
        mac.update(V)
        mac.update(b'\x01')
        mac.update(x)
        mac.update(h1)
        K = mac.digest()
        #print("K=",binascii.hexlify(K))

        #g: V = HMAC_K(V)
        V = hmac.new(K,V,digestmod=sha256).digest()
        #print("V=",binascii.hexlify(V))
        k=0
        bitmask = (1<<self.qlen)-1
        while (k==0) | (k>=self.q):
            T = b''
            while self.bytes_bitlen(T)<self.qlen:
                V = hmac.new(K,V,digestmod=sha256).digest()
                T = T+V
                #print("T=",binascii.hexlify(T))

            k = self.bits2int(T)
            #print("k=%x"%k)

            #K = HMAC_K(V || 0x00)
            K = hmac.new(K,V+b'\x00',digestmod=sha256).digest()
            #print("K=",binascii.hexlify(K))
            #V = HMAC_K(V)
            V = hmac.new(K,V,digestmod=sha256).digest()
            #print("V=",binascii.hexlify(V))

        return k

    def sign(self, message, private_key):
        calc = ModArithmetic(self.bp_order)
        ln = self.bp_order.bit_length()
        e = sha256(message).digest()
        if self.deterministicECDSA:
            z = self.bits2int(e)
        else:
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
            #side channel counter measure, of course a fresh random shall be taken each time rather than a constant
            rnd = 0x12345678
            s = calc.mul(s,rnd)
            k = calc.mul(k,rnd)
            s = calc.div(s,k)
        return (r,s)

class ECDSA_SignatureVerifier(object):
    def __init__(self, group_operation, base_point, bp_order):
        self.deterministicECDSA=True
        self._ = group_operation
        self.modulus = group_operation.modulus()
        self.bp = base_point
        self.bp_order = bp_order
        self.q = bp_order
        self.qlen = self.bp_order.bit_length()
        self.rlen = ((self.qlen+7)//8)*8
        assert(self._.IsOnCurve(base_point))
        if self._.arithmetic != Mod2Arithmetic:
            assert(self._.IsNeutralPoint(self._.ScalarMul(bp_order,base_point)))

    def bits2int(self,b):
        v = int.from_bytes(b,byteorder='big')
        vlen = len(b)*8
        if(vlen > self.qlen):
            v = v>>(vlen-self.qlen)
        return v

    def verify(self,message, signature, public_key):
        calc = ModArithmetic(self.bp_order)
        (r,s) = signature
        assert(r>0)
        assert(r<self.bp_order)
        assert(s>0)
        assert(s<self.bp_order)
        assert(not self._.IsNeutralPoint(public_key))
        assert(self._.IsOnCurve(public_key))
        if self._.arithmetic != Mod2Arithmetic:
            assert(self._.IsNeutralPoint(self._.ScalarMul(self.bp_order,public_key)))
        e = sha256(message).digest()
        ln = self.bp_order.bit_length()
        #z = ECDSA_SignatureGenerator.bytes_to_int(e[0:ln])
        if self.deterministicECDSA:
            z = self.bits2int(e)
        else:
            z = ECDSA_SignatureGenerator.bytes_to_int(e[0:ln])
        w = calc.inv(s)
        u1 = calc.mul(z,w)
        u2 = calc.mul(r,w)
        #u1G = self._.ScalarMul(u1,self.bp)
        #u2Qa = self._.ScalarMul(u2,public_key)
        #p = self._.PointAdd(u1G,u2Qa)
        p = self._.ShamirsTrick(u1,self.bp,u2,public_key)
        out = calc.sub(p.x,r)
        if out != 0:
            raise(Exception('signature verification failed'))



def testMod2ModArithmetic():
    f = 0x11B
    a = 0x53
    b = 0xCA
    calc = Mod2Arithmetic(f)
    s = calc.mul(a,b)
    print("%x"%s)
    assert(s == 0x01)
    s = calc.mul(0x02,0x87)
    print("%x"%s)
    assert(s == 0x15)

    g=(1<<10)|(1<<9)|(1<<8)|(1<<6)|(1<<5)|(1<<4)|1
    h=(1<<9)|(1<<6)|(1<<5)|(1<<3)|(1<<2)|1
    (d,s,t) = calc.egcd(g,h)
    print("%x"%d)
    print("%x"%s)
    print("%x"%t)
    assert(d==((1<<3)|(1<<1)|1))
    assert(s==((1<<4)))
    assert(t==((1<<5)|(1<<4)|(1<<3)|(1<<2)|(1<<1)|1))

    assert(0xCA==calc.inv(0x53))


def TestGroupOperationMod2():
    g = 2
    f = (1<<4)|(1<<1)|1
    acc=1
    calc = Mod2Arithmetic(f)
    gp = [1]
    for i in range(1,16):
        acc = calc.mul(acc,g)
        gp.append(acc)
        print("%x"%acc)
    assert(gp[4]==3)
    assert(gp[15]==1)
    assert(gp[3]==8)
    assert(gp[5]==6)
    print("%x"%calc.div(gp[3],gp[5]))
    assert(gp[13]==calc.div(gp[3],gp[5]))
    assert(gp[7]==calc.add(gp[5],gp[13]))
    assert(gp[0] ==calc.add(calc.mul(gp[7],gp[7]),calc.add(gp[7],gp[4])))
    assert(gp[13]==calc.add(calc.mul(gp[5],gp[5]),calc.add(gp[7],gp[0])))
    a = gp[4]
    curve = Curve(1,0,a,0,1)
    calc = GroupOperation(curve,f,Mod2Arithmetic)
    p = Point(gp[5],gp[3])
    r = calc.PointAdd(p,p)
    print(r)
    assert(r.x==gp[0])
    assert(r.y==gp[13])
    q = Point(gp[9],gp[13])
    r = calc.PointAdd(p,q)
    print(r)
    assert(r.x==gp[3])
    assert(r.y==gp[13])


def TestGroupOperation(CalcClass):
    curve = Curve(0,0,0,2,2)
    calc = CalcClass(curve,17)
    p = Point(5,1)
    acc = calc.O
    for i in range(1,23):
        acc = calc.PointAdd(acc,p)
        print("%2d*p="%i,acc,calc.ScalarMul(i,p))
        assert(acc == calc.ScalarMul(i,p))
        assert(calc.IsOnCurve(acc))

def RFC6979_testvector(CalcClass):
    message = bytes("sample", 'utf-8')

    #Curve B163 (FIPS 186-4)
    test_Fp_curve(
        calc_class=CalcClass,
        curve=Curve(1,0,1,0,0x20a601907b8c953ca1481eb10512f78744a3205fd),
        modulus=(1<<163)|(1<<7)|(1<<6)|(1<<3)|1, #q in RFC6979_testvector
        arithmetic=Mod2Arithmetic,
        base_point=Point(
            0x3f0eba16286a2d57ea0991168d4994637e8343e36,
            0x0d51fbc6c71a0094fa2cdd545b11c5c0c797324f1),
        base_point_order=5846006549323611672814742442876390689256843201587,
        private_key=0x35318FC447D48D7E6BC93B48617DDDEDF26AA658F,
        public_key=Point(
            0x126CF562D95A1D77D387BA75A3EA3A1407F23425A,
            0x7D7CB5273C94DA8CA93049AFDA18721C24672BD71),
        message=message,
        sig=(0x134E00F78FC1CB9501675D91C401DE20DDF228CDC,
             0x373273AEC6C36CB7BAFBB1903A5F5EA6A1D50B624))

    #Curve B233 (FIPS 186-4)
    test_Fp_curve(
        calc_class=CalcClass,
        curve=Curve(1,0,1,0,0x066647ede6c332c7f8c0923bb58213b333b20e9ce4281fe115f7d8f90ad),
        modulus=(1<<233)|(1<<74)|1, #q in RFC6979_testvector
        arithmetic=Mod2Arithmetic,
        base_point=Point(
            0x0fac9dfcbac8313bb2139f1bb755fef65bc391f8b36f8f8eb7371fd558b,
            0x1006a08a41903350678e58528bebf8a0beff867a7ca36716f7e01f81052),
        base_point_order=6901746346790563787434755862277025555839812737345013555379383634485463 ,
        private_key=0x07ADC13DD5BF34D1DDEEB50B2CE23B5F5E6D18067306D60C5F6FF11E5D3,
        public_key=Point(
            0x0FB348B3246B473AA7FBB2A01B78D61B62C4221D0F9AB55FC72DB3DF478,
            0x1162FA1F6C6ACF7FD8D19FC7D74BDD9104076E833898BC4C042A6E6BEBF),
        message=message,
        sig=(0x0A797F3B8AEFCE7456202DF1E46CCC291EA5A49DA3D4BDDA9A4B62D5E0D,
             0x01F6F81DA55C22DA4152134C661588F4BD6F82FDBAF0C5877096B070DC2))

    #curve: NIST B-283
    test_Fp_curve(
        calc_class=CalcClass,
        curve=Curve(1,0,1,0,0x27b680ac8b8596da5a4af8a19a0303fca97fd7645309fa2a581485af6263e313b79a2f5),
        modulus=(1<<283)|(1<<12)|(1<<7)|(1<<5)|1, #q in RFC6979_testvector
        arithmetic=Mod2Arithmetic,
        base_point=Point(
            0x5f939258db7dd90e1934f8c70b0dfec2eed25b8557eac9c80e2e198f8cdbecd86b12053,
            0x3676854fe24141cb98fe6d4b20d02b4516ff702350eddb0826779c813f0df45be8112f4),
        base_point_order=7770675568902916283677847627294075626569625924376904889109196526770044277787378692871,
        private_key=0x14510D4BC44F2D26F4553942C98073C1BD35545CEABB5CC138853C5158D2729EA408836,
        public_key=Point(
            0x17E3409A13C399F0CA8A192F028D46E3446BCFFCDF51FF8A905ED2DED786E74F9C3E8A9,
            0x47EFCBCC31C01D86D1992F7BFAC0277DBD02A6D289274099A2C0F039C8F59F318371B0E),
        message=message,
        sig=(0x29FD82497FB3E5CEF65579272138DE59E2B666B8689466572B3B69A172CEE83BE145659,
             0x05A89D9166B40795AF0FE5958201B9C0523E500013CA12B4840EA2BC53F25F9B3CE87C0))

    #Curve K163 (FIPS 186-4)
    test_Fp_curve(
        calc_class=CalcClass,
        curve=Curve(1,0,1,0,1),
        modulus=(1<<163)|(1<<7)|(1<<6)|(1<<3)|1, #q in RFC6979_testvector
        arithmetic=Mod2Arithmetic,
        base_point=Point(
            0x2fe13c0537bbc11acaa07d793de4e6d5e5c94eee8,
            0x289070fb05d38ff58321f2e800536d538ccdaa3d9),
        base_point_order=5846006549323611672814741753598448348329118574063,
        private_key=0x09A4D6792295A7F730FC3F2B49CBC0F62E862272F,
        public_key=Point(
            0x79AEE090DB05EC252D5CB4452F356BE198A4FF96F,
            0x782E29634DDC9A31EF40386E896BAA18B53AFA5A3),
        message=message,
        sig=(0x113A63990598A3828C407C0F4D2438D990DF99A7F,
             0x1313A2E03F5412DDB296A22E2C455335545672D9F))

    #curve: NIST K-233
    test_Fp_curve(
        calc_class=CalcClass,
        curve=Curve(1,0,0,0,1),
        modulus=(1<<233)|(1<<74)|1, #q in RFC6979_testvector
        arithmetic=Mod2Arithmetic,
        base_point=Point(
            0x17232ba853a7e731af129f22ff4149563a419c26bf50a4c9d6eefad6126,
            0x1db537dece819b7f70f555a67c427a8cd9bf18aeb9b56e0c11056fae6a3),
        base_point_order=3450873173395281893717377931138512760570940988862252126328087024741343,
        private_key=0x103B2142BDC2A3C3B55080D09DF1808F79336DA2399F5CA7171D1BE9B0,
        public_key=Point(
            0x0682886F36C68473C1A221720C2B12B9BE13458BA907E1C4736595779F2,
            0x1B20639B41BE0927090999B7817A3B3928D20503A39546044EC13A10309),
        message=message,
        sig=(0x38AD9C1D2CB29906E7D63C24601AC55736B438FB14F4093D6C32F63A10,
             0x647AAD2599C21B6EE89BE7FF957D98F684B7921DE1FD3CC82C079624F4))

    #curve: NIST K-283
    test_Fp_curve(
        calc_class=CalcClass,
        curve=Curve(1,0,0,0,1),
        modulus=(1<<283)|(1<<12)|(1<<7)|(1<<5)|1, #q in RFC6979_testvector
        arithmetic=Mod2Arithmetic,
        base_point=Point(
            0x503213f78ca44883f1a3b8162f188e553cd265f23c1567a16876913b0c2ac2458492836,
            0x1ccda380f1c9e318d90f95d07e5426fe87e45c0e8184698e45962364e34116177dd2259),
        base_point_order=3885337784451458141838923813647037813284811733793061324295874997529815829704422603873  ,
        private_key=0x06A0777356E87B89BA1ED3A3D845357BE332173C8F7A65BDC7DB4FAB3C4CC79ACC8194E,
        public_key=Point(
            0x25330D0A651D5A20DC6389BC02345117725640AEC3C126612CE444EDD19649BDECC03D6,
            0x505BD60A4B67182474EC4D1C668A73140F70504A68F39EFCD972487E9530E0508A76193),
        message=message,
        sig=(0x19E90AA3DE5FB20AED22879F92C6FED278D9C9B9293CC5E94922CD952C9DBF20DF1753A,
             0x135AA7443B6A25D11BB64AC482E04D47902D017752882BD72527114F46CF8BB56C5A8C3))

    #curve: NIST K-409
    test_Fp_curve(
        calc_class=CalcClass,
        curve=Curve(1,0,0,0,1),
        modulus=(1<<409)|(1<<87)|1, #q in RFC6979_testvector
        arithmetic=Mod2Arithmetic,
        base_point=Point(
            0x060f05f658f49c1ad3ab1890f7184210efd0987e307c84c27accfb8f9f67cc2c460189eb5aaaa62ee222eb1b35540cfe9023746,
            0x1e369050b7c4e42acba1dacbf04299c3460782f918ea427e6325165e9ea10e3da5f6c42e9c55215aa9ca27a5863ec48d8e0286b),
        base_point_order=330527984395124299475957654016385519914202341482140609642324395022880711289249191050673258457777458014096366590617731358671,
        private_key=0x29C16768F01D1B8A89FDA85E2EFD73A09558B92A178A2931F359E4D70AD853E569CDAF16DAA569758FB4E73089E4525D8BBFCF,
        public_key=Point(
            0x0CF923F523FE34A6E863D8BA45FB1FE6D784C8F219C414EEF4DB8362DBBD3CA71AEB28F568668D5D7A0093E2B84F6FAD759DB42,
            0x13B1C374D5132978A1B1123EBBE9A5C54D1A9D56B09AFDB4ADE93CCD7C4D332E2916F7D4B9D18578EE3C2E2DE4D2ECE0DE63549),
        message=message,
        sig=(0x49EC220D6D24980693E6D33B191532EAB4C5D924E97E305E2C1CCFE6F1EAEF96C17F6EC27D1E06191023615368628A7E0BD6A9,
             0x1A4AB1DD9BAAA21F77C503E1B39E770FFD44718349D54BA4CF08F688CE89D7D7C5F7213F225944BE5F7C9BA42B8BEE382F8AF9))

    #curve: NIST K-571
    test_Fp_curve(
        calc_class=CalcClass,
        curve=Curve(1,0,0,0,1),
        modulus=(1<<571)|(1<<10)|(1<<5)|(1<<2)|1, #q in RFC6979_testvector
        arithmetic=Mod2Arithmetic,
        base_point=Point(
            0x26eb7a859923fbc82189631f8103fe4ac9ca2970012d5d46024804801841ca44370958493b205e647da304db4ceb08cbbd1ba39494776fb988b47174dca88c7e2945283a01c8972,
            0x349dc807f4fbf374f4aeade3bca95314dd58cec9f307a54ffc61efc006d8a2c9d4979c0ac44aea74fbebbb9f772aedcb620b01a7ba7af1b320430c8591984f601cd4c143ef1c7a3),
        base_point_order=1932268761508629172347675945465993672149463664853217499328617625725759571144780212268133978522706711834706712800825351461273674974066617311929682421617092503555733685276673,
        private_key=0x0C16F58550D824ED7B95569D4445375D3A490BC7E0194C41A39DEB732C29396CDF1D66DE02DD1460A816606F3BEC0F32202C7BD18A32D87506466AA92032F1314ED7B19762B0D22,
        public_key=Point(
            0x6CFB0DF7541CDD4C41EF319EA88E849EFC8605D97779148082EC991C463ED32319596F9FDF4779C17CAF20EFD9BEB57E9F4ED55BFC52A2FA15CA23BC62B7BF019DB59793DD77318,
            0x1CFC91102F7759A561BD8D5B51AAAEEC7F40E659D67870361990D6DE29F6B4F7E18AE13BDE5EA5C1F77B23D676F44050C9DBFCCDD7B3756328DDA059779AAE8446FC5158A75C227),
        message=message,
        sig=(0x1604BE98D1A27CEC2D3FA4BD07B42799E07743071E4905D7DCE7F6992B21A27F14F55D0FE5A7810DF65CF07F2F2554658817E5A88D952282EA1B8310514C0B40FFF46F159965168,
             0x18249377C654B8588475510F7B797081F68C2F8CCCE49F730353B2DA3364B1CD3E984813E11BB791824038EA367BA74583AB97A69AF2D77FA691AA694E348E15DA76F5A44EC1F40))

    #curve: NIST P-192
    test_Fp_curve(
        calc_class=CalcClass,
        curve=Curve(0,0,0,-3,0x64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1),
        modulus=6277101735386680763835789423207666416083908700390324961279, #q in RFC6979_testvector
        arithmetic=ModArithmetic,
        base_point=Point(
            0x188da80eb03090f67cbf20eb43a18800f4ff0afd82ff1012,
            0x07192b95ffc8da78631011ed6b24cdd573f977a11e794811),
        base_point_order=6277101735386680763835789423176059013767194773182842284081,
        private_key=0x6FAB034934E4C0FC9AE67F5B5659A9D7D1FEFD187EE09FD4,
        public_key=Point(
            0xAC2C77F529F91689FEA0EA5EFEC7F210D8EEA0B9E047ED56,
            0x3BC723E57670BD4887EBC732C523063D0A7C957BC97C1C43),
        message=message,
        sig=(0x4B0B8CE98A92866A2820E20AA6B75B56382E0F9BFD5ECB55,
        0xCCDB006926EA9565CBADC840829D8C384E06DE1F1E381B85))

    #curve: NIST P-224
    test_Fp_curve(
        calc_class=CalcClass,
        curve=Curve(0,0,0,-3,0xb4050a850c04b3abf54132565044b0b7d7bfd8ba270b39432355ffb4),
        modulus=26959946667150639794667015087019630673557916260026308143510066298881, #q in RFC6979_testvector
        arithmetic=ModArithmetic,
        base_point=Point(
            0xb70e0cbd6bb4bf7f321390b94a03c1d356c21122343280d6115c1d21,
            0xbd376388b5f723fb4c22dfe6cd4375a05a07476444d5819985007e34),
        base_point_order=26959946667150639794667015087019625940457807714424391721682722368061,
        private_key=0xF220266E1105BFE3083E03EC7A3A654651F45E37167E88600BF257C1,
        public_key=Point(
            0x00CF08DA5AD719E42707FA431292DEA11244D64FC51610D94B130D6C,
            0xEEAB6F3DEBE455E3DBF85416F7030CBD94F34F2D6F232C69F3C1385A),
        message=message,
        sig=(0x61AA3DA010E8E8406C656BC477A7A7189895E7E840CDFE8FF42307BA,
             0xBC814050DAB5D23770879494F9E0A680DC1AF7161991BDE692B10101))

    #curve: NIST P-256
    test_Fp_curve(
        calc_class=CalcClass,
        curve=Curve(0,0,0,-3,0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b),
        modulus=115792089210356248762697446949407573530086143415290314195533631308867097853951, #q in RFC6979_testvector
        arithmetic=ModArithmetic,
        base_point=Point(
            0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296,
            0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5),
        base_point_order=115792089210356248762697446949407573529996955224135760342422259061068512044369,
        private_key=0xC9AFA9D845BA75166B5C215767B1D6934E50C3DB36E89B127B8A622B120F6721,
        public_key=Point(
            0x60FED4BA255A9D31C961EB74C6356D68C049B8923B61FA6CE669622E60F29FB6,
            0x7903FE1008B8BC99A41AE9E95628BC64F2F1B20C2D7E9F5177A3C294D4462299),
        message=message,
        sig=(0xEFD48B2AACB6A8FD1140DD9CD45E81D69D2C877B56AAF991C34D0EA84EAF3716,
             0xF7CB1C942D657C41D436C7A1B6E29F65F3E900DBB9AFF4064DC4AB2F843ACDA8))

     #curve: NIST P-384
    test_Fp_curve(
        calc_class=CalcClass,
        curve=Curve(0,0,0,-3,0xb3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef),
        modulus=39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319 , #q in RFC6979_testvector
        arithmetic=ModArithmetic,
        base_point=Point(
            0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7,
            0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f),
        base_point_order=39402006196394479212279040100143613805079739270465446667946905279627659399113263569398956308152294913554433653942643,
        private_key=0x6B9D3DAD2E1B8C1C05B19875B6659F4DE23C3B667BF297BA9AA47740787137D896D5724E4C70A825F872C9EA60D2EDF5,
        public_key=Point(
            0xEC3A4E415B4E19A4568618029F427FA5DA9A8BC4AE92E02E06AAE5286B300C64DEF8F0EA9055866064A254515480BC13,
            0x8015D9B72D7D57244EA8EF9AC0C621896708A59367F9DFB9F54CA84B3F1C9DB1288B231C3AE0D4FE7344FD2533264720),
        message=message,
        sig=(0x21B13D1E013C7FA1392D03C5F99AF8B30C570C6F98D4EA8E354B63A21D3DAA33BDE1E888E63355D92FA2B3C36D8FB2CD,
             0xF3AA443FB107745BF4BD77CB3891674632068A10CA67E3D45DB2266FA7D1FEEBEFDC63ECCD1AC42EC0CB8668A4FA0AB0))

     #curve: NIST P-521
    test_Fp_curve(
        calc_class=CalcClass,
        curve=Curve(0,0,0,-3,0x051953eb9618e1c9a1f929a21a0b68540eea2da725b99b315f3b8b489918ef109e156193951ec7e937b1652c0bd3bb1bf073573df883d2c34f1ef451fd46b503f00),
        modulus=6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151, #q in RFC6979_testvector
        arithmetic=ModArithmetic,
        base_point=Point(
            0xc6858e06b70404e9cd9e3ecb662395b4429c648139053fb521f828af606b4d3dbaa14b5e77efe75928fe1dc127a2ffa8de3348b3c1856a429bf97e7e31c2e5bd66,
            0x11839296a789a3bc0045c8a5fb42c7d1bd998f54449579b446817afbd17273e662c97ee72995ef42640c550b9013fad0761353c7086a272c24088be94769fd16650),
        base_point_order=6864797660130609714981900799081393217269435300143305409394463459185543183397655394245057746333217197532963996371363321113864768612440380340372808892707005449  ,
        private_key=0x0FAD06DAA62BA3B25D2FB40133DA757205DE67F5BB0018FEE8C86E1B68C7E75CAA896EB32F1F47C70855836A6D16FCC1466F6D8FBEC67DB89EC0C08B0E996B83538,
        public_key=Point(
            0x1894550D0785932E00EAA23B694F213F8C3121F86DC97A04E5A7167DB4E5BCD371123D46E45DB6B5D5370A7F20FB633155D38FFA16D2BD761DCAC474B9A2F5023A4,
            0x0493101C962CD4D2FDDF782285E64584139C2F91B47F87FF82354D6630F746A28A0DB25741B5B34A828008B22ACC23F924FAAFBD4D33F81EA66956DFEAA2BFDFCF5),
        message=message,
        sig=(0x1511BB4D675114FE266FC4372B87682BAECC01D3CC62CF2303C92B3526012659D16876E25C7C1E57648F23B73564D67F61C6F14D527D54972810421E7D87589E1A7,
             0x04A171143A83163D6DF460AAF61522695F207A58B95C0644D87E52AA1A347916E4F7A72930B1BC06DBE22CE3F58264AFD23704CBB63B29B931F7DE6C9D949A7ECFC))


def RFC5639_testvector(CalcClass):
    message = bytes("sample", 'utf-8')

    #curve: brainpoolP160r1
    test_ecdh(
        calc_class=CalcClass,
        curve=Curve(0,0,0,0x340E7BE2A280EB74E2BE61BADA745D97E8F7C300,0x1E589A8595423412134FAA2DBDEC95C8D8675E58),
        modulus=0xE95E4A5F737059DC60DFC7AD95B3D8139515620F, #p in RFC5639_testvector
        arithmetic=ModArithmetic,
        base_point=Point(
            0xBED5AF16EA3F6A4F62938C4631EB5AF7BDBCDBC3,
            0x1667CB477A1A8EC338F94741669C976316DA6321),
        base_point_order=0xE95E4A5F737059DC60DF5991D45029409E60FC09, #q in RFC5639_testvector
        private_key_a=0x8BF7BC5CBE8AC8B34940C2C5652D6AE4EC9F53CE ,
        public_key_a=Point(
            0x651DA24C7FF64DD863F8F650E53F07B8EC943C39,
            0xD1C68C656E44034D0DAD60A1589FD49594E7C2A4),
        private_key_b=0xB6F7160F0DE61CCEAAC528A32BD7AD942E8017B2,
        public_key_b=Point(
            0xDF9F259AEA6DFA1F28B16B8FEC52044CC1DFBA35,
            0xDF72AEA65A5E3EF69166DA161ABE00FC9C81C4D0),
        shared_secret=Point(
            0xD78792AC4CBE3390DDD6557060066BC25579CA97,
            0x3A3DAB50421585FB9DE9D87BB3BBBAFE3379A571)
        )

    #curve: brainpoolP256r1
    test_ecdh(
        calc_class=CalcClass,
        curve=Curve(0,0,0,0x7D5A0975FC2C3057EEF67530417AFFE7FB8055C126DC5C6CE94A4B44F330B5D9,0x26DC5C6CE94A4B44F330B5D9BBD77CBF958416295CF7E1CE6BCCDC18FF8C07B6),
        modulus=0xA9FB57DBA1EEA9BC3E660A909D838D726E3BF623D52620282013481D1F6E5377, #p in RFC5639_testvector
        arithmetic=ModArithmetic,
        base_point=Point(
            0x8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262,
            0x547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997),
        base_point_order=0xA9FB57DBA1EEA9BC3E660A909D838D718C397AA3B561A6F7901E0E82974856A7, #q in RFC5639_testvector
        private_key_a=0x81DB1EE100150FF2EA338D708271BE38300CB54241D79950F77B063039804F1D ,
        public_key_a=Point(
            0x44106E913F92BC02A1705D9953A8414DB95E1AAA49E81D9E85F929A8E3100BE5,
            0x8AB4846F11CACCB73CE49CBDD120F5A900A69FD32C272223F789EF10EB089BDC),
        private_key_b=0x55E40BC41E37E3E2AD25C3C6654511FFA8474A91A0032087593852D3E7D76BD3 ,
        public_key_b=Point(
            0x8D2D688C6CF93E1160AD04CC4429117DC2C41825E1E9FCA0ADDD34E6F1B39F7B,
            0x990C57520812BE512641E47034832106BC7D3E8DD0E4C7F1136D7006547CEC6A),
        shared_secret=Point(
            0x89AFC39D41D3B327814B80940B042590F96556EC91E6AE7939BCE31F3A18BF2B,
            0x49C27868F4ECA2179BFD7D59B1E3BF34C1DBDE61AE12931648F43E59632504DE)
        )



def test_ecdh(calc_class,curve,modulus,arithmetic,base_point,base_point_order,private_key_a,public_key_a,private_key_b,public_key_b,shared_secret):
    print("qlen=",modulus.bit_length())
    group_op = calc_class(curve,modulus,arithmetic)
    signer = ECDSA_SignatureGenerator(group_op,base_point,base_point_order)
    gen_public_key = signer.ComputePublicKey(private_key_a)
    print("\tpublic_key_a=",gen_public_key)
    assert(public_key_a==gen_public_key)
    gen_public_key = signer.ComputePublicKey(private_key_b)
    print("\tpublic_key_b=",gen_public_key)
    assert(public_key_b==gen_public_key)
    secret_a = group_op.ScalarMul(private_key_a,public_key_b)
    print("\tsecret_a=",secret_a)
    secret_b = group_op.ScalarMul(private_key_b,public_key_a)
    print("\tsecret_b=",secret_b)
    assert(secret_b==secret_a)
    print("EDCH successful")

def test_Fp_curve(calc_class,curve,modulus,arithmetic,base_point,base_point_order,private_key,public_key,message,sig):
    print("qlen=",modulus.bit_length())
    signer = ECDSA_SignatureGenerator(calc_class(curve,modulus,arithmetic),base_point,base_point_order)
    gen_public_key = signer.ComputePublicKey(private_key)
    print("\tpublic_key=",gen_public_key)
    assert(public_key==gen_public_key)
    gen_sig = signer.sign(message,private_key)
    print("\tsignature=",gen_sig)
    assert(sig==gen_sig)
    verifier = ECDSA_SignatureVerifier(calc_class(curve,modulus,arithmetic),base_point,base_point_order)
    verifier.verify(message,gen_sig,gen_public_key)

def ECDSA_basictest(CalcClass):
    curve = Curve(0,0,0,2,2)
    signer = ECDSA_SignatureGenerator(CalcClass(curve,17),Point(5,1),19)
    private_key = 12
    public_key = signer.ComputePublicKey(private_key)
    print("public_key=",public_key)
    message = bytes("hello world!" , 'utf-8')
    sig = signer.sign(message,private_key)
    print("signature=",sig)
    verifier = ECDSA_SignatureVerifier(CalcClass(curve,17),Point(5,1),19)
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
    testMod2ModArithmetic()
    TestGroupOperationMod2()
    calc = GroupOperation
    TestGroupOperation(calc)
    ECDSA_basictest(calc)
    RFC6979_testvector(calc)
    RFC5639_testvector(calc)
