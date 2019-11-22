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

    def div1(self,a,b):
        return self.mul(a,self.inv(b))

    @staticmethod
    def even(a):
        return 0==(a & 1)

    @staticmethod
    def odd(a):
        return (a & 1)

    def div2(self,y,x):
        print("div")
        print("y=0x%x"%y)
        print("x=0x%x"%x)
        print("m=0x%x"%self.irreducible_poly)
        assert(y<self.irreducible_poly)
        assert(x<self.irreducible_poly)

        if x==0:
            return y
        a = x
        b = self.irreducible_poly
        u = y
        v = 0
        print("a=%3x, b=%3x, u=%3x, v=%3x"%(a,b,u,v))
        while a != b:
            if self.even(a):
                a = a>>1
                if self.even(u):
                    u = u>>1
                else:
                    u = (u^self.irreducible_poly)>>1
            elif self.even(b):
                b = b>>1
                if self.even(v):
                    v = v>>1
                else:
                    v = (v^self.irreducible_poly)>>1
            elif a>b:
                a = (a^b)>>1
                u = u ^ v
                if self.even(u):
                    u = u>>1
                else:
                    u = (u^self.irreducible_poly)>>1
            else:
                b = (a^b)>>1
                v = u^v
                if self.even(v):
                    v = v>>1
                else:
                    v = (v^self.irreducible_poly)>>1
            print("a=%3x, b=%3x, u=%3x, v=%3x"%(a,b,u,v))
        assert(a==1)
        assert(b==1)
        print("u=0x%x"%u)
        return u

    def div3(self,y,x):
        #print("div")
        #print("y=0x%x"%y)
        #print("x=0x%x"%x)
        #print("m=0x%x"%self.irreducible_poly)
        assert(y<self.irreducible_poly)
        assert(x<self.irreducible_poly)

        if x==0:
            return y
        a = x
        b = self.irreducible_poly
        u = y
        v = 0
        ca = a.bit_length()
        cb = self.irreducible_poly.bit_length()
        ca = cb-1
        #print("a=%3x, b=%3x, u=%3x, v=%3x,ca=%3d, cb=%3d"%(a,b,u,v,ca,cb))
        #while (ca > 1)|(cb > 1):
        while a!=1:
            a_even = self.even(a)
            b_even = self.even(b)
            u_even = self.even(u)
            v_even = self.even(v)
            ca_gt_cb = ca > cb
            if (not a_even) & (not b_even):
                uv_adder = u ^ v
                if ca_gt_cb:
                    a = (a^b)>>1
                    ca -= 1
                    u = uv_adder
                    if self.even(u):
                        u = u>>1
                    else:
                        u = (u^self.irreducible_poly)>>1
                else:
                    b = (a^b)>>1
                    cb -= 1
                    v = uv_adder
                    if self.even(v):
                        v = v>>1
                    else:
                        v = (v^self.irreducible_poly)>>1
            else:
                if a_even:
                    a = a>>1
                    ca -= 1
                    if u_even:
                        u = (u)>>1
                    else:
                        u = (u^self.irreducible_poly)>>1

                if b_even:
                    b = b>>1
                    cb -= 1
                    if v_even:
                        v = (v)>>1
                    else:
                        v = (v^self.irreducible_poly)>>1
            #print("a=%3x, b=%3x, u=%3x, v=%3x,ca=%3d, cb=%3d"%(a,b,u,v,ca,cb))
            assert(ca>0)
            assert(cb>0)
            #print("ca=%3d, cb=%3d"%(ca,cb))
        assert(a==1)
        #assert(b==1)
        #print("u=0x%x"%u)
        return u


    def div(self,y,x):
        #print("div")
        #print("y=0x%x"%y)
        #print("x=0x%x"%x)
        #print("m=0x%x"%self.irreducible_poly)
        assert(y<self.irreducible_poly)
        assert(x<self.irreducible_poly)

        if x==0:
            return y
        a = x
        b = self.irreducible_poly
        u = y
        v = 0
        ca = self.irreducible_poly.bit_length()-1
        cb = self.irreducible_poly.bit_length()
        ca = 4095
        cb = 4096
        assert(self.irreducible_poly.bit_length()<4096)
        #print("a=%3x, b=%3x, u=%3x, v=%3x,ca=%3d, cb=%3d"%(a,b,u,v,ca,cb))
        while a!=1:
            a_even = self.even(a)
            b_even = self.even(b)
            u_even = self.even(u)
            v_even = self.even(v)
            ca_gt_cb = ca > cb
            add_both = (not a_even) & (not b_even)
            update_au = (add_both & ca_gt_cb) | a_even
            update_bv = (add_both & (not ca_gt_cb)) | b_even

            next_uv = 0
            if update_au | add_both:
                next_uv ^= u
            if update_bv | add_both:
                next_uv ^= v
            if self.even(next_uv):
                next_uv = next_uv>>1
            else:
                next_uv = (next_uv^self.irreducible_poly)>>1

            next_ab = 0
            if update_au | add_both:
                next_ab ^= a
            if update_bv | add_both:
                next_ab ^= b
            if self.even(next_ab):
                next_ab = next_ab>>1
            else:
                next_ab = (next_ab^self.irreducible_poly)>>1

            if update_au:
                a = next_ab
                ca -= 1
                u = next_uv
            if update_bv:
                b = next_ab
                cb -= 1
                v = next_uv

            #print("a=%3x, b=%3x, u=%3x, v=%3x,update_au=%d,update_bv=%d,ca=%3d, cb=%3d"%(a,b,u,v,update_au,update_bv,ca,cb))
            assert(ca>0)
            assert(cb>0)
            #print("ca=%3d, cb=%3d"%(ca,cb))
        assert(a==1)
        #assert(b==1)
        #print("u=0x%x"%u)
        return u

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

    def plain_mul1(self, a, b):
        o = 0
        for i in range(0,a.bit_length()):
            if(1 & (a>>i)):
                o ^= b<<i
        return o

    def mul1(self, a, b):
        o = self.plain_mul(a,b)
        return self.mod(o,self.irreducible_poly)

    def square(self, a):
        out = 0
        for i in range(0,a.bit_length()):
            out ^= (a & (1<<i))<<i
        while out.bit_length() >= self.irreducible_poly.bit_length():
            shift = out.bit_length() - self.irreducible_poly.bit_length()
            out ^= self.irreducible_poly << shift
        return out

    plain_mul_s_max = 0
    plain_mul_cnt = 0
    def plain_mul(self, a, b):
        Mod2Arithmetic.plain_mul_s_max = max(Mod2Arithmetic.plain_mul_s_max,a)
        Mod2Arithmetic.plain_mul_cnt += 1
        #print_verilog("plain_mul_s",a)
        #print_verilog("plain_mul_x",b)
        p = 0
        for i in range(a.bit_length()-1,-1,-1):
            p=p<<1
            if a & (1<<i):
                p ^= b
        #print_verilog("plain_mul_o",p)
        return p

    mod_mul_s_max = 0
    mod_mul_cnt = 0
    def mul(self, a, b):
        Mod2Arithmetic.mod_mul_s_max = max(Mod2Arithmetic.mod_mul_s_max,a)
        Mod2Arithmetic.mod_mul_cnt += 1
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

    def __init__(self, curve,modulus,arithmetic = Mod2Arithmetic):
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
            print("a.x=%x"%a.x)
            print("a.y=%x"%a.y)
            s = self._.div(a.y,a.x)
            print("s=%x"%s)
            l = self._.add(a.x,s)
            print("l=%x"%l)
            ox = self._.square(l)
            print("ox=%x"%ox)
            ox = self._.add(ox,l)
            print("ox=%x"%ox)
            ox = self._.add(ox,self.curve.a2)
            print("ox=%x"%ox)
            oy = self._.square(a.x)
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
            ox = self._.square(l)
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
        x2 = self._.square(p.x)
        x3 = self._.mul(p.x,x2)
        a2x2 = self._.mul(self.curve.a2,x2)
        right = self._.mul(self.curve.a4,p.x)
        right = self._.add(right,self.curve.a6)
        right = self._.add(right,x3)
        right = self._.add(right,a2x2)
        y2 = self._.square(p.y)
        a1xy = self._.mul(self.curve.a1,p.x)
        a1xy = self._.mul(a1xy,p.y)
        a3y = self._.mul(self.curve.a3,p.y)
        left = self._.add(y2,a1xy)
        left = self._.add(left,a3y)
        #print("right=",right)
        #print("left =",left)
        return left==right


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
        print_verilog("z",z)
        r=0
        s=0
        while s==0:
            while r==0:
                k = self.generate_k(e,private_key)
                print("k=%x"%k)
                print_verilog("k",k)
                kG = self._.ScalarMul(k,self.bp)
                print_verilog("kG_x",kG.x)
                r = kG.x % self.bp_order
                print_verilog("r",r)
            s = calc.mul(r,private_key)
            print_verilog("s1",s)
            s = calc.add(s,z)
            print_verilog("s2",s)
            #side channel counter measure, of course a fresh random shall be taken each time rather than a constant
            rnd = 0x12345678
            s = calc.mul(s,rnd)
            k = calc.mul(k,rnd)
            s = calc.div(s,k)
            print_verilog("s",s)
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

    s = calc.plain_mul(0x17232ba853a7e731af129f22ff4149563a419c26bf50a4c9d6eefad6126,0x1db537dece819b7f70f555a67c427a8cd9bf18aeb9b56e0c11056fae6a3)
    print("%x"%s)
    assert(s == 0x188cbf3002f0fceb9905f445cefa0c661ddaf85088e234de26a723857e03aa103f0259eaf60e30a4cf2404f9d08ceea4ad7a9200e8f6eab0f80aa)

    s = calc.plain_mul(
        0xF2C295F27A96B9435935807A7359F67FD014F9A8C9EE2589E13F0CC8B6630CA6,
        0x876E46A6F24CE78C4D904AD897ECC395A14F3DFE78E803FC10D5A8DF4C632923)
    print("%x"%s)
    assert(s == 0x7b9a69879baeac0dc1c9bb962d704bfaa0b5919ab00d1a2227b397211e40747634978dc6bc0cbde09bc5b93a7a8594f15272fcaf1aad7238b8fbcc5d8c93d72a)

    s = calc.plain_mul(
        0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF)
    print("%x"%s)
    assert(s == 0x55555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555)

    modulus=(1<<233)|(1<<74)|1 #q in RFC6979_testvector
    calc = Mod2Arithmetic(modulus)
    s = calc.mul(0x17232ba853a7e731af129f22ff4149563a419c26bf50a4c9d6eefad6126,0x1db537dece819b7f70f555a67c427a8cd9bf18aeb9b56e0c11056fae6a3)
    print("%x"%s)
    assert(s == 0x404c43af73958b87742ff9e35ec83a50fb77c1d266fa5b7e749ddd12ca)



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

    curve=Curve(1,0,0,0,1)
    modulus=(1<<233)|(1<<74)|1 #q in RFC6979_testvector
    base_point=Point(
        0x17232ba853a7e731af129f22ff4149563a419c26bf50a4c9d6eefad6126,
        0x1db537dece819b7f70f555a67c427a8cd9bf18aeb9b56e0c11056fae6a3)
    base_point_order=3450873173395281893717377931138512760570940988862252126328087024741343
    calc = GroupOperation(curve,modulus,Mod2Arithmetic)
    r=calc.PointAdd(base_point,base_point)
    print("r.x=0x%x"%r.x)
    print("r.y=0x%x"%r.y)
    assert(r==Point(
        0x17232ba853a7e731af129f22ff4149563a419c26bf50a4c9d6eefad6126,
        0x1db537dece819b7f70f555a67c427a8cd9bf18aeb9b56e0c11056fae6a3))

def RFC6979_testvector(CalcClass):
    message = bytes("sample", 'utf-8')

    #Curve B163 (FIPS 186-4)
    test_curve(
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
    test_curve(
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
    test_curve(
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
    test_curve(
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
    test_curve(
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
    test_curve(
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
    test_curve(
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
    test_curve(
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

def print_verilog(name,val):
    print("wire [%d-1:0] %s = %d'h%x;"%(val.bit_length(), name, val.bit_length(), val))

def test_curve(calc_class,curve,modulus,arithmetic,base_point,base_point_order,private_key,public_key,message,sig):
    print("qlen=",modulus.bit_length())
    print_verilog("modulus",modulus)
    print_verilog("base_point_x",base_point.x)
    print_verilog("base_point_y",base_point.y)
    print_verilog("base_point_order",base_point_order)
    print_verilog("private_key",private_key)
    print_verilog("public_key_x",public_key.x)
    print_verilog("public_key_y",public_key.y)
    signer = ECDSA_SignatureGenerator(calc_class(curve,modulus,arithmetic),base_point,base_point_order)
    gen_public_key = signer.ComputePublicKey(private_key)
    print("\tpublic_key=",gen_public_key)
    assert(public_key==gen_public_key)
    gen_sig = signer.sign(message,private_key)
    print("\tsignature=",gen_sig)
    assert(sig==gen_sig)
    verifier = ECDSA_SignatureVerifier(calc_class(curve,modulus,arithmetic),base_point,base_point_order)
    verifier.verify(message,gen_sig,gen_public_key)

if __name__ == "__main__":
    testMod2ModArithmetic()
    TestGroupOperationMod2()
    #exit()
    calc = GroupOperation

    #curve: NIST K-233
    test_curve(
        calc_class=calc,
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
        message=bytes("sample", 'utf-8'),
        sig=(0x38AD9C1D2CB29906E7D63C24601AC55736B438FB14F4093D6C32F63A10,
             0x647AAD2599C21B6EE89BE7FF957D98F684B7921DE1FD3CC82C079624F4))

    RFC6979_testvector(calc)
    print("Mod2Arithmetic.plain_mul_s_max=0x%x"%Mod2Arithmetic.plain_mul_s_max)
    print("Mod2Arithmetic.plain_mul_cnt=%d"%Mod2Arithmetic.plain_mul_cnt)
    print("Mod2Arithmetic.mod_mul_s_max=0x%x"%Mod2Arithmetic.mod_mul_s_max)
    print("Mod2Arithmetic.mod_mul_cnt=%d"%Mod2Arithmetic.mod_mul_cnt)
