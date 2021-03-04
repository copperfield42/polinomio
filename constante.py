"""
Modulo para constantes arbitraria, para uso en calculos simbolicos
"""
from constanteclass import ConstanteABC, ConstanteBase, sqrt, evaluar_constante, evaluar_formula

import typing
#from abc import abstractmethod
from fractions import Fraction
from numbers import Number, Real, Complex, Integral
#from itertools import chain

from powclass import PowClass

__all__ = ['Constante', 'ConstanteABC', 'ConstanteBase', 'sqrt', 'evaluar_constante', 'evaluar_formula']

class Constante(PowClass, ConstanteBase):
    '''Constante numerica arbitraria de forma: (a+mC**e) con a,m,e conocidos'''

    def __add__(self,otro):
        """C+X"""
        #print("add",self,otro)
        if isinstance(otro,ConstanteABC):
            X,Y = sorted([self,otro],key=lambda x:x.name)
            a1,m1,e1,C1,name1 = X.partes()
            a2,m2,e2,C2,name2 = Y.partes()
            C = C1
            e = e1
            name = name1
            if (C1 and C2 and C1==C2 and e1==e2) or (not C1 and not C2 and name1==name2 and e1==e2):
                m = m1 + m2
                a = a1 + a2
            else:
                m = m1
                a = a1 + Y
            if not m:
                return a
            return self.__class__(name, a=a, m=m, e=e, C=C)
        else:
            return self.__class__(self.name, e=self.e, m=self.m, C=self.C, a=self.a+otro)
        raise NotImplementedError

    def __mul__(self,otro):
        """C * X"""
        #print("mul",self,otro)
        if otro:
            if isinstance(otro,ConstanteABC):
                X,Y = sorted([self,otro],key=lambda x:x.name)
                #print("mul",self,otro,"->",X,Y)
                a1,m1,e1,C1,name1 = X.partes()
                a2,m2,e2,C2,name2 = Y.partes()
                cons = self.__class__
                XX = cons(name1,C=C1,e=e1,m=m1)*a2
                YY = cons(name2,C=C2,e=e2,m=m2)*a1
                AA = a1*a2
                M  = m1*m2
                if (C1 and C2 and C1==C2) or (not C1 and not C2 and name1==name2):
                    e = e1+e2
                    if e:
                        return cons(name1,C=C1,e=e,m=M) + XX + YY + AA
                    else:
                        return XX + YY + M + AA
                XYA = XX + YY + AA
                if e1 == e2 and C1 and C2:
                    e = e1
                    return ((C1*C2)**e)*M + XYA
                return cons(name1,C=C1,e=e1,m=m1*cons(name2,C=C2,e=e2,m=m2)) + XYA
            else:
                return self.__class__(self.name, C=self.C, e=self.e, m=self.m*otro, a=self.a*otro)
        else:
            return 0

    def __pow__(self,otro,modulo=None):
        """C**X
           pow(C,X,m)"""
        if isinstance(otro, ConstanteABC) or not isinstance(otro,Number):
            return NotImplemented
        if modulo is not None:
            if not (isinstance(otro,Integral) and isinstance(modulo,Integral) ):
                raise TypeError("pow() 3rd argument not allowed unless all arguments are integers")
        if not otro:
            return (1%modulo) if modulo is not None else (1)
        elif otro == 1:
            return (self%modulo) if modulo is not None else (+self)
        else:
            # (a+m(C)**e) ** x
            if int(otro) == otro:
                otro = int(otro)
            if otro == 0.5:
                otro = Fraction(1,2)
            a,m,e,C,name = self.partes()
            if a:
                # (a+m(C)**e) ** x
                if otro>0:
                    new = super().__pow__(otro,modulo)
                    if new is not NotImplemented:
                        return new
                else:
                    if modulo is not None:
                        raise ValueError("pow() 2nd argument cannot be negative when 3rd argument specified")
                return self.__class__(self,e=otro)
            else:
                # (m(C)**e) ** x
                if m != 1:
                    m = sqrt(m) if otro == 0.5 else m**otro                            
                if C:
                    if modulo is None:
                        return (C**(e*otro))*m
                    else:
                        return pow(C,e*otro,modulo)*(m%modulo)
                else:
                    return self.__class__(name,e=e*otro,m=m)
                raise NotImplementedError("(m(C)**e) ** x")
        raise NotImplementedError

    def _division(self,otro,div):
        """C/X y/o C//X """
        if not isinstance(otro,Number):
            return NotImplemented
        if otro:
            if self == otro:
                return 1
            if isinstance(otro,ConstanteABC):
                a1,m1,e1,C1,name1 = self.partes()
                a2,m2,e2,C2,name2 = otro.partes()
                cons = self.__class__
                if a2: # (a1+m1(C1)**e1) / (a2+m2(C2)**e2)
                    return self * (otro)**(-1)
                else:  # (a1+m1(C1)**e1) / (m2(C2)**e2)
                    if C2:
                        C3 = C2**(-e2)
                    else:
                        C3 = cons(name2,e=-e2)
                    A = C3*div(a1,m2)
                    if C1:
                        return cons(C1,e=e1,m=div(m1,m2))*C3 + A
                    else:
                        return cons(name1,e=e1,m=div(m1,m2))*C3 + A
            else:
                return self.__class__(self.name, C=self.C, e=self.e, m=div(self.m,otro), a=div(self.a,otro) )            
        else:
            raise ZeroDivisionError



a,b,c = map(Constante,"abc")



