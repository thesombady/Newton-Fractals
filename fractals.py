import sympy
import numpy as np
import matplotlib.pyplot as plt
import scipy

class fractal2d:
    def __init__(self,function,derivitive=None,epsilon=1e-8,delta=0.001,iterations=100):
        if not callable(function):
            raise KeyError("The function isn't Callable")
        if derivitive!=None:
            if not callable(derivitive):
                self.derivitive=self.partialdiv
                print("The derivitive isn't Callable, using the partialdiv function instead.")
        if derivitive==None:
            self.function=function
            self.derivitive=self.partialdiv
        else:
            self.function=function
            self.derivitive=derivitive
        self.epsilon=epsilon
        self.delta=delta
        self.iterations=iterations
    
    def partialdiv(self,x,y):
        if not isinstance(x,(float,int)) and isinstance(y,(int,float)):
            raise TypeError("the input is not of correct type in partialdiv")
        Jacobian_matrix=(1/self.epsilon)*np.array([self.function(x+self.epsilon,y)-self.function(x,y),self.function(x,y+self.epsilon)-self.function(x,y)])
        return Jacobian_matrix


    def newtonmethod(self,guess):
        if not isinstance(guess,(list,tuple,np.array)):
            raise TypeError("The input in newtonmethod isn't of correct type!")
        guess=np.array(guess)
        for i in range(self.iterations):
            if np.linalg.det(self.derivitive(guess[0],guess[1]))==0:
                return None,self.iterations#Diverge
            xy=guess
            guess=guess-np.matmul(np.linalg.inv(self.derivitive(guess[0],guess[1])),self.function(guess[0],guess[1]))
            if abs(xy[0]-guess[0])<self.epsilon and abs(xy[1]-guess[1])<self.epsilon:
                return guess, i#Converge
        return None,self.iterations #Didn't converge, but not necissary diverged.

    def simple_newtonmethod(self,guess):
        if not isinstance(guess,(tuple,list,np.array)):
            raise TypeError("The input is of the wrong type in simple_newtonmethod!")
        guess=np.array(guess)
        Jacobian=self.derivitive(guess[0],guess[1])#Just the jacobian_matrix at point (guess,guess)
        if np.linalg.det(Jacobian)==0:
            return None, self.iterations#Diverges
        Jacobian_inv=np.linalg.inv(Jacobian)
        for i in range(self.iterations):
            xy=guess
            guess=guess-np.matmul(Jacobian_inv,self.function(guess[0],guess[1]))
            if abs(xy[0]-guess[0])<self.epsilon and abs(xy[1]-guess[1])<self.epsilon:
                return guess, i#Yay converged.
        return None,self.iterations
        

    def find_zeros(self,guess,which_newton=None):
        if not isinstance(guess,(list,tuple,np.array)):
            raise TypeError("The input is of wrong type!")
        if which_newton==False:
            self.newton=self.simple_newtonmethod
        else:
            self.newton=self.newtonmethod
        guess=np.array(guess)
        value,i=self.newton(guess[0],guess[1])
        if value==None:#From either newton or simple newton, means no convergence
            return -1,i #i being the number of iterations
        elif value!=None:
            scipy.zeros
            return

    def plot(self,N,a,b,c,d,which_newton=0):
        if which_newton==0:
            self.newton=self.newtonmethod
        elif which_newton==True:
            self.newton=self.simple_newtonmethod
        else:
            self.newton=self.newtonmethod
        pass

def f(x,y):
    return np.array([x**3-3*x*y**2-1,3*x**2*y-y**3])
print(fractal2d(f).partialdiv(0,0))#Print the partial derivitive of f at point (0,0)
print(fractal2d(f).newtonmethod((20,10)))#Prints either a convering or divering number
