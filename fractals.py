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
        self.zeros=np.array([[]])
    
    def partialdiv(self,x,y):
        if not isinstance(x,(float,int)) and isinstance(y,(int,float)):
            raise TypeError("the input is not of correct type in partialdiv")
        Jacobian_matrix=(1/self.epsilon)*np.array([self.function(x+self.epsilon,y)-self.function(x,y),self.function(x,y+self.epsilon)-self.function(x,y)])
        return np.transpose(Jacobian_matrix)


    def newtonmethod(self,guess):
        if not isinstance(guess,(list,tuple,np.ndarray,np.generic)):
            raise TypeError("The input in newtonmethod isn't of correct type!")
        guess=np.array(guess)
        for i in range(self.iterations):
            if np.linalg.det(self.derivitive(guess[0],guess[1]))==0:
                return None,self.iterations#Diverge
            xy=guess
            guess=guess-np.matmul(np.linalg.inv(self.derivitive(guess[0],guess[1])),self.function(guess[0],guess[1]))
            if abs(xy[0]-guess[0])<self.epsilon and abs(xy[1]-guess[1])<self.epsilon:
                return guess, i#Converge
        return None,self.iterations #Didn't converge, but not necissary diverged. Might have needed more repititions/iterations

    def simple_newtonmethod(self,guess):
        if not isinstance(guess,(tuple,list,np.ndarray,np.generic)):
            raise TypeError("The input is of the wrong type in simple_newtonmethod!")
        guess=np.array(guess)
        Jacobian=self.derivitive(guess[0],guess[1])#Just the jacobian_matrix at point (guess)
        if np.linalg.det(Jacobian)==0:
            return None, self.iterations#Diverges
        Jacobian_inv=np.linalg.inv(Jacobian)
        for i in range(self.iterations):
            xy=guess
            guess=guess-np.matmul(Jacobian_inv,self.function(guess[0],guess[1]))
            if abs(guess).all()>(1/self.epsilon):
                return None,self.iterations
            if abs(xy[0]-guess[0])<self.epsilon and abs(xy[1]-guess[1])<self.epsilon:
                return guess, i#Yay converged.
        return None,self.iterations#Didn't converge but not neccissary diverged.
        

    def find_zeros(self,guess,which_newton=None):
        if not isinstance(guess,(list,tuple,np.ndarray,np.generic)):
            raise TypeError("The input is of wrong type!")
        if which_newton==False:
            self.newton=self.simple_newtonmethod
        else:
            self.newton=self.newtonmethod
        guess=np.array(guess)
        value,i=self.newton(guess)
        if value is None:#From either newton or simple newton, means no convergence
            return -1,i #i being the number of iterations
        
        if self.zeros.size == 0: # If not root existed before this, add this value, "value" as a root
            self.zeros = np.array([value])

        q = (self.zeros-value)**2
        distance = q[:,0] + q[:,1]
        if (distance<self.delta).any() : #zero already exists
            return np.where(distance<self.delta)[0][0] , i
        else:  # value dose not exist yet
            self.zeros=np.reshape(np.append(self.zeros,value),(-1,2)) #add value to zeros
            return self.zeros.size/2-1 , i

    def callfind_zeros(self,x,y,which_newton=None):
        guess=(x,y)
        return self.find_zeros(guess,which_newton)

    def plot(self,N,a,b,c,d,which_newton=None):
        X,Y=np.meshgrid(np.linspace(a,b,N),np.linspace(c,d,N),indexing='ij')
        v_zeroes=np.vectorize(self.callfind_zeros)
        A,B=v_zeroes(X,Y,which_newton)#B is the number of iterations
        plt.pcolor(B,cmap='viridis')
        plt.show()
        return A,B
    
    def plot2(self,N,a,b,c,d,which_newton=None):
        X,Y=np.meshgrid(np.linspace(a,b,N),np.linspace(c,d,N),indexing='ij')
        v_zeroes=np.vectorize(self.callfind_zeros)
        A,B=v_zeroes(X,Y,which_newton)#B is the number of iterations
        #B = ((B-1)%10)/10
        plt.pcolormesh(B,cmap='viridis')
        plt.show()
        return A,B


def f(x,y):
    return np.array([x**3-3*x*y**2-1,3*x**2*y-y**3])
#print(fractal2d(f).partialdiv(0,0))#Print the partial derivitive of f at point (0,0)
#print(fractal2d(f).newtonmethod((20,10)))#Prints either a convering or divering number
#print(fractal2d(f).plot(100,-10,10,-10,10))
#fractal2d(f).plot(200,-1,1,-1,1)
def h(x,y):
    return np.array([x**8-28*x**6*y**2+70*x**4*y**4+15*x**4-28*x**2*y**6-90*x**2*y**2+y**8+15*y**4-16,8*x**7*y-56*x**5*y**3+56*x**3*y**5+60*x**3*y-8*x*y**7-60*x*y**3])
#fractal2d(h).plot(200,-1,1,-1,1)
def g(x,y):
    return np.array([x**3-3*x*y**2-2*x-2,3*x**2*y-y**3-2*y])

#fractal2d(h).plot2(1000,-10,10,-10,10)
fractal2d(h).plot2(100,-10,10,-10,10)