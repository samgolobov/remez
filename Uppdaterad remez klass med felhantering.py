from numpy import *
from matplotlib.pyplot import *
from numpy.polynomial import Polynomial
from scipy import optimize

class RemezApprox:
    def __init__(self,f):
        if not inspect.isfunction(f):
            raise TypeError(f'Expected argument should be of type function, not {type(f)}')
        self.func=f

    def __createPolynom(self,x):
            n=self.n
            M=zeros((n+2,n+2)) #Skapa en n+2 x n+2 matris fylld med nollor
            for i in range(n+2): #Första kolumnen bara ettor? Fyll matrisen rad för rad, förutom sista kolumnen.
                    for j in range(n+1):
                        M[i][j]=x[i]**j

            for i in range(n+2): #Sista kolumnen, fylld med 1 och -1.
                M[i][-1]=(-1)**i
            F=self.func(x) #Kända funktionsvärden för testfunktionen i värdena x0,x1,x2 ... i vårt intervall.
            self.Pcoeff= linalg.solve(M,F) #Löser linjära ekvationsystemet. Får ut koefficienterna för polynomet, där sista elementet är h.  
            p=Polynomial(self.Pcoeff[:-1]) #Polynomial skapar ett polynom objekt med värdena i Pcoeff förutom sista värdet 'h'.
            return p
    
    def remez(self,n,a,b, maxIteration = 50, tolerance = 1e-6):
        if not isinstance(n, int):
            raise TypeError(f'The polynomial degree needs to be of type Int, not {type(n)}')
        if not isinstance(a, (int, float)):
            raise TypeError(f'Reference value needs to be of type Int/Float, not {type(a)}')
        if not isinstance(b, (int, float)):
            raise TypeError(f'Reference value needs to be of type Int/Float, not {type(b)}')
        if not isinstance(maxIteration, int):
            raise TypeError(f'Number of max iterations needs to be of type Int, not {type(maxIteration)}')
        if not isinstance(tolerance, (int, float)):
            raise TypeError(f'The tolerance needs to be of type Int/Float, not {type(tolerance)}')

        self.tol = tolerance
        self.maxiter = maxIteration
        self.n=n

        x=array(linspace(a,b,n+2))
        self.x=x

        t=linspace(a,b,int(1e4))
        self.t=t
        
        f = self.func
        p = self.__createPolynom(x)

        counter=0
        Pcoeff=self.Pcoeff

        while abs(Pcoeff[-1]-max(abs(p(t)-f(t))))>self.tol and counter<self.maxiter:
            counter+=1
            self.p = self.__createPolynom(x) #Skapar nytt polynom med efter bytt ut en punkt i referensen.
            L = [(abs(f(i)-p(i))) for i in t] #Används för att hitta ett max i ett tätare intervall istället för derivera. 
            max_value = max(L)
            my = t[L.index(max_value)]

        #Steg 4.
            if my<x[0]:
                if sign(f(my)-p(my)) == sign(f(x[0])-p(x[0])):
                    x[0]=my
                else: 
                    x[-1]=my

            if my>x[-1]:
                if sign(f(my)-p(my)) == sign(f(x[-1])-p(x[-1])):
                    x[-1]=my
                else:
                    x[0]=my

            else:
                for i in range(len(x)):
                    if (my<=x[i]) :
                        x_i=i-1
                        break
                if sign(f(my)-p(my)) == sign(f(x[x_i])-p(x[x_i])):
                    x[x_i] = my
                else:
                    x[x_i+1] = my 

        # self.polynomial = p(x) # Kommer inte till använding i koden?
        # self.counter=counter # Kommer inte till använding i koden?

    def plotRemez(self):
        t=self.t
        figure(0)
        plot(t,self.p(t))
        plot(t,self.func(t))

    def plotError(self):
        x=self.x
        t=self.t
        figure(1)
        plot(t,(abs(self.p(t)-self.func(t)))) #Plottar error funktionen.
        yscale('log')
        plot(x,1+0*x,'*')
        
def f(x): #definera en testfunktion
    return sin(x)
G=RemezApprox(f)
G.remez(10,0,10) #degree,a,b
G.plotRemez()
G.plotError()

