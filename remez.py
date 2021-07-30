from numpy import *
from matplotlib.pyplot import *
from numpy.polynomial import Polynomial
from scipy import optimize



n=2 #graden på polynomet
a=0 #startvärde
b=1 #slutvärde

x=array(linspace(a,b,n+2)) # vårt intervall mellan a och b, n+2 punkter
x1=array(linspace(0,1,100))# xaxel med tätare xvärden för att testa funktionen


def testf(x): #definera en testfunktion
    return exp(x)
X=zeros((n+2,n+2)) #Skapa en n+2 x n+2 matris fylld med nollor
for i in range(n+2): #Första kolumnen bara ettor? Fyll matrisen rad för rad, förutom sista kolumnen.
    for j in range(n+1):
        X[i][j]=x[i]**j

for i in range(n+2): #Sista kolumnen, fylld med 1 och -1.
    X[i][-1]=(-1)**i
F=testf(x) #Kända funktionsvärden för testfunktionen i värdena x0,x1,x2 ... i vårt intervall


Pcoeff= linalg.inv(X)@(F) #Löser linjära ekvationsystemet. Får ut koefficienterna för polynomet, där sista elementet är h. 



p=Polynomial(Pcoeff[:-1]) #Polynomial skapar ett polynom objekt med värdena i Pcoeff förutom sista värdet 'h'.
figure(0)
plot(x1,p(x1))
plot(x1,testf(x1))
figure(1)
plot(x1,abs(p(x1)-testf(x1)))

#I steg 3 vill vi hitta min eller max till error funktionen abs(f(x)-p(x)). Derivera den och hitta den punkt my där skillnaden är störst.

f2=lambda x:abs((p.deriv())(x)-exp(x))
plot(x1,f2(x1))
figure(2)
z=optimize.newton(f2,1)
print(f'root = {z}')
print(f'f2(z) = {f2(z)}')
#Alternativt istället enbart hitta maximum, eftersom det står i uppgiften att det inte behövs optimeras.
my=max((abs(p(x1)-testf(x1))))


print(f'h={abs(Pcoeff[-1])}')
print(p)
print(f'my={my}')
max((abs(p(x)-testf(x))))

#Steg 1
#List comprehension. Mina punkter x0,x1,x2,x3 i testintervallet (a,b)
#ekvation ett 
#(sista punkten som h/epsilon?) n+2 värden.
#variabler p0,p1,p2,p3 ... sökta koefficienter för mitt polynom i en kolumnmatris. givna värden från funktionen som skall approximeras f(x0),f(x1),f(x2)
# i en columnn. Lös linjära systemet för p0,p1,p2.
#Steg 2

class RemezApprox:
    def __init__(self,f):
        self.func=f
    def remez(self,n):
        a=0
        b=1
        self.degree=n
        x=array(linspace(a,b,n+2))
        
        X=zeros((n+2,n+2))
        for i in range(n+2): #Första kolumnen bara ettor?
            for j in range(n+1):
                X[i][j]=x[i]**j
        for i in range(n+2):
            X[i][-1]=(-1)**i
        F=self.func(x)
        Pcoeff= linalg.inv(X)@(F)
        self.p=Polynomial(Pcoeff[:-1])
    def plot(self,xval):
        figure(0)
        plot(xval,self.p(xval))
        plot(xval,self.func(xval))

fz=lambda x:exp(x)
E = RemezApprox(fz)

E.remez(4)

E.plot(x1)

test
