import string
import numpy as np
from matplotlib import pyplot as diagram
from numpy.polynomial import polynomial as P

#dane wejściowe
tablicaX = np.array ([-1, 0, 1, 2])
tablicaY = np.array ([5, 6, 4, 7])

#eliminacja Gaussa dla macierzy kwadratowych
def Gauss (A, b):
	n = len(A)
	for i in range(n-1):
		#szukamy mocnych przekątnych
		wiersz = i
		for j in range(i+1, n):
			if abs(A[j][i]) > abs(A[wiersz][i]):
				wiersz = j
		#zamianka
		if wiersz != i:
			b[i], b[wiersz] = b[wiersz], b[i]
			A[i], A[wiersz] = A[wiersz], A[i]
		#tango down! Eliminujemy dolny trójkącik
		for j in range(i+1, n):
			mnoznik = A[j][i] / A[i][i]
			for k in range(i+1, n):
				A[j][k] = A[j][k] - mnoznik * A[i][k]
			b[j] = b[j] - mnoznik * b[i]
	#rozwiązywańsko w górę
	x = [0] * n
	for i in range(n-1, -1, -1):
		x[i] = b[i]
		for j in range(i+1, n):
			x[i] = x[i] - A[i][j] * x[j]
		x[i] = x[i] / A[i][i]
	return np.flip(x)

def AproksymacjaNStopnia (tablicaX, tablicaY, N):
	#metoda 2-go stopnia ma mieć 3 rzędy i 3 kolumny, itd.
	A=np.empty([N+1,N+1])
	b=np.empty([N+1])
	#wypełnianie macierzy od tyłu (od najmniejszych potęg x do najwyższych)
	for i in range (N+1):
		k=i
		b[-i-1]=np.sum((tablicaX**k)*tablicaY)
		for j in range (N+1):
			A[-i-1][-j-1] = np.sum(tablicaX**k)
			k=k+1
	#A=([np.sum(tablicaX**2),np.sum(tablicaX)],[np.sum(tablicaX),len(tablicaX)])
	#b=([np.sum(tablicaX*tablicaY), np.sum(tablicaY)])
	return (Gauss(A,b))

def InterpolacjaNewtona (tablicaX, tablicaY):
	stopien=len(tablicaX)-1

print ('Dane wejściowe:')
for i in range(len(tablicaX)):
	print (string.ascii_uppercase[i],': (',tablicaX[i],',',tablicaY[i],')')

#współczynniki wielomianu
wspolczynnikiWielomianu = P.polyfit(tablicaX,tablicaY, len(tablicaX)-1)

#aproksymacja1 = P.polyfit(tablicaX,tablicaY, 1) #te biblioteki robią za ciebie wszystko
#aproksymacja2 = P.polyfit(tablicaX,tablicaY, 2) #naprawdę wszystko

#wyświetlanie współczynników
print ('Wielomian interpolujący:')
print (np.polynomial.Polynomial(np.round(wspolczynnikiWielomianu,decimals=3)))

N=3
print ('Wielomian aproksymujący ',N,'-go stopnia:')
print(np.polynomial.Polynomial(np.round(AproksymacjaNStopnia(tablicaX, tablicaY, N),decimals=3)))

N=2
print ('Wielomian aproksymujący ',N,'-go stopnia:')
print(np.polynomial.Polynomial(np.round(AproksymacjaNStopnia(tablicaX, tablicaY, N),decimals=3)))

N=1
print ('Wielomian aproksymujący ',N,'-go stopnia:')
print(np.polynomial.Polynomial(np.round(AproksymacjaNStopnia(tablicaX, tablicaY, N),decimals=3)))

N=2
aproksymacja = AproksymacjaNStopnia(tablicaX, tablicaY, N)
x=np.arange (tablicaX[0], tablicaX[-1]+0.01, 0.01)

#konfiguracja diagramu
diagram.grid(linestyle='--') #styl siatki w tle
diagram.title ('Wielomian interpolacji Newtona') #tytuł
diagram.gca().set_aspect('equal') #zachowanie proporcji na osiach XY
diagram.plot(tablicaX, tablicaY,'r.',label='punkty wejściowe') #punkty, kolor czerwony, użycie kropki
diagram.plot(x, P.polyval(x, wspolczynnikiWielomianu), 'g', label='wielomian interpolacji')
diagram.plot(x, P.polyval(x, aproksymacja), 'b', label='wielomian aproksymujący')
diagram.legend(loc='upper left') #legenda w lewym górnym rogu
diagram.show()
