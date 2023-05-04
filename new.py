"""
Autor Marek Gorski, nr indeksu 155647, grupa D1, semestr II
Wydzial Informatyki, Zarzadzania i Transportu, rok akademicki 2022/2023
"""

import string
import numpy as np
from matplotlib import pyplot as diagram
from numpy.polynomial import polynomial as P

#dane wejściowe
tablicaX = ([-1, 0, 1, 2])
tablicaY = ([5, 6, 4, 7])

print ('Dane wejściowe:')
for i in range(len(tablicaX)):
	print (string.ascii_uppercase[i],': (',tablicaX[i],',',tablicaY[i],')')

#współczynniki wielomianu
wspolczynnikiWielomianu = P.polyfit(tablicaX,tablicaY, len(tablicaX)-1)

aproksymacja1 = P.polyfit(tablicaX,tablicaY, 1)
aproksymacja2 = P.polyfit(tablicaX,tablicaY, 2)

#wyświetlanie współczynników
print ('Wielomian interpolujący:')
print (np.polynomial.Polynomial(wspolczynnikiWielomianu.round(decimals=3)))

print ('Wielomian aproksymujący pierwszego stopnia:')
print (np.polynomial.Polynomial(aproksymacja1.round(decimals=3)))

print ('Wielomian aproksymujący drugiego stopnia:')
print (np.polynomial.Polynomial(aproksymacja2.round(decimals=3)))

x=np.arange (tablicaX[0], tablicaX[-1]+0.01, 0.01)

#konfiguracja diagramu
diagram.grid(linestyle='--') #styl siatki w tle
diagram.title ('Wielomian interpolacji Newtona') #tytuł
diagram.gca().set_aspect('equal') #zachowanie proporcji na osiach XY
diagram.plot(tablicaX, tablicaY,'r.',label='punkty wejściowe') #punkty, kolor czerwony, użycie kropki
diagram.plot(x, P.polyval(x, wspolczynnikiWielomianu), 'g', label='wielomian interpolacji') #dane z zakresu x i wartości do funkcji za pomocą polyval
diagram.plot(x, P.polyval(x, aproksymacja1), 'b', label='wielomian aproksymujący')
diagram.legend(loc='upper left') #legenda w lewym górnym rogu
diagram.show()
