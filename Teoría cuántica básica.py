'TEORÍA CUÁNTICA BÁSICA'
'Librería comoutación cuántica'
#Nota:Es necesario importar la calculadora de números complejos y numpy.
import numpy as np
from numpy import linalg as LA
import ComplejosCal
import math

'--------------------- Probabilidad de encontrar una partícula en una posición------------------------'

#Probabilidad de la particula.(Permite calcular la probabilidad de encontrar una particula en una posición dada)
#La probabilidad se realiza  dividiendo la posición al cuadrado entre el módulo al cuadrado del vector 
def probabilidadPosicion(vk,pos):
    norma = ComplejosCal.normaVector(vk)
    probabilidad = ComplejosCal.modulo(vk[pos])**2/(norma)**2
    print('La probabilidad de la particula de estar en la posición',pos,'es:')
    print(probabilidad)

#Probabilidad de la particula.(Permite calcular la probabilidad de encontrar una particula todas las posiciones)
#La probabilidad se realiza  dividiendo la posición al cuadrado entre el módulo al cuadrado del vector 
def probabilidadPosiciones(vk):
    norma = ComplejosCal.normaVector(vk)
    for i in range(len(vk)):
        vk[i] = ComplejosCal.modulo(vk[i])**2/(norma)**2
        
    print('El vector de probabilidad es:')
    print(vk)
    print('')
    
    for j in range(len(vk)):
        grafico= '*' * round(vk[j]*100)
        print(j,grafico)

'---------------------------------Amplitud de transición------------------------------------------------'

#Amplitud de transición.(Permite hallar la probabilidad de transitar de un vector a otro después de hacer la observación).
#La amplitud se halla con el producto interno de los dos vectores normalizados
def amplitudTransicion(v1,v2):
    normalizado1 = normalizarVector(v1)
    normalizado2 = normalizarVector(v2)
    interno = ComplejosCal.productoInterno(normalizado1,normalizado2)
    print('La amplitud de transición de los vectores es:')
    print(interno)

#Nomalizar vector.(Permite normalizar un vector)
#La normalización se halla dividiendo cada número complejo del vector entre la norma del vector
def normalizarVector(vk):
    vkNormalizado = []
    norma = ComplejosCal.normaVector(vk)
    for i in range(len(vk)):
        vkNormalizado.append([])
        vkNormalizado[i] = divisionKet(vk[i],norma)
    return vkNormalizado

#División para normalizar.(Permite hallar un número complejo dividido por la norma del vector)
#La división se halla con el número complejo y la norma del vector
def divisionKet(num,norma):
    lista = []
    real = num[0]/norma
    imaginario = num[1]/norma
    lista.append(real)
    lista.append(imaginario)
    return lista

'------------------------------------Valor esperado y varianza-------------------------------------------'

#Valor esperado y varianza.(Permite hallar el valor esperado y la varianza de un observable con respecto a un estado)
def varianza_valor(m1,vk):
    verificar = ComplejosCal.matrizHermitiana(m1)
    if (verificar == True):
        print('El valor esperado del observable es:')
        print(valorEsperado(m1,vk))
        print('')
        print('La varianza del observable es:')
        print(varianza(m1,vk))
        
    else:
        print('El observable no es una matriz hermitiana')
        
#Valor esperado.(Permite hallar el valor esperado de un observable con respecto a un estado)
#El valor esperado se halla con el producto interno de la multiplicación (la matriz y el vector) y el vector
def valorEsperado(m1,vk):
    vNor = normalizarVector(vk)
    accion = ComplejosCal.accionMatriz(m1,vNor)
    return (ComplejosCal.productoInterno(accion,vNor))

#Varianza.(Permite hallar la varianza de un observable con respecto a un estado)
#La varianza se halla encontrando el valor esperado de la resta del observable y (el valor esperado por la matriz identidad) al cuadrado
def varianza(m1,vk):
    valor = valorEsperado(m1,vk)
    vNor = normalizarVector(vk)
    resultado = ComplejosCal.escalarMatriz(valor,ComplejosCal.matrizIdentidad(len(m1),len(m1[0])))
    resta = ComplejosCal.restaMatriz(m1,resultado)
    restaCuadrada= ComplejosCal.productoMatriz(resta,resta)
    return(valorEsperado(restaCuadrada,vNor))

'------------------------------------Valores y vectores propios-------------------------------------------'

#Valores y vectores propios.(Permite calcular los valores y vectores propios de un observable)
#Los valores y vectores propios se hallan por medio de la calculadora numpy 
def vPropios(m1):
    A = np.array(m1)
    valores, vectores = LA.eig(A)
    print('Vectores propios')
    for i in vectores:
      print(i)
    print('')
    print('Valores propios')
    print(valores)

