from django.shortcuts import render
from .forms import *
import csv
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier
import numpy as np

# Create your views here.

def classificar_proteina(request):
        rato = ''
        veredicto = ''
        teste = []
        if request.method == 'POST':
                form = ContactCourse(request.POST)
                if form.is_valid():
                        rato = form.cleaned_data['sequencia']
                        rato = rato.replace('\n', '')
                        rato = rato.replace('\r', '')

                        # DICIONÁRIO DE VALORES DE HIDROPATIA DOS 20 AMINOÁCIDOS
                        # Fonte: Kyte e Doolittle (1982)
                        # ----------------------------------------------------------
                        kdha = {}
                        kdha['I'] = 4.5
                        kdha['V'] = 4.2
                        kdha['L'] = 3.8
                        kdha['F'] = 2.8
                        kdha['C'] = 2.5
                        kdha['M'] = 1.9
                        kdha['A'] = 1.8
                        kdha['G'] = -0.4
                        kdha['T'] = -0.7
                        kdha['S'] = -0.8
                        kdha['W'] = -0.9
                        kdha['Y'] = -1.3
                        kdha['P'] = -1.6
                        kdha['H'] = -3.2
                        kdha['E'] = -3.5
                        kdha['Q'] = -3.5
                        kdha['D'] = -3.5
                        kdha['N'] = -3.5
                        kdha['K'] = -3.9
                        kdha['R'] = -4.5

                        # ----------------------------------------------------------

                        # HIDROPATIA TOTAL => SOMA DAS HIDROPATIAS DE TODOS OS AMINOÁCIDOS
                        # QUE COMPÕE A PROTEÍNA
                        # ----------------------------------------------------------
                        def hidropatia_total(sequencia_aminoacidos):
                                soma = 0
                                for i in range(len(sequencia_aminoacidos)):
                                        soma += kdha[sequencia_aminoacidos[i]]
                                return soma

                        # ----------------------------------------------------------

                        # HIDROPATIA C TERMINAL => SOMA DAS HIDROPATIAS DE TODOS OS AMINOÁCIDOS
                        # QUE COMPÕE A REGIÃO DO C TERMINAL (25 ULTIMOS)
                        # ----------------------------------------------------------
                        def hidropatia_c_terminal(sequencia_aminoacidos):
                                soma = 0
                                cterminal = sequencia_aminoacidos[-25::]
                                for i in range(len(cterminal)):
                                        soma += kdha[cterminal[i]]
                                return soma

                        # ----------------------------------------------------------

                        # CARGA CALCULADA APENAS COM OS 25 ULTIMOS AMINOÁCIDOS DA SEQUENCIA
                        # AMINOACIDOS H R K SOMA MAIS UM NA CARGA
                        # AMINOACIDOS E D SUBTRAI UM DA CARGA
                        # ----------------------------------------------------------
                        def carga_c_terminal(sequencia_aminoacidos):
                                soma = 0
                                cterminal = sequencia_aminoacidos[-25::]
                                for i in range(len(cterminal)):
                                        if (cterminal[i] == 'H' or cterminal[i] == 'R' or cterminal[i] == 'K'):
                                                soma += 1
                                        if cterminal[i] == 'E' or cterminal[i] == 'D':
                                                soma -= 1
                                return soma

                        # ----------------------------------------------------------

                        # QUANTIDADE CALCULADA APENAS COM OS 25 ULTIMOS AMINOÁCIDOS DA SEQUENCIA
                        # AMINOACIDOS H R K SOMA MAIS UM NA CARGA
                        # AMINOACIDOS POLARES BASICOS
                        # INDICA A QUANTIDADE DE ALCALINOS NO C TERMINAL
                        # ----------------------------------------------------------
                        def aminoacidos_basicos_c_terminal(sequencia_aminoacidos):
                                soma = 0
                                cterminal = sequencia_aminoacidos[-25::]
                                for i in range(len(cterminal)):
                                        if (cterminal[i] == 'H' or cterminal[i] == 'R' or cterminal[i] == 'K'):
                                                soma += 1
                                return soma

                        # ----------------------------------------------------------

                        # SINAL RECONHECIDO POR UMA SEQUENCIA
                        # TEM Q TER PKKKRKV
                        # E PELO MENOS 4 KR
                        # SINAL QUE INDICA O CAMINHO DO NUCLEO PARA A PROTEÍNA
                        # ----------------------------------------------------------
                        def sinal_localizacao_nuclear(sequencia_aminoacidos):
                                x = sequencia_aminoacidos.find('PKKKRKV')
                                if x != -1:
                                        return 1
                                x = sequencia_aminoacidos.count('KR')
                                if x >= 4:
                                        return 1
                                else:
                                        return 0

                                # x = sequencia_aminoacidos.find('PKKKRKV')
                                # y = sequencia_aminoacidos.count('KR')
                                # print(x)
                                # print(y)
                                # if ((x != -1) or (y >= 4)):
                                #   return 1

                        # else:
                        #   return 0
                        # ----------------------------------------------------------

                        # TEM QUE ENCONTRAR NO MINIMO 2 SEQUENCIAS DE 7 AMINOACIDOS
                        # PRIMEIRO E O QUARTO HIDROFOBICOS (POSITIVOS)
                        # SEGUNDO TERCEIRO QUINTO SEXTO SETIMO HIDROFILICOS (NEGATIVOS)
                        # ----------------------------------------------------------
                        def dominio_espiral_enrolada(sequencia_aminoacidos):
                                inicio = 0
                                fim = 7
                                ocorrencia = 0
                                while fim <= len(sequencia_aminoacidos):
                                        soma = 0
                                        for i in range(inicio, fim):
                                                posicao = i - inicio
                                                valor = kdha[sequencia_aminoacidos[i]]
                                                if (posicao == 0 or posicao == 3) and valor > 0:
                                                        soma += 1
                                                if (
                                                        posicao == 1 or posicao == 2 or posicao == 4 or posicao == 5 or posicao == 6) and valor < 0:
                                                        soma += 1
                                        if soma == 7:
                                                ocorrencia += 1
                                        inicio += 1
                                        fim += 1
                                if ocorrencia >= 2:
                                        return 1
                                else:
                                        return 0
                                # soma = 0
                                # for i in range(len(sequencia_aminoacidos)):
                                #   if(sequencia_aminoacidos[i] == 'I' or sequencia_aminoacidos[i] == 'V' or sequencia_aminoacidos[i] == 'L' or sequencia_aminoacidos[i] == 'C' or sequencia_aminoacidos[i] == 'M' or sequencia_aminoacidos[i] == 'A'):
                                # print(sequencia_aminoacidos[i])
                                #      if(sequencia_aminoacidos[i+1] != 'I' and sequencia_aminoacidos[i+1] != 'V' and sequencia_aminoacidos[i+1] != 'L' and sequencia_aminoacidos[i+1] != 'C' and sequencia_aminoacidos[i+1] != 'M' and sequencia_aminoacidos[i+1] != 'A'):
                                # print(sequencia_aminoacidos[i+1])
                                #         if(sequencia_aminoacidos[i+2] != 'I' and sequencia_aminoacidos[i+2] != 'V' and sequencia_aminoacidos[i+2] != 'L' and sequencia_aminoacidos[i+2] != 'C' and sequencia_aminoacidos[i+2] != 'M' and sequencia_aminoacidos[i+2] != 'A'):
                                # print(sequencia_aminoacidos[i+2])
                                #            if(sequencia_aminoacidos[i+3] == 'I' or sequencia_aminoacidos[i+3] == 'V' or sequencia_aminoacidos[i+3] == 'L' or sequencia_aminoacidos[i+3] == 'C' or sequencia_aminoacidos[i+3] == 'M' or sequencia_aminoacidos[i+3] == 'A'):
                                # print(sequencia_aminoacidos[i+3])
                                #               if(sequencia_aminoacidos[i+4] != 'I' and sequencia_aminoacidos[i+4] != 'V' and sequencia_aminoacidos[i+4] != 'L' and sequencia_aminoacidos[i+4] != 'C' and sequencia_aminoacidos[i+4] != 'M' and sequencia_aminoacidos[i+4] != 'A'):
                                # print(sequencia_aminoacidos[i+4])
                                #                  if(sequencia_aminoacidos[i+5] != 'I' and sequencia_aminoacidos[i+5] != 'V' and sequencia_aminoacidos[i+5] != 'L' and sequencia_aminoacidos[i+5] != 'C' and sequencia_aminoacidos[i+5] != 'M' and sequencia_aminoacidos[i+5] != 'A'):
                                # print(sequencia_aminoacidos[i+5])
                                #                     if(sequencia_aminoacidos[i+6] != 'I' and sequencia_aminoacidos[i+6] != 'V' and sequencia_aminoacidos[i+6] != 'L' and sequencia_aminoacidos[i+6] != 'C' and sequencia_aminoacidos[i+6] != 'M' and sequencia_aminoacidos[i+6] != 'A'):
                                # print(sequencia_aminoacidos[i+6])
                                #                        soma = soma + 1
                                # if(soma >= 2):
                                #   return 1
                                # else:
                                #   return 0

                        # ----------------------------------------------------------

                        # TEM QUE ENCONTRAR NO MINIMO 1 SEQUENCIA COMEÇANDO COM C
                        # SEGUIDO DE I OU V OU L OU A OU P OU M
                        # SEGUIDO DE I OU V OU L OU A OU P OU M
                        # SEGUIDO DE A OU M OU S OU L
                        # ----------------------------------------------------------

                        def dominio_prenilacao(sequencia_aminoacidos):

                                inicio = 0
                                fim = 4
                                resposta = 0
                                sequencia_aminoacidos = sequencia_aminoacidos[-25::]
                                while fim <= len(sequencia_aminoacidos):
                                        if sequencia_aminoacidos[inicio] == 'C' and (
                                                sequencia_aminoacidos[inicio + 1] == 'I' or sequencia_aminoacidos[
                                                inicio + 1] == 'V' or sequencia_aminoacidos[inicio + 1] == 'L' or
                                                sequencia_aminoacidos[inicio + 1] == 'A' or sequencia_aminoacidos[
                                                        inicio + 1] == 'P' or sequencia_aminoacidos[
                                                        inicio + 1] == 'M') and (
                                                sequencia_aminoacidos[inicio + 2] == 'I' or sequencia_aminoacidos[
                                                inicio + 2] == 'V' or sequencia_aminoacidos[inicio + 2] == 'L' or
                                                sequencia_aminoacidos[inicio + 2] == 'A' or sequencia_aminoacidos[
                                                        inicio + 2] == 'P' or sequencia_aminoacidos[
                                                        inicio + 2] == 'M') and (
                                                sequencia_aminoacidos[inicio + 3] == 'A' or sequencia_aminoacidos[
                                                inicio + 3] == 'M' or sequencia_aminoacidos[inicio + 3] == 'S' or
                                                sequencia_aminoacidos[inicio + 3] == 'L'):
                                                resposta = 1
                                        inicio += 1
                                        fim += 1
                                return resposta
                                # soma = 0
                                # for i in range(len(sequencia_aminoacidos)):
                                #    if sequencia_aminoacidos[i] == 'C' and (sequencia_aminoacidos[i+1] == 'I' or sequencia_aminoacidos[i+1] == 'V' or sequencia_aminoacidos[i+1] == 'L' or sequencia_aminoacidos[i+1] == 'A' or sequencia_aminoacidos[i+1] == 'P' or sequencia_aminoacidos[i+1] == 'M') and (sequencia_aminoacidos[i+2] == 'I' or sequencia_aminoacidos[i+2] == 'V' or sequencia_aminoacidos[i+2] == 'L' or sequencia_aminoacidos[i+2] == 'A' or sequencia_aminoacidos[i+2] == 'P' or sequencia_aminoacidos[i+2] == 'M') and (sequencia_aminoacidos[i+3] == 'A' or sequencia_aminoacidos[i+3] == 'M' or sequencia_aminoacidos[i+3] == 'S' or sequencia_aminoacidos[i+3] == 'L'):
                                #        soma = 1
                                # if(soma == 1):
                                #    return 1
                                # else:
                                #    return 0

                        # ----------------------------------------------------------

                        # TEM QUE SER CAPAZ DE GERAR DUPLA HELICE
                        # TEM Q TER PREDOMINANCIA DE R L S A
                        # TEM Q TER ESCASSEZ DE E D
                        # ----------------------------------------------------------

                        def sinal_localizacao_mitocondrial(sequencia_aminoacidos):
                                # espiral = espiraldupla(sequencia_aminoacidos)
                                soma = 0
                                acidos = 0
                                nterminal = sequencia_aminoacidos[0:20]
                                for i in range(len(nterminal)):
                                        if sequencia_aminoacidos[i] == 'R' or sequencia_aminoacidos[i] == 'L' or \
                                                sequencia_aminoacidos[i] == 'S' or sequencia_aminoacidos[i] == 'A':
                                                soma += 1
                                        if sequencia_aminoacidos[i] == 'H' or sequencia_aminoacidos[i] == 'R' or \
                                                sequencia_aminoacidos[i] == 'K':
                                                acidos += 1
                                if acidos >= 2 and soma > 5:
                                        return 1
                                else:
                                        return 0

                                # acidos = 0
                                # soma = 0
                                # espiral = dominio_espiral_enrolada(sequencia_aminoacidos)
                                # for i in range(25):
                                #    if sequencia_aminoacidos[i] == 'R' or sequencia_aminoacidos[i] == 'L' or sequencia_aminoacidos[i] == 'S' or sequencia_aminoacidos[i] == 'A':
                                #        soma += 1
                                #    if sequencia_aminoacidos[i] == 'H' or sequencia_aminoacidos[i] == 'R' or sequencia_aminoacidos[i] == 'K':
                                #        acidos += 1
                                # if acidos >= 2 and soma > 5 and espiral == 1:

                        # if acidos >= 2 and soma > 5:
                        #    return 1
                        # else:
                        #    return 0

                        # ----------------------------------------------------------
                        # ABRE O ARQUIVO QUE CONTÉM AS SEQUENCIAS DAS PROTEÍNAS QUE VÃO
                        # COMPOR A BASE DE TREINAMENTO (SOMENTE PARA LEITURA)
                        #arq_base_treinamento = open('base_treinamento_classificador.txt', 'r')
                        arq_base_treinamento = rato

                        # ABRE O ARQUIVO QUE CONTÉM AS SEQUENCIAS DAS PROTEÍNAS QUE VÃO
                        # COMPOR A BASE DE CASOS DE TESTE (SOMENTE PARA LEITURA)
                        # arq_base_caso_teste = open('base_casos_teste.txt','r')

                        # ABRE O ARQUIVO QUE VAI RECEBER OS ATRIBUTOS QUE SERÃO RETIRADOS
                        # DA SEQUÊNCIA DE AMINOÁCIDOS (APAGA O CONTEÚDO ARQUIVO) (TREINAMENTO)
                        base_treinamento = open('base_treinamento_classificador.csv', 'w', newline='')

                        # ABRE O ARQUIVO QUE VAI RECEBER OS ATRIBUTOS QUE SERÃO RETIRADOS
                        # DA SEQUÊNCIA DE AMINOÁCIDOS (APAGA O CONTEÚDO ARQUIVO) (CASOS DE TESTE)
                        # base_caso_teste = open('base_casos_teste.csv','w',newline='')
                        # ----------------------------------------------------------

                        arqcompc = ''

                        for line in arq_base_treinamento:
                                arqcompc += line

                        sequencias = []

                        inicio = False

                        fim = False

                        arqcompc.replace('\n', '')
                        arqcompc.replace('\r', '')

                        # print arqcompc
                        sequencia = ''

                        for i in range(len(arqcompc)):
                                if arqcompc[i] == '>':
                                        inicio = not inicio
                                if arqcompc[i] == ']':
                                        fim = not fim
                                        if sequencia != '':
                                                sequencias.append(sequencia)
                                        sequencia = ''
                                if (inicio and fim) or (not inicio and not fim):
                                        if arqcompc[i] != '\n' and arqcompc[i] != ']':
                                                sequencia += arqcompc[i]
                        sequencias.append(sequencia)

                        # LENDO O ARQUIVO INTEIRO E ATRIBUINDO A VARIAVEL
                        # unica_string = arq_base_treinamento.read()

                        # ESSA VARIAVEL VAI RECEBER A SEQUENCIA DE AMINOÁCIDOS QUE COMPÕE A PROTEÍNA
                        # sequencia = ""
                        # VARIAVEL CONTADOR
                        # k = 0
                        # LISTA QUE VAI RECEBER TODAS AS SEQUÊNCIAS DE AMINOACIDOS ENCONTRADAS
                        # sequencias = []

                        # ESSE FOR VAI PASSAR EM TODOS OS CARACTERES DO ARQUIVO
                        # for i in range(len(unica_string)):
                        # SE O CARACTERE FOR UM @ INDICA O INICIO DE UMA SEQUENCIA DE AMINOÁCIDOS
                        # if (unica_string[i] == "@"):
                        # J & K VARIÁVEIS CONTADOR
                        #    j = i + 1
                        #   k = k + 1
                        # ENQUANTO A CADEIA DE AMINOÁCIDOS NAO ACABAR A VARIAVEL SEQUENCIA VAI CONCATENAR
                        # OS AMINOÁCIDOS QUE LER
                        # QUANDO CHEGAR NO # INDICA QUE O PADRÃO ACABOU
                        #   while (unica_string[j] != "#"):
                        #         sequencia += unica_string[j]
                        # CONTADOR PARA PERCORRER TODA A SEQUENCIA
                        #          j = j + 1
                        # INSERINDO A SEQUENCIA QUE FOI LIDA NA LISTA DE SEQUENCIAS
                        # sequencias.append(sequencia)
                        # ZERANDO A VARIAVEL SEQUENCIA PARA LER A PROXIMA CADEIA
                        # sequencia = ""

                        # CRIA UMA LISTA DE CARACTERISTICAS DAS PROTEINAS EFETORAS
                        lista_caracteristicas = []
                        # CRIA UMA INSTANCIA PARA ESCREVER NA BASE DE TREINAMENTO
                        writer = csv.writer(base_treinamento)
                        cabecalho = ['Hidropatia Total', 'Hidropatia Media', 'Hidropatia C-Terminal',
                                     'Carga C-Terminal', 'Aminoacidos Basicos C-Terminal', 'Dupla Helice',
                                     'Sinal Nuclear', 'Sinal Mitocondrial', 'Dominio Prenilacao']
                        writer.writerow(cabecalho)

                        caracteristicas_arredondadas = []
                        #caracteristicas_algoritmo = []

                        # UM FOR PARA PASSAR EM TODAS AS SEQUENCIAS DE AMINOACIDOS
                        for z in range(len(sequencias)):
                                # hidrototal = hidromedia = hidrocterminal = cargacterminal =
                                # aminocterminal = slm = sln = de = dp
                                # AQUI SÃO INSERIDOS OS 9 VALORES QUE SÃO EXTRAÍDOS DA SEQUENCIA
                                # 9 CARACTERISTICAS DAS PROTEINAS EFETORAS

                                # CHAMANDO A FUNÇÃO DE HIDROPATIA TOTAL E INSERINDO O VALOR NA LISTA
                                # DE CARACTERISTICAS
                                lista_caracteristicas.append(hidropatia_total(sequencias[z]))

                                # CHAMANDO A FUNÇÃO DE HIDROPATIA TOTAL E DIVIDINDO PELA QUANTIDADE
                                # DE AMINOÁCIDOS QUE COMPÕE A SEQUENCIA
                                lista_caracteristicas.append((hidropatia_total(sequencias[z]) / len(sequencias[z])))

                                # HIDROPATIA CALCULADA APENAS COM OS 25 ULTIMOS AMINOÁCIDOS DA SEQUENCIA
                                lista_caracteristicas.append(hidropatia_c_terminal(sequencias[z]))

                                # CARGA CALCULADA APENAS COM OS 25 ULTIMOS AMINOÁCIDOS DA SEQUENCIA
                                # AMINOACIDOS H R K SOMA MAIS UM NA CARGA
                                # AMINOACIDOS E D SUBTRAI UM DA CARGA
                                lista_caracteristicas.append(carga_c_terminal(sequencias[z]))

                                # QUANTIDADE CALCULADA APENAS COM OS 25 ULTIMOS AMINOÁCIDOS DA SEQUENCIA
                                # AMINOACIDOS H R K SOMA MAIS UM NA CARGA
                                # AMINOACIDOS POLARES BASICOS
                                # INDICA A QUANTIDADE DE ALCALINOS NO C TERMINAL
                                lista_caracteristicas.append(aminoacidos_basicos_c_terminal(sequencias[z]))

                                # TEM QUE ENCONTRAR NO MINIMO 2 SEQUENCIAS DE 7 AMINOACIDOS
                                # PRIMEIRO E O QUARTO HIDROFOBICOS (POSITIVOS)
                                # SEGUNDO TERCEIRO QUINTO SEXTO SETIMO HIDROFILICOS (NEGATIVOS)
                                lista_caracteristicas.append(dominio_espiral_enrolada(sequencias[z]))

                                # QUANTIDADE CALCULADA APENAS COM OS 25 ULTIMOS AMINOÁCIDOS DA SEQUENCIA
                                # AMINOACIDOS H R K SOMA MAIS UM NA CARGA
                                # AMINOACIDOS POLARES BASICOS
                                # INDICA A QUANTIDADE DE ALCALINOS NO C TERMINAL
                                lista_caracteristicas.append(sinal_localizacao_nuclear(sequencias[z]))

                                # TEM QUE SER CAPAZ DE GERAR DUPLA HELICE
                                # TEM Q TER PREDOMINANCIA DE R L S A
                                # TEM Q TER ESCASSEZ DE E D
                                lista_caracteristicas.append(sinal_localizacao_mitocondrial(sequencias[z]))

                                # TEM QUE ENCONTRAR NO MINIMO 1 SEQUENCIA COMEÇANDO COM C
                                # SEGUIDO DE I OU V OU L OU A OU P OU M
                                # SEGUIDO DE I OU V OU L OU A OU P OU M
                                # SEGUIDO DE A OU M OU S OU L
                                lista_caracteristicas.append(dominio_prenilacao(sequencias[z]))

                                # ARREDONDA TODOS OS VALORES DA LISTA PARA DUAS CASAS DECIMAIS
                                for b in range(len(lista_caracteristicas)):
                                        print(lista_caracteristicas[b])
                                        caracteristicas_arredondadas.append(round(lista_caracteristicas[b], ndigits=2))

                                caracteristicas_algoritmo = caracteristicas_arredondadas
                                # ESCREVE UMA LINHA COM TODOS OS ATRIBUTOS DA LISTA
                                writer.writerow(caracteristicas_arredondadas)
                                # ZERA AS LISTAS PARA EXTRAIR INFORMAÇÕES DE UMA NOVA PROTEÍNA
                                lista_caracteristicas[:] = []
                                #caracteristicas_arredondadas[:] = []

                        base_treinamento.close()
                        print(caracteristicas_algoritmo)
                        base = pd.read_csv("treinamento.csv")

                        #base_teste = pd.read_csv("treinamento.csv")
                        teste = caracteristicas_arredondadas
                        teste = np.asarray(teste)
                        teste = teste.reshape(1, -1)

                        previsores = base.iloc[:, 0:9].values
                        classe = base.iloc[:, 9].values

                        classificador = KNeighborsClassifier(n_neighbors=7, metric="minkowski", p=3)

                        classificador.fit(previsores, classe)

                        resultado = classificador.predict(teste)
                        veredicto = ''
                        if(resultado == 1):
                                veredicto = 'Proteína  Classificada Como Efetora!'
                        else:
                                veredicto = 'Proteína Classificada Como Não Efetora!'


        else:
                form = ContactCourse()
        return render(request, 'classificar_proteina.html', {'form': form, 'veredicto': veredicto , 'teste': teste})