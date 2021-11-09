
###Alineamiento con msa

library(msa)
secuencias<-system.file("examples", "exampleAA.fasta", package="msa") ##Estas son secuencias precargadas
secuencias<-readAAStringSet(secuencias) ##Leerlas Biostring
secuencias

##Función para alineammiento
Alineamiento1 <- msa(secuencias)

##Para ver el alineamiento.
print(Alineamiento1, show="complete")

##Revisar lo de Latex

##Permite realizar alineamineto con ClustalW, aunque por default ya lo hacia 
AlineaminetoClustalW<- msa(secuencias, "ClustalW")

#Con Omega
AlineaminetoClustalO<- msa(secuencias, "ClustalOmega")

#Con Muscle
AlineaminetoMu<- msa(secuencias, "Muscle")

##Acá, para imprimir los resultados
print(AlineaminetoMu, show="alignment")

print(AlineaminetoMu, show="alignment", showConsensus=FALSE)

print(AlineaminetoMu, show=c("alignment", "version")) 

print(AlineaminetoMu, show="standardParams") #Muestra los parámetros score de gaps

print(AlineaminetoMu, show="algParams")

print(AlineaminetoMu, show=c("call", "version"))


######################
####Ejercicios msa
######################

###Tenemos 10 secuencias concatenadas de un receptor en distintas especies

TLR2concat<-readDNAStringSet("/Users/suhanjuarez/bioinformatik/Alineamiento/TLR2seq.fasta")
TLR2concat

##Aquí les colocamos nombres, sencillos
Idnomb<-c("Duck","Dog","ZebrFish","Human","Seal","Macaca","NES","Mouse","Pig","Dolphin")
names(TLR2concat)<-Idnomb
TLR2concat

Nucleo<-alphabetFrequency(TLR2concat) #Vamoss a revisar la frecuencia de nucleótidos
Nucleo

##Gráfica general de la frecuencia de nucleótidos
boxplot(Nucleo[,1:4])
title(main="Frecuencia de nucleótidos en TLR2 de 10 especies", xlab="Base", ylab="Frecuencia de las bases")

###Para una especie
pie(Nucleo[1,1:4], col= rainbow(4)) + title(main = "Fecuencia de nucleótidos en TLR2 en el pato")

####Vamos a traducir a proteínas

TLR2AA<- translate(TLR2concat)
TLR2AA

###Vamos a ver la frecuencia de aminoácidos

Aminoácid<-alphabetFrequency(TLR2AA)
Aminoácid

##Gráfica general de la frecuencia de aminoácidos
boxplot(Aminoácid[,1:20])
title(main="Frecuencia de aminnoácidos en TLR2 de 10 especies", xlab="Aminoácido", ylab="Frecuencia")

####Alineamientos#############

###Alineamiento con ClustalW

TLR2AliClustalW<- msa(TLR2AA, "ClustalW")
print(TLR2AliClustalW, show="alignment")
print(TLR2AliClustalW, show="standardParams")

###Alineamiento con ClustalOmega

TLR2AliClustalO<- msa(TLR2AA, "ClustalOmeg")
print(TLR2AliClustalO, show="alignment")
print(TLR2AliClustalO, show="standardParams")

###Alineamiento con Muscle

TLR2AliMuscle<- msa(TLR2AA, "Muscle")  ####Revisar algo raro *
print(TLR2AliMuscle, show="alignment")
print(TLR2AliMuscle, show="standardParams")


#####Árbol filogenético#### ###REVISAR DETALLES DESDE AQUÍ, OBJETO ES PO CADA ALINEAMIENTO

TLR2tree<-msaConvert(TLR2AA, type=c("seqinr::alignment")) ##convertir el alineamiento a un objeto para hacer un árbol

##Primero se debe calcular la matriz de distancias entre las secuencias

library(seqinr) #Abrir este paquete
Distancias<-dist.alignment(TLR2tree, matrix= c("identity"),gap = FALSE)
Distancias

as.matrix(Distancias)[2:5,"Human",drop=FALSE] ##Creamos otra matriz 
#### Se extraen datos de algunos organismos respesto a uno, en este caso con Human

###Abrimos el siguuente paquete:

library(ape)

###Nos va a permitir crear un árbol

TreeTLR2<-nj(Distancias)
#La función nj realiza una estimación del árbol "neighbor-joining"

###Ahora vamos a visualizarlo
plot(TreeTLR2, main="Árbol filogenético de las secuencias de TLR2")

################################AHORA ÁRBOLES DE LOS TRES ALINEAMIENTOS

##ÁRBOL CLUSTALW
TLR2treeCLU<-msaConvert(TLR2AliClustalW, type=c("seqinr::alignment"))
DistClu<-dist.alignment(TLR2treeCLU, matrix= c("identity"),gap = FALSE)
TreeTLR2CLU<-nj(DistClu)
plot(TreeTLR2CLU, main="Árbol filogenético de las secuencias de TLR2 (Clustalw)")


##ÁRBOL MUSCLE
TLR2treeM<-msaConvert(TLR2AliMuscle, type=c("seqinr::alignment"))
DistM<-dist.alignment(TLR2treeM, matrix= c("identity"),gap = FALSE)
TreeTLR2Mu<-nj(DistM)
plot(TreeTLR2Mu, main="Árbol filogenético de las secuencias de TLR2 (Muscle)")

##ÁRBOL ClustalOmega
TLR2treeOM<-msaConvert(TLR2AliClustalO, type=c("seqinr::alignment"))
DistOM<-dist.alignment(TLR2treeOM, matrix= c("identity"),gap = FALSE)
TreeTLR2OM<-nj(DistOM)
DistOM ###Algo pasa aquí "NaN"


plot(TreeTLR2OM, main="Árbol filogenético de las secuencias de TLR2 (Clustal Omega)")
