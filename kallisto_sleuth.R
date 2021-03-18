#########################
BiocManager::install("rhdf5")
library(rhdf5)

install.packages("devtools")

library(devtools)

devtools::install_github("pachterlab/sleuth")

library("sleuth")

BiocManager::install("edgeR")
library(edgeR)

BiocManager::install("ensembldb")

library(ensembldb)
#Instalar y cargar las librerias necesarias para el codigo, ademas de que estuve como loco tratando de instalar Kallisto y despues
#de un rato lo logré (':
#Se puede usar la funcion setwd para cambiar el directorio donde se va a trabajar


#Esta función permite mapear, a partir de la base de datos de
tx2gene <- function(){
  
  #     Dataset you want to use. To see the different datasets available within a biomaRt yo$
  #     host
  
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = "ensembl_transcript_id",
                       ens_gene = "ensembl_gene_id", ext_gene = "external_gene_name")
  return(t2g)
}

#La funcion sirve para ver la informacion de los genes con los que se va a trabajar
t2g <- tx2gene()


#Seleccionar el directorio en el que se encuentran las muestras para que funcione y de resultados
base_dir<-"Archivo/"

#All samples
#samples <- paste0("sample_", c("14BE01A_R","14BE02A_R","14BE08A_R","17BE01A_R",
#                                "17BE02A_R","17BE03A_R"))

#Selected samples

samples <- paste0("sample", c("1","2","3","4","5","6"))  #Se seleccionan las muestras que se van a usar, con base en los nombres de las carpetas
#o archivos de la carpeta que los contiene

kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))

#All samples
#s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("2014","2014","2014", "2017",
#
#"2017","2017"), stringsAsFactors=FALSE)

#Selected samples

s2c <- data.frame(path=kal_dirs, sample=samples, muestras = c("7","7","8","8","9","9"), stringsAsFactors=FALSE)
#Aca se hace una cosa similar a la anterior, la diferencia es que son muestras de referencia, se tienen que repetir para analizar como dos de las muestras
#anteriores se comportan o expresan bajo condiciones iguales.
#Analizar como dos muestras se expresan bajo la misma condicion
#Se genera un data.frame en el cual se ingresan los resultados de los analisis, el numero de muestras debe ser el mismo que el de las muestras
#de referencia para que la matriz o data frame sea simetrica... creo

so <- sleuth_prep(s2c, ~muestras, target_mapping = t2g,extra_bootstrap_summary = TRUE)#Se compara las muestras con las muestras de referencia
#A partir de aca me salen varios warnings y errores que no impiden que el codigo funcione, no se bien como resolverlos o que
#es lo que quieren decir, pero se puede continuar corriendo el script din problema, al menos eso creo yo

so <- sleuth_fit(so) #aparece un error "fitting meassurements error models" 
so <- sleuth_wt(so, which_beta="muestras8") #en esta parte me sugiere opciones del which_beta, supongo que se ajustan bien al analisis
sleuth_live(so)
#Esta ultima parte hace un pop de una ventana llamada shiny, esta ventana ofrece herramientas de visualizacion de los datos bastante
#feten (chido jaja), salen varias opciones con las que se puede  jugar y de aqui se obtiene la test table que se usa en la siguiente parte


#Analisis disponibles para visualizacion

setwd("C:/Users/luis_/OneDrive/Documentos")
resultados<-read.table("test_table.csv",sep=",",
                       header=TRUE)
significativos<-which(resultados$qval<0.1)
significativos<-resultados[significativos,]
upregulated<-which(significativos$b>0)
upregulated<-significativos[upregulated,]
downregulated<-which(significativos$b<0)
downregulated<-significativos[downregulated,]
#Sinceramente a esta parte no le entendi muy bien, solo se que se usa la tabla de la parte anterior para generar los resultados
#y la parte de abajo es para generar archivos de texto que contienen esos resultados (?), no estoy muy seguro si son los resultados

write.table(upregulated,file="C:/Users/luis_/OneDrive/Documentos/aversisale1.txt",sep="\t")
write.table(downregulated,file="C:/Users/luis_/OneDrive/Documentos/aversisale2.txt",sep="\t")

upregulated$ext_gene
#Este ni idea de que sea :c