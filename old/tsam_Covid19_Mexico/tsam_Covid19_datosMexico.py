import sys, zipfile
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '../')

#get_ipython().magic(u'matplotlib inline')
from tsam_Covid19_baseCode import *

localSource = '/Users/curandero/data/Salud_GobiernoMexico/diccionario_datos_covid19/'
urlSource = 'http://187.191.75.115/gobmx/salud/datos_abiertos/'
urlDatos=urlSource+'datos_abiertos_covid19.zip'
urlDiccionario = urlSource+'diccionario_datos_covid19.zip'
fDescriptores= localSource+'Descriptores_0412.xlsx'
fCatalogos = localSource+'Catalogos_0412.xlsx'

datosMexico=pd.read_csv(urlDatos, compression='zip',encoding='latin-1')
descriptores = pd.read_excel(fDescriptores)
entidades = pd.read_excel(fCatalogos,sheet_name=['Catálogo de ENTIDADES'])
mexReference="""Data obtained from %s and %s"""%(datosAbiertosCovid19Mexico,diccionarioDatos)
print(mexReference)

if localData>0:
    fDescriptores= '$HOME/data/Salud_GobiernoMexico/diccionario_datos_covid19/Descriptores_0412.xlsx'
    fCatalogos ='$HOME/data/Salud_GobiernoMexico/diccionario_datos_covid19/Catalogos_0412.xlsx'
    fDatos ='$HOME/data/Salud_GobiernoMexico/datos_abiertos_covid1920200418.csv'
    descriptores = pd.read_csv(z.read('Descriptores_0412.xlsx'))
    catalogos = pd.read_csv(z.read('Catalogos_0412.xlsx'))
    datosMexico=pd.read_csv(fDatos, compression='zip',encoding='latin-1')

estados = entidades['Catálogo de ENTIDADES']
estadosAb = entidades.iloc['ABREVIATURA'].to_numpy()
estadosNombre = entidades.iloc['ENTIDAD_FEDERATIVA'].to_numpy()
estadosClave = entidades.iloc['CLAVE_ENTIDAD'].to_numpy()

descInt =
