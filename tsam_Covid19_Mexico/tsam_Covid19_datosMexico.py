import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '../')

#get_ipython().magic(u'matplotlib inline')
from tsam_Covid19_baseCode import *

urlDatosAbiertosCovid19Mexico='http://187.191.75.115/gobmx/salud/datos_abiertos/datos_abiertos_covid19.zip'
urlDiccionarioDatos = 'http://187.191.75.115/gobmx/salud/datos_abiertos/diccionario_datos_covid19.zip'
r = requests.get(diccionarioDatos)
z = zipfile.ZipFile(StringIO.StringIO(r.content))
descriptores = pandas.read_csv(z.read('Descriptores_0412.xlsx'))
catalogos = pandas.read_csv(z.read('Catalogos_0412.xlsx'))
diccionarioMexico=pd.read_csv(diccionarioDatos, compression='zip',encoding='latin-1')

datosMexico=pd.read_csv(urlDatosAbiertosCovid19Mexico, compression='zip',encoding='latin-1')
diccionarioMexico=pd.read_csv(diccionarioDatos, compression='zip',encoding='latin-1')
mexReference="""Data obtained from %s and %s"""%(datosAbiertosCovid19Mexico,diccionarioDatos)
print(mexReference)


# Extracting indices of interest
fDescriptores= '$HOME/data/Salud_GobiernoMexico/diccionario_datos_covid19/Descriptores_0412.xlsx'
fCatalogos ='$HOME/data/Salud_GobiernoMexico/diccionario_datos_covid19/Catalogos_0412.xlsx'
