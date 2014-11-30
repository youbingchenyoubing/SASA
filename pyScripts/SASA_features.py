import glob
import csv
import numpy
import re

def convert(s):
    if s.isdigit(): s = int(s)
    return s

def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    # convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


numpy.seterr(invalid='raise')

'''
IZ = 182
KZ = 211
LZ = 180
MZ = 204
NZ = 158
PZ = 143
QZ = 189
RZ = 241
SZ = 122
TZ = 146
VZ = 160
WZ = 259
YZ = 229
GZ = 85
AZ = 113
CZ = 140
DZ = 151
EZ = 183
FZ = 218
HZ = 194
'''

ILE = 182
LYS = 211
LEU = 180
MET = 204
ASN = 158
PRO = 143
GLN = 189
ARG = 241
SER = 122
THR = 146
VAL = 160
TRP = 259
TYR = 229
GLY = 85
ALA = 113
CYS = 140
ASP = 151
GLU = 183
PHE = 218
HIS = 194
HID = 195
HIE = 195
HIP = 196
CYM = 139
masterlist = []
masterlist2 = []
cabecalho = ['Res. num','Res. name']
cabecalho1 = ['','']
cabecalho2 =  []
ficheiro1 = open('result.csv', 'wb')
ficheiro2 = open('medias.csv', 'wb')
result = csv.writer(ficheiro1, 'excel')
result2 = csv.writer(ficheiro2, 'excel')
#result.writerow(('Numero','Residuo'))
nomes = glob.glob('res_sasa'+'*.dat')
nomes1 = sorted_nicely(nomes)
last_row = []
#nomes.sort()
for x in nomes1:
    if x not in glob.glob('res_sasa'+'A*.dat') and x not in glob.glob('res_sasa'+'B*.dat'):
        cabecalho.append(x)
        cabecalho.append(' ')
        cabecalho.append(' ')
        cabecalho.append(' ')
        cabecalho.append(' ')
        cabecalho.append(' ')
        cabecalho.append(' ')
        cabecalho.append(' ')
        cabecalho.append(' ')
        cabecalho.append(' ')
        cabecalho.append(' ')
        cabecalho.append(' ')
for x in nomes1:
    if x not in glob.glob('res_sasa'+'A*.dat') and x not in glob.glob('res_sasa'+'B*.dat'):
        cabecalho1.append('SASA comp')
        cabecalho1.append('SASA mon')
        cabecalho1.append('SASA del')
        cabecalho1.append('SASA rel')
        cabecalho1.append('SASA comp/res')
        cabecalho1.append('SASA mon/res')
        cabecalho1.append('SASA del/res')
        cabecalho1.append('SASA rel/res')
        cabecalho1.append('SASA comp/ave')
        cabecalho1.append('SASA mon/ave')
        cabecalho1.append('SASA del/ave')
        cabecalho1.append('SASA rel/ave')

cabecalho2.append('SASA comp')
cabecalho2.append(' ')
cabecalho2.append('SASA mon')
cabecalho2.append(' ')
cabecalho2.append('SASA del')
cabecalho2.append(' ')
cabecalho2.append('SASA rel')
cabecalho2.append(' ')
cabecalho2.append('SASA comp/res')
cabecalho2.append(' ')
cabecalho2.append('SASA mon/res')
cabecalho2.append(' ')
cabecalho2.append('SASA del/res')
cabecalho2.append(' ')
cabecalho2.append('SASA rel/res')
cabecalho2.append(' ')
cabecalho2.append('SASA comp/ave')
cabecalho2.append(' ')
cabecalho2.append('SASA mon/ave')
cabecalho2.append(' ')
cabecalho2.append('SASA del/ave')
cabecalho2.append(' ')
cabecalho2.append('SASA rel/ave')
cabecalho2.append(' ')



cabecalho.append(' ')
cabecalho.append('Medias')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append('Deviations')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')
cabecalho.append(' ')

cabecalho1.append(' ')
cabecalho1.append('SASA comp')
cabecalho1.append('SASA mon')
cabecalho1.append('SASA del')
cabecalho1.append('SASA rel')
cabecalho1.append('SASA comp/res')
cabecalho1.append('SASA mon/res')
cabecalho1.append('SASA del/res')
cabecalho1.append('SASA rel/res')
cabecalho1.append('SASA comp/ave')
cabecalho1.append('SASA mon/ave')
cabecalho1.append('SASA del/ave')
cabecalho1.append('SASA rel/ave')
cabecalho1.append(' ')
cabecalho1.append('SASA comp')
cabecalho1.append('SASA mon')
cabecalho1.append('SASA del')
cabecalho1.append('SASA rel')
cabecalho1.append('SASA comp/res')
cabecalho1.append('SASA mon/res')
cabecalho1.append('SASA del/res')
cabecalho1.append('SASA rel/res')
cabecalho1.append('SASA comp/ave')
cabecalho1.append('SASA mon/ave')
cabecalho1.append('SASA del/ave')
cabecalho1.append('SASA rel/ave')

ave_ILE = 0
ave_LYS = 0
ave_LEU = 0
ave_MET = 0
ave_ASN = 0
ave_PRO = 0
ave_GLN = 0
ave_ARG = 0
ave_SER = 0
ave_THR = 0
ave_VAL = 0
ave_TRP = 0
ave_TYR = 0
ave_GLY = 0
ave_ALA = 0
ave_CYS = 0
ave_ASP = 0
ave_GLU = 0
ave_PHE = 0
ave_HIS = 0
ave_HID = 0
ave_HIE = 0
ave_HIP = 0
ave_CYM = 0

temp_ILE = []
temp_LYS = []
temp_LEU = []
temp_MET = []
temp_ASN = []
temp_PRO = []
temp_GLN = []
temp_ARG = []
temp_SER = []
temp_THR = []
temp_VAL = []
temp_TRP = []
temp_TYR = []
temp_GLY = []
temp_ALA = []
temp_CYS = []
temp_ASP = []
temp_GLU = []
temp_PHE = []
temp_HIS = []
temp_HID = []
temp_HIE = []
temp_HIP = []
temp_CYM = []

for x in nomes1:
    if x not in glob.glob('res_sasa'+'A*.dat') and x not in glob.glob('res_sasa'+'B*.dat'):
        temp_ILE = []
        temp_LYS = []
        temp_LEU = []
        temp_MET = []
        temp_ASN = []
        temp_PRO = []
        temp_GLN = []
        temp_ARG = []
        temp_SER = []
        temp_THR = []
        temp_VAL = []
        temp_TRP = []
        temp_TYR = []
        temp_GLY = []
        temp_ALA = []
        temp_CYS = []
        temp_ASP = []
        temp_GLU = []
        temp_PHE = []
        temp_HIS = []
        temp_HID = []
        temp_HIE = []
        temp_HIP = []
        temp_CYM = []
        print x
        acomp = []
        bcomp = []
        ccomp = []
        compfile = open(x)
        for line in compfile:
            b = line.strip('\n').split(' ')
            acomp.append(str(int(b[0])+1))
            bcomp.append(b[1])
            ccomp.append(float(b[2]))
            try:
                globals()['temp_'+b[1]].append(float(b[2]))
            except KeyError:
                pass
        compfile.close()
        z = x[0:8]+'A'+x[8:]
        print z
        amon = []
        bmon = []
        cmon = []
        monfile = open(z)
        for line in monfile:
            b = line.strip('\n').split(' ')
            amon.append(str(int(b[0])+1))
            bmon.append(b[1])
            cmon.append(float(b[2]))
        monfile.close()
        y = x[0:8]+'B'+x[8:]
        print y
        monfile2 = open(y)
        for line in monfile2:
            b = line.strip('\n').split(' ')
            amon.append(str(int(len(amon))+1))
            bmon.append(b[1])
            cmon.append(float(b[2]))
        monfile2.close()

        try:
            if len(temp_ILE)!=0:
                ave_ILE = numpy.average(temp_ILE)
        except FloatingPointError:
            print temp_ILE
            print 'temp_ILE'

        try:
            if len(temp_LYS)!=0:
                ave_LYS = numpy.average(temp_LYS)
        except FloatingPointError:
            print temp_LYS
            print 'temp_LYS'

        try:
            if len(temp_LEU)!=0:
                ave_LEU = numpy.average(temp_LEU)
        except FloatingPointError:
            print temp_LEU
            print 'temp_LEU'

        try:
            if len(temp_MET)!=0:
                ave_MET = numpy.average(temp_MET)
        except FloatingPointError:
            print temp_MET
            print 'temp_MET'

        try:
            if len(temp_ASN)!=0:
                ave_ASN = numpy.average(temp_ASN)
        except FloatingPointError:
            print temp_ASM
            print 'temp_ASM'

        try:
            if len(temp_PRO)!=0:
                ave_PRO = numpy.average(temp_PRO)
        except FloatingPointError:
            print temp_PRO
            print 'temp_PRO'

        try:
            if len(temp_GLN)!=0:
                ave_GLN = numpy.average(temp_GLN)
        except FloatingPointError:
            print temp_GLN
            print 'temp_GLN'

        try:
            if len(temp_ARG)!=0:
                ave_ARG = numpy.average(temp_ARG)
        except FloatingPointError:
            print temp_ARG
            print 'temp_ARG'

        try:
            if len(temp_SER)!=0:
                ave_SER = numpy.average(temp_SER)
        except FloatingPointError:
            print temp_SER
            print 'temp_SER'

        try:
            if len(temp_THR)!=0:
                ave_THR = numpy.average(temp_THR)
        except FloatingPointError:
            print temp_THR
            print 'temp_THR'

        try:
            if len(temp_VAL)!=0:
                ave_VAL = numpy.average(temp_VAL)
        except FloatingPointError:
            print temp_VAL
            print 'temp_VAL'

        try:
            if len(temp_TRP)!=0:
                ave_TRP = numpy.average(temp_TRP)
        except FloatingPointError:
            print temp_TRP
            print 'temp_TRP'

        try:
            if len(temp_TYR)!=0:
                ave_TYR = numpy.average(temp_TYR)
        except FloatingPointError:
            print temp_TYR
            print 'temp_TYR'

        try:
            if len(temp_GLY)!=0:
                ave_GLY = numpy.average(temp_GLY)
        except FloatingPointError:
            print temp_GLY
            print 'temp_GLY'

        try:
            if len(temp_ALA)!=0:
                ave_ALA = numpy.average(temp_ALA)
        except FloatingPointError:
            print temp_ALA
            print 'temp_ALA'

        try:
            if len(temp_CYS)!=0:
                ave_CYS = numpy.average(temp_CYS)
        except FloatingPointError:
            print temp_CYS
            print 'temp_CYS'

        try:
            if len(temp_ASP)!=0:
                ave_ASP = numpy.average(temp_ASP)
        except FloatingPointError:
            print temp_ASP
            print 'temp_ASP'

        try:
            if len(temp_GLU)!=0:
                ave_GLU = numpy.average(temp_GLU)
        except FloatingPointError:
            print temp_GLU
            print 'temp_GLU'

        try:
            if len(temp_PHE)!=0:
                ave_PHE = numpy.average(temp_PHE)
        except FloatingPointError:
            print temp_PHE
            print 'temp_PHE'

        try:
            if len(temp_HIS)!=0:    
                ave_HIS = numpy.average(temp_HIS)
        except FloatingPointError:
            print temp_HIS
            print 'temp_HIS'

        try:
            if len(temp_HID)!=0:
                ave_HID = numpy.average(temp_HID)
        except FloatingPointError:
            print temp_HID
            print 'temp_HID'

        try:
            if len(temp_HIE)!=0:
                ave_HIE = numpy.average(temp_HIE)
        except FloatingPointError:
            print temp_HIE
            print 'temp_HIE'

        try:
            if len(temp_HIP)!=0:
                ave_HIP = numpy.average(temp_HIP)
        except FloatingPointError:
            print temp_HIP
            print 'temp_HIP'

        try:
            if len(temp_CYM)!=0:
                ave_CYM = numpy.average(temp_CYM)
        except FloatingPointError:
            print temp_CYM
            print 'temp_CYM'
        print len(acomp), len(amon)
        delsasa = []
        relsasa = []
        compdiv = []
        mondiv = []
        sasadiv = []
        reldiv = []
        compdiv1 = []
        mondiv1 = []
        sasadiv1 = []
        reldiv1 = []
        for h in range(len(acomp)):
            #delsasa.append(str(float(ccomp[h])-float(cmon[h])))
            delsasa.append(float(ccomp[h])-float(cmon[h]))
            if float(cmon[h]) != 0:
                relsasa.append((float(ccomp[h])-float(cmon[h]))/(float(cmon[h])))
            else:
                relsasa.append(float('0'))

            try:
                compdiv.append(float(ccomp[h])/globals()[bcomp[h]])
            except KeyError:
                compdiv.append(float('0'))
            except ZeroDivisionError:
                compdiv.append(float('0'))
            try:
                mondiv.append(float(cmon[h])/globals()[bcomp[h]])
            except KeyError:
                mondiv.append(float('0'))
            except ZeroDivisionError:
                mondiv.append(float('0'))
            try:
                sasadiv.append(float(delsasa[h])/globals()[bcomp[h]])
            except KeyError:
                sasadiv.append(float('0'))
            except ZeroDivisionError:
                sasadiv.append(float('0'))
            try:
                reldiv.append(float(relsasa[h])/globals()[bcomp[h]]*1000)
            except KeyError:
                reldiv.append(float('0'))
            except ZeroDivisionError:
                reldiv.append(float('0'))
            try:
                compdiv1.append(float(ccomp[h])/globals()['ave_'+bcomp[h]])
            except KeyError:
                compdiv1.append(float('0'))
            except ZeroDivisionError:
                compdiv1.append(float('0'))
            except FloatingPointError:
                if ccomp[h] == 0 and globals()['ave_'+bcomp[h]] == 0:
                    compdiv1.append(float('0'))
            try:
                mondiv1.append(float(cmon[h])/globals()['ave_'+bcomp[h]])
            except KeyError:
                mondiv1.append(float('0'))
            except ZeroDivisionError:
                mondiv1.append(float('0'))
            except FloatingPointError:
                if ccomp[h] == 0 and globals()['ave_'+bcomp[h]] == 0:
                    mondiv1.append(float('0'))
            try:
                sasadiv1.append(float(delsasa[h])/globals()['ave_'+bcomp[h]])
            except KeyError:
                sasadiv1.append(float('0'))
            except ZeroDivisionError:
                sasadiv1.append(float('0'))
            except FloatingPointError:
                if ccomp[h] == 0 and globals()['ave_'+bcomp[h]] == 0:
                    sasadiv1.append(float('0'))
            try:
                reldiv1.append(float(relsasa[h])/globals()['ave_'+bcomp[h]]*1000)
            except KeyError:
                reldiv1.append(float('0'))
            except ZeroDivisionError:
                reldiv1.append(float('0'))
            except FloatingPointError:
                if ccomp[h] == 0 and globals()['ave_'+bcomp[h]] == 0:
                    reldiv1.append(float('0'))

        
        print len(acomp)
        ave10 = ' '
        ave11 = ' '
        #acomp.append(' ')
        #bcomp.append(' ')
        ave = numpy.average(ccomp)
        #ccomp.append(ave)
        ave1 = numpy.average(cmon)
        #cmon.append(ave)
        ave2 = numpy.average(delsasa)
        #delsasa.append(ave)
        ave3 = numpy.average(relsasa)
        #relsasa.append(ave)
        ave4 = numpy.average(compdiv)
        #compdiv.append(ave)
        ave5 = numpy.average(mondiv)
        #mondiv.append(ave)
        ave6 = numpy.average(sasadiv)
        #sasadiv.append(ave)
        ave12 = numpy.average(reldiv)
        #reldiv.append(ave)
        ave7 = numpy.average(compdiv1)
        #ave7 = 'ERRO'
        #compdiv1.append(ave)
        ave8 = numpy.average(mondiv1)
        #ave8 = 'ERRO'
        #mondiv1.append(ave)
        ave9 = numpy.average(sasadiv1)
        #ave9 = 'ERRO'
        #sasadiv1.append(ave)
        ave13 = numpy.average(reldiv1)
        #ave13 = 'ERRO'
        #reldiv1.append(ave)
        if x == 'res_sasa1.dat':
            last_row.append(ave10)
            last_row.append(ave11)
            last_row.append(ave)
            last_row.append(ave1)
            last_row.append(ave2)
            last_row.append(ave3)
            last_row.append(ave4)
            last_row.append(ave5)
            last_row.append(ave6)
            last_row.append(ave12)
            last_row.append(ave7)
            last_row.append(ave8)
            last_row.append(ave9)
            last_row.append(ave13)
        else:
            last_row.append(ave)
            last_row.append(ave1)
            last_row.append(ave2)
            last_row.append(ave3)
            last_row.append(ave4)
            last_row.append(ave5)
            last_row.append(ave6)
            last_row.append(ave12)
            last_row.append(ave7)
            last_row.append(ave8)
            last_row.append(ave9)
            last_row.append(ave13)
        if x == 'res_sasa1.dat':
            masterlist.append(acomp)
            masterlist.append(bcomp)
            masterlist.append(ccomp)
            masterlist.append(cmon)
            masterlist.append(delsasa)
            masterlist.append(relsasa)
            masterlist.append(compdiv)
            masterlist.append(mondiv)
            masterlist.append(sasadiv)
            masterlist.append(reldiv)
            masterlist.append(compdiv1)
            masterlist.append(mondiv1)
            masterlist.append(sasadiv1)
            masterlist.append(reldiv1)
        else:
            masterlist.append(ccomp)
            masterlist.append(cmon)
            masterlist.append(delsasa)
            masterlist.append(relsasa)
            masterlist.append(compdiv)
            masterlist.append(mondiv)
            masterlist.append(sasadiv)
            masterlist.append(reldiv)
            masterlist.append(compdiv1)
            masterlist.append(mondiv1)
            masterlist.append(sasadiv1)
            masterlist.append(reldiv1)
sep = []
for i in range(len(masterlist[0])):
    sep.append('')

masterlist.append(sep)

med1 = []
med2 = []
med3 = []
med4 = []
med5 = []
med6 = []
med7 = []
med8 = []
med9 = []
med10 = []
med11 = []
med12 = []

std1 = []
std2 = []
std3 = []
std4 = []
std5 = []
std6 = []
std7 = []
std8 = []
std9 = []
std10 = []
std11 = []
std12 = []

#print masterlist[0]

for i in range(len(masterlist[0])):
    #for x in range(len(masterlist[:-1])/10):
    #for x in range(len(masterlist)/10):
    #    temp_med.append(masterlist[(10*x)+2][i])
    #print i, masterlist[0][i]
    for l in range(12):
        vals_med = []
        #for x in range(len(masterlist)/12):
        for x in range(len(masterlist[:-2])/12):
            #print x
            vals_med.append(masterlist[(12*x)+l+2][i])
            #vals_med.append('ERRO')
            #vals_med.append(0)
        temp_med = numpy.average(vals_med)
        #temp_med = 'ERRO'
        temp_std = numpy.std(vals_med)
        #temp_std = 'ERRO'
        globals()['med'+str(l+1)].append(temp_med)
        globals()['std'+str(l+1)].append(temp_std)


masterlist.append(med1)
masterlist.append(med2)
masterlist.append(med3)
masterlist.append(med4)
masterlist.append(med5)
masterlist.append(med6)
masterlist.append(med7)
masterlist.append(med8)
masterlist.append(med9)
masterlist.append(med10)
masterlist.append(med11)
masterlist.append(med12)
masterlist.append(sep)
masterlist.append(std1)
masterlist.append(std2)
masterlist.append(std3)
masterlist.append(std4)
masterlist.append(std5)
masterlist.append(std6)
masterlist.append(std7)
masterlist.append(std8)
masterlist.append(std9)
masterlist.append(std10)
masterlist.append(std11)
masterlist.append(std12)

masterlist2.append(med1)
masterlist2.append(std1)
masterlist2.append(med2)
masterlist2.append(std2)
masterlist2.append(med3)
masterlist2.append(std3)
masterlist2.append(med4)
masterlist2.append(std4)
masterlist2.append(med5)
masterlist2.append(std5)
masterlist2.append(med6)
masterlist2.append(std6)
masterlist2.append(med7)
masterlist2.append(std7)
masterlist2.append(med8)
masterlist2.append(std8)
masterlist2.append(med9)
masterlist2.append(std9)
masterlist2.append(med10)
masterlist2.append(std10)
masterlist2.append(med11)
masterlist2.append(std11)
masterlist2.append(med12)
masterlist2.append(std12)


result.writerow(cabecalho)
result.writerow(cabecalho1)
result2.writerow(cabecalho2)
for i in range(len(masterlist[0])):
    temp_row = []
    for x in range(len(masterlist)):
        try:
            temp_row.append(masterlist[x][i])
        except IndexError:
            temp_row.append('NADA')
    result.writerow(temp_row)

#print len(masterlist2[0])

for i in range(len(masterlist2[0])):
    temp_row2 = []
    for x in range(len(masterlist2)):
        try:
            temp_row2.append(masterlist2[x][i])
        except IndexError:
            temp_row2.append('NADA')
    result2.writerow(temp_row2)

result.writerow(last_row)
ficheiro1.close()
ficheiro2.close()
