import mpcorbfile
import pathlib
from tqdm import tqdm
import numpy as np
from scipy.optimize import minimize
from astropy.constants import GM_sun, au
from astroquery.jplhorizons import Horizons


def D_Nesvorny(a1,a2,e1,e2,i1,i2,O1,O2,vp1,vp2):
    ka = 5./4.
    ke = 2.
    ki = 2.
    kO = 1.e-6
    kvp = 1.e-6

    n = np.sqrt(GM_sun.value/a1**3)
    D2  = ka*((a1-a2)/a1)**2 + ke*(e1-e2)**2 + ki*(np.sin(i1)-np.sin(i2))**2 \
        + kO*(O1-O2) + kvp*(vp1-vp2)**2
    D = np.sqrt(D2)*n*a1

    return D

def Drummond_criterion(q1,q2,e1,e2,i1,i2,O1,O2,o1,o2):

    def comp_beta(i,omega):
        return np.arcsin(np.sin(i)*np.sin(omega))

    def comp_lambda(i,Omega,omega):
        lamb = Omega + np.arctan(np.cos(i)*np.tan(omega))
        if np.cos(omega) < 0:
            lamb = lamb + np.pi

        return lamb

    beta_1 = comp_beta(i1,o1)
    beta_2 = comp_beta(i2,o2)

    lamb_1 = comp_lambda(i1,O1,o1)
    lamb_2 = comp_lambda(i2,O2,o2)

    I = np.arccos(np.cos(i1)*np.cos(i2) \
        + np.sin(i1)*np.sin(i2)*np.cos(O1-O2))

    cos_lamb_12 = np.cos(lamb_1)*np.cos(lamb_2) + np.sin(lamb_1)*np.sin(lamb_2)

    THETA = np.arccos(
        np.sin(beta_1)*np.sin(beta_2) +
        np.cos(beta_1)*np.cos(beta_2)*cos_lamb_12
    )


    D2 = ((e1 - e2)/(e1 + e2))**2 + ((q1 - q2)/(q1 + q2))**2 \
        + (I/np.pi)**2 + (((e1 + e2)/2.)*(THETA/np.pi))**2

    D = np.sqrt(D2)

    return D




#
#
#
# def Drummond_criterion(eA,eB,qA,qB,iA,iB,OA,OB,oA,oB):
#
#     def comp_beta(i,omega):
#         return np.arcsin(np.sin(i)*np.sin(omega))
#
#     def comp_lambda(i,Omega,omega):
#         lamb = Omega + np.arctan(np.cos(i)*np.tan(omega))
#         if np.cos(omega) < 0:
#             lamb = lamb + np.pi
#
#         return lamb
#
#     beta_A = comp_beta(iA,oA)
#     beta_B = comp_beta(iB,oB)
#
#     lamb_A = comp_lambda(iA,OA,oA)
#     lamb_B = comp_lambda(iB,OB,oB)
#
#     I_BA = np.arccos(np.cos(iA)*np.cos(iB) \
#         + np.sin(iA)*np.sin(iB)*np.cos(OA-OB))
#
#     cos_lamb_BA = np.cos(lamb_A)*np.cos(lamb_B) + np.sin(lamb_A)*np.sin(lamb_B)
#
#     theta_BA = np.arccos(
#         np.sin(beta_A)*np.sin(beta_B) +
#         np.cos(beta_A)*np.cos(beta_B)*cos_lamb_BA
#     )
#
#
#     D2 = ((qB - qA)/(qB + qA))**2 + ((eB - eA)/(eB + eA))**2 \
#         + (I_BA/np.pi)**2 + (((eB + eA)/2.)*(theta_BA/np.pi))**2
#
#     D = np.sqrt(D2)
#
#     return D





# READ FILE
#mpc = mpcorbfile.mpcorb_file('MPCORB.DAT')
mpc = mpcorbfile.mpcorb_file('input/MPCORB-11-MAY.DAT')

N = len(mpc.bodies)

dtype = [('name', 'U20'), ('index', 'i4'), ('Drummond', 'f4')]
data = np.zeros((3,N))
id_ast = ['' for _ in range(N)]

FB =[1.337640, 0.284638, 26.195974, 338.578091, 550.646025]
a2 = FB[0]
q2 = FB[0]*(1.-FB[1])
e2 = FB[1]
i2 = np.radians(FB[2])%(2.*np.pi)
O2 = np.radians(FB[3])%(2.*np.pi)
o2 = np.radians(FB[4])%(2.*np.pi)
vp2 = (O2 + o2)%(2.*np.pi)

for i in tqdm(range(N)):
    a1 = mpc.bodies[i]['a']
    e1 = mpc.bodies[i]['e']
    q1 = a1*(1.-e1)
    i1 = np.radians(mpc.bodies[i]['i'])%(2.*np.pi)
    O1 = np.radians(mpc.bodies[i]['Node'])%(2.*np.pi)
    o1 = np.radians(mpc.bodies[i]['Peri'])%(2.*np.pi)
    vp1 = (O1 + o1)%(2.*np.pi)

    if i < 793066:
        id_ast[i] = mpc.bodies[i]['packed_designation']
    else:
        id_ast[i] = mpc.bodies[i]['designation']

    data[0, i] = i

    data[1,i] = Drummond_criterion(q1,q2,e1,e2,i1,i2,O1,O2,o1,o2)
    data[2,i] = D_Nesvorny(a1,a2,e1,e2,i1,i2,O1,O2,vp1,vp2)



a1 = mpc.bodies[i]['a']
e1 = mpc.bodies[i]['e']
q1 = a1*(1.-e1)
i1 = np.radians(mpc.bodies[i]['i'])%(2.*np.pi)
O1 = np.radians(mpc.bodies[i]['Node'])%(2.*np.pi)
o1 = np.radians(mpc.bodies[i]['Peri'])%(2.*np.pi)
vp1 = (O1 + o1)%(2.*np.pi)




print(Drummond_criterion(q1,q2,e1,e2,i1,i2,O1,O2,o1,o2))






sort_D = np.argsort(data[1, :])
sort_N = np.argsort(data[2, :])
Dru = data[:, sort_D]
Ner = data[:, sort_N]


redo_N = 100
rd_Dru = np.zeros((3,redo_N))
rd_Ner = np.zeros((3,redo_N))


jd_date = 2460003.2604167
for i in tqdm(range(redo_N)):
    rd_Dru[0,i] = Dru[0,i]
    rd_Dru[1,i] = Dru[1,i]

    obj = Horizons(id=f'{id_ast[int(rd_Dru[0,i])]}', id_type='smallbody', location='@sun', epochs=jd_date)
    try:
        elements = obj.elements()
    except ValueError:
        obj = Horizons(id=f'{int(rd_Dru[0,i]+1)}', id_type='smallbody', location='@sun', epochs=jd_date)
        elements = obj.elements()

    #print(f'{id_ast[int(rd_Dru[0,i])],int(rd_Dru[0,i])}', mpc.bodies[int(rd_Dru[0,i])]['packed_designation'], mpc.bodies[int(rd_Dru[0,i])]['designation'])


    a1 = elements["a"].value[0]
    e1 = elements["e"].value[0]
    q1 = a1*(1.-e1)
    i1 = np.radians(elements["incl"].value[0])%(2.*np.pi)
    O1 = np.radians(elements["Omega"].value[0])%(2.*np.pi)
    o1 = np.radians(elements["w"].value[0])%(2.*np.pi)
    vp1 = (O1 + o1)%(2.*np.pi)
    rd_Dru[2,i] = Drummond_criterion(q1,q2,e1,e2,i1,i2,O1,O2,o1,o2)



obj = Horizons(id=f'K18P10L', id_type='smallbody', location='@sun', epochs=jd_date)
elements = obj.elements()
a1 = elements["a"].value[0]
e1 = elements["e"].value[0]
q1 = a1*(1.-e1)
i1 = np.radians(elements["incl"].value[0])%(2.*np.pi)
O1 = np.radians(elements["Omega"].value[0])%(2.*np.pi)
o1 = np.radians(elements["w"].value[0])%(2.*np.pi)
vp1 = (O1 + o1)%(2.*np.pi)
print(Drummond_criterion(q1,q2,e1,e2,i1,i2,O1,O2,o1,o2))

sort_D = np.argsort(rd_Dru[2, :])

Dru_N = rd_Dru[:, sort_D]



for i in range(10):
    print(i, mpc.bodies[int(Dru[0,i])]['Name'], Dru[1,i], mpc.bodies[int(Dru_N[0,i])]['Name'], Dru_N[1,i],Dru_N[2,i])



