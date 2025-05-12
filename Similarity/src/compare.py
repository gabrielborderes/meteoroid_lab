import mpcorbfile
#import rebound
import pathlib
from tqdm import tqdm
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from astropy.constants import GM_sun, au
from astroquery.jplhorizons import Horizons
import re
import os


def D_NesvornyDaniela(a1,a2,e1,e2,i1,i2,O1,O2,vp1,vp2):


    n = np.sqrt(GM_sun.value/a1**3)
    deltaA2 = (abs(a2 - a1)/a1)**2
    deltaE2 = (abs(e2 - e1))**2
    deltaI2 = (abs(np.sin(np.radians(i2)) - np.sin(np.radians(i1))))**2
    deltaw2 = (abs(vp2 - vp1))**2
    deltaW2 = (abs((O2) - (O1)))**2
    ka,ke,ki,kw,kW = 5/4,2,2,1.e-6,1.e-6
    dna = deltaA2*ka + deltaE2*ke + deltaI2*ki + deltaw2*kw + deltaW2*kW
    d = np.sqrt(dna)*n*a1
    return d



def D_Nesvorny(a1,a2,e1,e2,i1,i2,O1,O2,vp1,vp2):
    ka = 5./4.
    ke = 2.
    ki = 2.
    kO = 1.e-6
    kvp = 1.e-6

    n = np.sqrt(GM_sun.value/a1**3)
    A  = ka*((a1-a2)/a1)**2 + ke*(e1-e2)**2 + ki*(np.sin(i1)-np.sin(i2))**2 \
        + kO*(O1-O2) + kvp*(vp1-vp2)**2
    D2 = A*(n*a1)**2
    D = np.sqrt(D2)

    return D




def D_criterion(eA,eB,qA,qB,iA,iB,OA,OB,oA,oB):

    half_piece_1 = 2.*(np.sin(.5*iA)*np.cos(.5*iB)-np.cos(.5*iA)*np.sin(.5*iB))
    half_piece_1 = half_piece_1**2

    two_sim_O = 2.*(np.sin(.5*OA)*np.cos(.5*OB)-np.cos(.5*OA)*np.sin(.5*OB))
    quarter_piece_1 = np.sin(iA)*np.sin(iB)*two_sim_O**2

    OoA = OA + oA
    OoB = OB + oB
    two_sim_Oo = 2.*(np.sin(.5*OoA)*np.cos(.5*OoB)-np.cos(.5*OoA)*np.sin(.5*OoB))
    quarter_piece_2 = (((eA+eB)*.5)*two_sim_Oo)**2

    half_piece_2 = quarter_piece_1 + quarter_piece_2

    piece = half_piece_1 + half_piece_2
    D2 = ((qA - qB)/(au.value*1.e-3))**2 + (eA - eB)**2 + piece

    D = np.sqrt(D2)

    return D


def Drummond_criterion(eA,eB,qA,qB,iA,iB,OA,OB,oA,oB):

    def comp_beta(i,omega):
        return np.arcsin(np.sin(i)*np.sin(omega))

    def comp_lambda(i,Omega,omega):
        lamb = Omega + np.arctan(np.cos(i)*np.tan(omega))
        if np.cos(omega) < 0:
            lamb = lamb + np.pi

        return lamb

    beta_A = comp_beta(iA,oA)
    beta_B = comp_beta(iB,oB)

    lamb_A = comp_lambda(iA,OA,oA)
    lamb_B = comp_lambda(iB,OB,oB)

    I_BA = np.arccos(np.cos(iA)*np.cos(iB) \
        + np.sin(iA)*np.sin(iB)*np.cos(OA-OB))

    cos_lamb_BA = np.cos(lamb_A)*np.cos(lamb_B) + np.sin(lamb_A)*np.sin(lamb_B)

    theta_BA = np.arccos(
        np.sin(beta_A)*np.sin(beta_B) +
        np.cos(beta_A)*np.cos(beta_B)*cos_lamb_BA
    )


    D2 = ((qB - qA)/(qB + qA))**2 + ((eB - eA)/(eB + eA))**2 \
        + (I_BA/np.pi)**2 + (((eB + eA)/2.)*(theta_BA/np.pi))**2

    D = np.sqrt(D2)
    #print("EU>>>",((eB - eA)/(eB + eA))**2, ((qB - qA)/(qB + qA))**2,(I_BA/np.pi)**2,(((eB + eA)/2.)*(theta_BA/np.pi))**2)
    #print("EU>>>",D)
    return D


def calcularD_Drummond(mede,se,q2,q1,medI,sI,medW,sW,medw,sw):


  lambda1 = sW + (np.arctan(np.cos(sI)*np.tan(sw)))
  lambda1 = lambda1 + np.pi if np.cos(sw) < 0 else lambda1
  lambda2 = medW + (np.arctan(np.cos(medI)*np.tan(medw)))
  lambda2 = lambda2 + np.pi if np.cos(medw) < 0 else lambda2
  beta1 =  (np.arcsin(np.sin(sI)*np.sin(sw)))
  beta2 =  (np.arcsin(np.sin(medI)*np.sin(medw)))
  theta1 = np.arccos(np.sin(beta1)*np.sin(beta2)+np.cos(beta1)*np.cos(beta2)*np.cos(lambda1-lambda2))
  incG1 = np.arccos(np.cos(sI)*np.cos(medI)+np.sin(sI)*np.sin(medI)*np.cos(sW - medW))

  deltaE2 = ((se - mede )/(mede + se))**2
  deltaQ2 = ((q1-q2)/(q1+q2))**2
  I2 = (np.degrees(incG1)/180.)**2
  f2 =(((mede+se)/2)*(np.degrees(theta1)/180))**2

  d = np.sqrt(deltaE2+deltaQ2+I2+f2)
  #print("Daniela>>>",deltaE2, deltaQ2,I2,f2,incG1)
  #print("Daniela>>>",d)
  return d


def OE_2_position(a, e, i, Omega, omega, f):

    r = a * (1 - e**2) / (1 + e * np.cos(f))
    varpi = omega + f

    x = r * (np.cos(Omega) * np.cos(varpi) - np.sin(Omega) * np.sin(varpi) * np.cos(i))
    y = r * (np.sin(Omega) * np.cos(varpi) + np.cos(Omega) * np.sin(varpi) * np.cos(i))
    z = r * (np.sin(varpi) * np.sin(i))

    return np.array([x, y, z])


def func_dist(f, a, e, i, Omega, omega, EP):
    POO = OE_2_position(a, e, i, Omega, omega, f)
    dist = np.sqrt((POO[0] - EP[0]) ** 2 + (POO[1] - EP[1]) ** 2 + (POO[2] - EP[2]) ** 2)

    return dist

def read_meteor(file_in):
    with open(file_in, 'r') as file:
            content = file.read()

            jd_match = re.search(r'JD\s*=\s*([\d\.]+)', content)
            axis_match = re.search(r'AXIS\s*,?\s*=\s*([\d\.]+)', content)
            daxis_match = re.search(r'AXISERR\s*,?\s*=\s*([\d\.]+)', content)
            ecc_match = re.search(r'ECC\s*=\s*([\d\.]+)', content)
            decc_match = re.search(r'ECCERR\s*=\s*([\d\.]+)', content)
            inc_match = re.search(r'INCL\s*=\s*([\d\.]+)', content)
            dinc_match = re.search(r'INCLERR\s*=\s*([\d\.]+)', content)
            Om_match = re.search(r'NODE\s*=\s*([\d\.]+)', content)
            dOm_match = re.search(r'NODEERR\s*=\s*([\d\.]+)', content)
            om_match = re.search(r'ARGUP\s*=\s*([\d\.]+)', content)
            dom_match = re.search(r'ARGUPERR\s*=\s*([\d\.]+)', content)

            if jd_match:
                jd = float(jd_match.group(1))
            if axis_match:
                a = float(axis_match.group(1))* au.value * 1.e-3
            if daxis_match:
                da = float(daxis_match.group(1))* au.value * 1.e-3
            if ecc_match:
                e = float(ecc_match.group(1))
            if decc_match:
                de = float(decc_match.group(1))
            if inc_match:
                i = np.radians(float(inc_match.group(1)))
            if dinc_match:
                di = np.radians(float(dinc_match.group(1)))
            if Om_match:
                Omg = np.radians(float(Om_match.group(1)))
            if dOm_match:
                dOmg = np.radians(float(dOm_match.group(1)))
            if om_match:
                omg = np.radians(float(om_match.group(1)))
            if dom_match:
                domg = np.radians(float(dom_match.group(1)))


    pos = np.array([
        a,
        da,
        e,
        de,
        i,
        di,
        Omg,
        dOmg,
        omg,
        domg
    ])


    return jd, pos






    TB = time_back/(365.25*24.*3600.)
    TF = time_forw/(365.25*24.*3600.)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    fig.suptitle(fig_title, fontsize=16)

    axes[0, 0].plot(TB, back_q[0]/(au.value*1.e-3), color=clr, linewidth=2)
    axes[0, 0].plot(TF, forw_q[0]/(au.value*1.e-3), color=clr, linewidth=2)

    for j in range(1, n_clones + 1):
        axes[0, 0].plot(TB, back_q[j]/(au.value*1.e-3), color=clr, alpha=0.2)
        axes[0, 0].plot(TF, forw_q[j]/(au.value*1.e-3), color=clr, alpha=0.2)

    #axes[1, 0].set_title('Plot 3')
    axes[0, 0].set_xlabel('time (years)')
    axes[0, 0].set_ylabel('Perihelion distance (AU)')
    #axes[1, 0].legend()

    # axes[0, 1].plot(TB, data[0], color=clr, linewidth=2, label='Curva principal')

    axes[0, 1].plot(TB, back_D_show[0], color=clr, linewidth=2)
    axes[0, 1].plot(TF, forw_D_show[0], color=clr, linewidth=2)

    for j in range(1, n_clones + 1):
        axes[0, 1].plot(TB, back_D_show[j], color=clr, alpha=0.2)
        axes[0, 1].plot(TF, forw_D_show[j], color=clr, alpha=0.2)

    axes[0, 1].set_title('Meteor Shower')
    axes[0, 1].set_xlabel('time (years)')
    axes[0, 1].set_ylabel('D criterion')
    #axes[1, 1].legend()


    axes[1, 0].plot(TB, back_D_HY[0], color=clr, linewidth=2)
    axes[1, 0].plot(TF, forw_D_HY[0], color=clr, linewidth=2)

    for j in range(1, n_clones + 1):
        axes[1, 0].plot(TB, back_D_HY[j], color=clr, alpha=0.2)
        axes[1, 0].plot(TF, forw_D_HY[j], color=clr, alpha=0.2)

    axes[1, 0].set_title("2010 HY22")
    axes[1, 0].set_xlabel('time (years)')
    axes[1, 0].set_ylabel('D criterion')
    #axes[1, 1].legend()


    axes[1, 1].plot(TB, back_D_OA[0], color=clr, linewidth=2)
    axes[1, 1].plot(TF, forw_D_OA[0], color=clr, linewidth=2)

    for j in range(1, n_clones + 1):
        axes[1, 1].plot(TB, back_D_OA[j], color=clr, alpha=0.2)
        axes[1, 1].plot(TF, forw_D_OA[j], color=clr, alpha=0.2)

    axes[1, 1].set_title("2015 OA22")
    axes[1, 1].set_xlabel('time (years)')
    axes[1, 1].set_ylabel('D criterion')
    #axes[1, 1].legend()
    #axes[0, 0].plt.yscale('log')
    axes[1, 0].set_yscale('log')
    axes[0, 1].set_yscale('log')
    axes[1, 1].set_yscale('log')
    plt.tight_layout()

    plt.savefig("log_" + fig_name)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    fig.suptitle(fig_title, fontsize=16)

    axes[0, 0].plot(TB, back_q[0]/(au.value*1.e-3), color=clr, linewidth=2)
    axes[0, 0].plot(TF, forw_q[0]/(au.value*1.e-3), color=clr, linewidth=2)

    for j in range(1, n_clones + 1):
        axes[0, 0].plot(TB, back_q[j]/(au.value*1.e-3), color=clr, alpha=0.2)
        axes[0, 0].plot(TF, forw_q[j]/(au.value*1.e-3), color=clr, alpha=0.2)

    #axes[1, 0].set_title('Plot 3')
    axes[0, 0].set_xlabel('time (years)')
    axes[0, 0].set_ylabel('Perihelion distance (AU)')
    #axes[1, 0].legend()

    # axes[0, 1].plot(TB, data[0], color=clr, linewidth=2, label='Curva principal')

    axes[0, 1].plot(TB, back_D_show[0], color=clr, linewidth=2)
    axes[0, 1].plot(TF, forw_D_show[0], color=clr, linewidth=2)

    for j in range(1, n_clones + 1):
        axes[0, 1].plot(TB, back_D_show[j], color=clr, alpha=0.2)
        axes[0, 1].plot(TF, forw_D_show[j], color=clr, alpha=0.2)

    axes[0, 1].set_title('Meteor Shower')
    axes[0, 1].set_xlabel('time (years)')
    axes[0, 1].set_ylabel('Drummond criterion')
    #axes[1, 1].legend()


    axes[1, 0].plot(TB, back_D_HY[0], color=clr, linewidth=2)
    axes[1, 0].plot(TF, forw_D_HY[0], color=clr, linewidth=2)

    for j in range(1, n_clones + 1):
        axes[1, 0].plot(TB, back_D_HY[j], color=clr, alpha=0.2)
        axes[1, 0].plot(TF, forw_D_HY[j], color=clr, alpha=0.2)

    axes[1, 0].set_title("2010 HY22")
    axes[1, 0].set_xlabel('time (years)')
    axes[1, 0].set_ylabel('Drummond criterion')
    #axes[1, 1].legend()


    axes[1, 1].plot(TB, back_D_OA[0], color=clr, linewidth=2)
    axes[1, 1].plot(TF, forw_D_OA[0], color=clr, linewidth=2)

    for j in range(1, n_clones + 1):
        axes[1, 1].plot(TB, back_D_OA[j], color=clr, alpha=0.2)
        axes[1, 1].plot(TF, forw_D_OA[j], color=clr, alpha=0.2)

    axes[1, 1].set_title("2015 OA22")
    axes[1, 1].set_xlabel('time (years)')
    axes[1, 1].set_ylabel('Drummond criterion')
    #axes[1, 1].legend()
    #axes[0, 0].plt.yscale('log')

    plt.tight_layout()

    plt.savefig(fig_name)



# READ FILE
#mpc = mpcorbfile.mpcorb_file('MPCORB.DAT')
mpc = mpcorbfile.mpcorb_file('MPCORB-11-MAY.DAT')

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

    data[1,i] = Drummond_criterion(e1,e2,q1,q2,i1,i2,O1,O2,o1,o2)
    data[2,i] = D_Nesvorny(a1,a2,e1,e2,i1,i2,O1,O2,vp1,vp2)



a1 = mpc.bodies[i]['a']
e1 = mpc.bodies[i]['e']
q1 = a1*(1.-e1)
i1 = np.radians(mpc.bodies[i]['i'])%(2.*np.pi)
O1 = np.radians(mpc.bodies[i]['Node'])%(2.*np.pi)
o1 = np.radians(mpc.bodies[i]['Peri'])%(2.*np.pi)
vp1 = (O1 + o1)%(2.*np.pi)




print(Drummond_criterion(e1,e2,q1,q2,i1,i2,O1,O2,o1,o2))






sort_D = np.argsort(data[1, :])
sort_N = np.argsort(data[2, :])
Dru = data[:, sort_D]
Ner = data[:, sort_N]


redo_N = 500
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

    print(f'{id_ast[int(rd_Dru[0,i])],int(rd_Dru[0,i])}', mpc.bodies[int(rd_Dru[0,i])]['packed_designation'], mpc.bodies[int(rd_Dru[0,i])]['designation'])


    a1 = elements["a"].value[0]
    e1 = elements["e"].value[0]
    q1 = a1*(1.-e1)
    i1 = np.radians(elements["incl"].value[0])%(2.*np.pi)
    O1 = np.radians(elements["Omega"].value[0])%(2.*np.pi)
    o1 = np.radians(elements["w"].value[0])%(2.*np.pi)
    vp1 = (O1 + o1)%(2.*np.pi)
    rd_Dru[2,i] = Drummond_criterion(e1,e2,q1,q2,i1,i2,O1,O2,o1,o2)



obj = Horizons(id=f'K18P10L', id_type='smallbody', location='@sun', epochs=jd_date)
elements = obj.elements()
a1 = elements["a"].value[0]
e1 = elements["e"].value[0]
q1 = a1*(1.-e1)
i1 = np.radians(elements["incl"].value[0])%(2.*np.pi)
O1 = np.radians(elements["Omega"].value[0])%(2.*np.pi)
o1 = np.radians(elements["w"].value[0])%(2.*np.pi)
vp1 = (O1 + o1)%(2.*np.pi)
print(Drummond_criterion(e1,e2,q1,q2,i1,i2,O1,O2,o1,o2))

sort_D = np.argsort(rd_Dru[2, :])

Dru_N = rd_Dru[:, sort_D]



for i in range(10):
    print(i, mpc.bodies[int(Dru[0,i])]['Name'], Dru[1,i], mpc.bodies[int(Dru_N[0,i])]['Name'], Dru_N[1,i],Dru_N[2,i])



