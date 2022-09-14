import scipy.special as sp
import numpy as np


def k_m(k_cal,K,cs2m=1,cs2p=1):
    return np.sqrt(cs2m*k_cal-K*(9*(cs2m**2)/2+16*(cs2m)+(97/6)+(k_cal*(9*(cs2m**3)+36*(cs2m**2)+59*(cs2m))+96*K)/(6*K*(1-cs2m)-4*k_cal*cs2m)), dtype=np.complex128)
def k_p(k_cal,K,cs2m=1,cs2p=1):
    to_ret=np.sqrt(cs2p*k_cal+K*(1/3-3*cs2p), dtype=np.complex128)
    return to_ret

#Defining the primordial power spectrum
def P(k, nt, K,cs2m=1,cs2p=1):
    m=(3*cs2m+1)/2;
    if K == +1:
        k_cal = k*(k+2)
    elif K == 0 or K == -1:
        k_cal = k**2
    k_mi = np.array(k_m(k_cal,K,cs2m,cs2p))
    k_pl = np.array(k_p(k_cal,K,cs2m,cs2p))

    # h2_1=sp.hankel2(1/2, m*k_pl*nt)
    # h2_pm=sp.hankel2(1/2-(1/m), k_mi*nt)
    #
    # h1_1=sp.hankel1(1/2, m*k_pl*nt)
    # h1_3=sp.hankel1(3/2, m*k_pl*nt)
    #
    # h2_3=sp.hankel2(3/2, m*k_pl*nt)
    # h2_nm=sp.hankel2(-1/2-(1/m), k_mi*nt)


    h2_1_l=np.log(sp.hankel2e(1/2, m*k_pl*nt))-1j*(m*k_pl*nt)
    h2_pm_l=np.log(sp.hankel2e(1/2-(1/m), k_mi*nt))-1j*(k_mi*nt)

    h1_1_l=np.log(sp.hankel1e(1/2, m*k_pl*nt))+1j*(m*k_pl*nt)
    h1_3_l=np.log(sp.hankel1e(3/2, m*k_pl*nt))+1j*(m*k_pl*nt)

    h2_3_l=np.log(sp.hankel2e(3/2, m*k_pl*nt))-1j*(m*k_pl*nt)
    h2_nm_l=np.log(sp.hankel2e(-1/2-(1/m), k_mi*nt))-1j*(k_mi*nt)


    c_k_1_log=np.log(k_pl)+h2_1_l+h2_pm_l
    c_k_2_log=np.log(k_mi)+h2_3_l+h2_nm_l

    d_k_1_log=np.log(k_pl)+h1_1_l+h2_pm_l
    d_k_2_log=np.log(k_mi)+h1_3_l+h2_nm_l


    # C_K=1j*(np.pi*np.sqrt(m)*nt/4)*(k_pl*h2_1*h2_pm+k_mi*h2_3*h2_nm)/np.power(k_pl,3/2)
    # D_K=-1j*(np.pi*np.sqrt(m)*nt/4)*(k_pl*h1_1*h2_pm+k_mi*h1_3*h2_nm)/np.power(k_pl,3/2)

    C_K_log=np.log(1j*(np.pi*np.sqrt(m)*nt/4))+ sp.logsumexp([c_k_1_log,c_k_2_log],axis=0)-(3/2)*np.log(k_pl)
    minD_K_log=np.log(1j*(np.pi*np.sqrt(m)*nt/4))+ sp.logsumexp([d_k_1_log,d_k_2_log],axis=0)-(3/2)*np.log(k_pl)

    diff_log=sp.logsumexp([C_K_log,minD_K_log],axis=0)
    if(np.isinf(np.power(cs2p,(3/2))*((k**3)*(np.exp(2*np.real(diff_log))))).any()):
        print(K)
        print(nt/(np.pi/(2*m)))
        # print(np.real(diff_log).max())
        print(np.real(diff_log).max(),k_pl[np.real(diff_log).argmax()],k_mi[np.real(diff_log).argmax()])
        print(np.real(diff_log).min(),k_pl[np.real(diff_log).argmin()],k_mi[np.real(diff_log).argmin()])
        print(k[np.real(diff_log).argmax()],k[np.real(diff_log).argmin()])

    # print('We have inaginary: ',np.isreal(k_mi).sum())
    val_to_ret=np.zeros(k.shape)+1e-50
    val_to_ret[np.isreal(k_mi)]=(np.power(cs2p,(3/2))*((k**3)*(np.exp(np.float128(2*np.real(diff_log))))))[np.isreal(k_mi)] *2/(np.sqrt(cs2m/cs2p) + np.sqrt(cs2p/cs2m))
    # val_to_ret=(np.power(cs2,(3/2))*((k**3)*(np.exp(np.float128(2*np.real(diff_log))))))
    return val_to_ret


def pps(universe,ks,kcom,nt,cs2m,cs2p,baseline=False):
    if(baseline):
        return universe.As*((ks/universe.kp)**(universe.ns-1))
    return universe.As*((ks/universe.kp)**(universe.ns-1))*P(kcom,nt,universe.K,cs2m,cs2p)
