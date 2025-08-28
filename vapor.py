import numpy as np 
from scipy.optimize import approx_fprime
class EqState:
    def __init__(self,pressao:float,temperatura:float,pressao_critica:float,temperatura_critica:float,fator_acentrico:float)->None:
        self.r = 83.14  # bar.cm^3/mol.K
        self.__tc = temperatura_critica
        self.__pc = pressao_critica
        self.p = pressao
        self.t = temperatura
        self.__w = fator_acentrico
    def __beta_pr(self)->float: # parâmetro beta de peng-robinson 
        omega = 0.07780
        pr = self.p/self.__pc
        tr = self.t/self.__tc
        return omega*(pr/tr)
    def __q_pr(self)->float: # parâmetro q de peng-robinson 
        omega = 0.07780
        psi = 0.45724
        tr = self.t/self.__tc
        alpha = (1+(0.37464+1.54226*self.__w-0.26992*(self.__w**2))*(1-(tr**0.5)))**2
        return (psi*alpha)/(omega*tr)
    def __beta_vdw(self)->float: # parâmetro beta de vw 
        omega = 1/8
        pr = self.p/self.__pc
        tr = self.t/self.__tc
        return omega*(pr/tr)
    def __q_vdw(self)->float: # parâmetro q de vw
        omega = 1/8
        psi = 27/64
        alpha = 1
        tr = self.t/self.__tc
        return (psi*alpha)/(omega*tr)
    def __beta_srk(self)->float: # parâmetro beta de srk
        omega = 0.08664
        pr = self.p/self.__pc
        tr = self.t/self.__tc
        return omega*(pr/tr)
    def __q_srk(self)->float: # parâmetro q de srk
        omega = 0.08664
        psi = 0.42748
        tr = self.t/self.__tc
        alpha = (1+(0.480+1.574*self.__w-0.176*(self.__w**2))*(1-(tr**0.5)))**2
        return (psi*alpha)/(omega*tr)
    def __beta_rk(self)->float: # parâmetro beta de rk
        omega = 0.08664
        pr = self.p/self.__pc
        tr = self.t/self.__tc
        return omega*(pr/tr)
    def __q_rk(self)->float: # parâmetro q de rk
        omega = 0.08664
        psi = 0.42748
        tr = self.t/self.__tc
        alpha = tr**(-0.5)
        return (psi*alpha)/(omega*tr)
    def __sigma_pr(self)->float: # parâmetro sigma de pr 
        return 1+(2**(0.5))
    def __epsilon_pr(self)->float: # parâmetro epsilon de pr 
        return 1-(2**(0.5))
    def __sigma_vdw(self)->float: # parâmetro sigma de vw
        return 0
    def __epsilon_vdw(self)->float: # parâmetro epsilon de vw 
        return 0
    def __sigma_rk(self)->float: # parâmetro sigma de rk
        return 1
    def __epsilon_rk(self)->float: # parâmetro epsilon de rk
        return 0
    def __sigma_srk(self)->float: # parâmetro sigma de srk
        return 1
    def __epsilon_srk(self)->float: # parâmetro epsilon de srk
        return 0
    #------------------------------------parâmetros das equações de estado---------------------------------------------
    #-----------------------------------------------pr-----------------------------------------------------------------
    def betapr(self)->float: # parâmetro beta pr
        beta_pr = self.__beta_pr()
        return beta_pr
    def qpr(self)->float: # parâmetro q pr 
        q_pr = self.__q_pr()
        return q_pr
    def sigmapr(self)->float: # parâmetro sigma pr 
        sigma_pr = self.__sigma_pr()
        return sigma_pr
    def espsilonpr(self)->float: # parâmetro epsilon pr 
        epsilon_pr = self.__epsilon_pr()
        return epsilon_pr
    #---------------------------------------------vw-------------------------------------------------------------------
    def betavw(self)->float: # parâmetro beta vw
        beta_vw = self.__beta_vdw()
        return beta_vw
    def qvw(self)->float: # parâmetro q vw
        q_vw = self.__q_vdw()
        return q_vw
    def sigmavw(self)->float: # parâmetro sigma vw
        sigma_vw = self.__sigma_vdw()
        return sigma_vw
    def espsilonvw(self)->float: # parâmetro epsilon vw
        epsilon_vw = self.__epsilon_vdw()
        return epsilon_vw
    #--------------------------------------------rk--------------------------------------------------------------------
    def betark(self)->float: # parâmetro beta rk
        beta_rk = self.__beta_rk()
        return beta_rk
    def qrk(self)->float: # parâmetro q rk
        q_rk = self.__q_rk()
        return q_rk
    def sigmark(self)->float: # parâmetro sigma rk
        sigma_rk = self.__sigma_rk()
        return sigma_rk
    def espsilonrk(self)->float: # parâmetro epsilon rk
        epsilon_rk = self.__epsilon_rk()
        return epsilon_rk 
    #---------------------------------------------srk------------------------------------------------------------------
    def betasrk(self)->float: # parâmetro beta srk
        beta_srk = self.__beta_srk()
        return beta_srk
    def qsrk(self)->float: # parâmetro q srk
        q_srk = self.__q_srk()
        return q_srk
    def sigmasrk(self)->float: # parâmetro sigma srk
        sigma_srk = self.__sigma_srk()
        return sigma_srk
    def espsilonsrk(self)->float: # parâmetro epsilon srk
        epsilon_srk = self.__epsilon_srk()
        return epsilon_srk 
    #------------------------------------RETORNO DE TC E PC e FATOR ACÊNTRICO ----------------------------------------------------------
    def temp_tc(self)->float: # retorno para a temperatura crítica
        return self.__tc
    def press_pc(self)->float: # retorno para a pressão crítica
        return self.__pc
    def w_c(self)->float: # retorno para o fator acêntrico 
        return self.__w
    
class vapor(EqState):
    def __init__(self, pressao:float, temperatura:float, pressao_critica:float, temperatura_critica:float, fator_acentrico:float,nmx:float,tol:float):
        self.nmx = nmx
        self.tol = tol
        super().__init__(pressao, temperatura, pressao_critica, temperatura_critica, fator_acentrico)
    #------------------------------------------------Z,V,PARA PR--------------------------------------------------------------------------------------
    def z_pr(self)->float: # retorna o coeficiente de compressibilidade para pr 
        erro = 1
        it = 0
        i = 1 
        while abs(erro)>self.tol and it<self.nmx: 
            i = 1 + (self.betapr()) - (self.qpr()*self.betapr())*((i-self.betapr())/((i+(self.espsilonpr()*self.betapr()))*(i+(self.sigmapr()*self.betapr()))))
            erro = erro - i
            it = it + 1
            #print(it,i)
        return i 
    def v_pr(self)->float: # retorna o volume molar em cm^3/mol para peng-robinson 
        volume = (self.z_pr()*self.r*self.t)/(self.p)
        return volume
    #-----------------------------------------------Z,V PARA VDW--------------------------------------------------------------------------------------
    def z_vdw(self)->float: # retorna o coeficiente de compressibilidade para vdw 
        erro = 1
        it = 0
        i = 1 
        while abs(erro)>self.tol and it<self.nmx: 
            i = 1 + (self.betavw()) - (self.qvw()*self.betavw())*((i-self.betavw())/((i+(self.espsilonvw()*self.betavw()))*(i+(self.sigmavw()*self.betavw()))))
            erro = erro - i
            it = it + 1
            #print(it,i)
        return i 
    def v_vdw(self)->float: # retorna o volume molar em cm^3/mol para vdw
        volume = (self.z_vdw()*self.r*self.t)/(self.p)
        return volume
    #---------------------------------------------Z,V PARA RK------------------------------------------------------------------------------------------
    def z_rk(self)->float: # retorna o coeficiente de compressibilidade para rk
        erro = 1
        it = 0
        i = 1 
        while abs(erro)>self.tol and it<self.nmx: 
            i = 1 + (self.betark()) - (self.qrk()*self.betark())*((i-self.betark())/((i+(self.espsilonrk()*self.betark()))*(i+(self.sigmark()*self.betark()))))
            erro = erro - i
            it = it + 1
            #print(it,i)
        return i 
    def v_rk(self)->float: # retorna o volume molar em cm^3/mol para rk
        volume = (self.z_rk()*self.r*self.t)/(self.p)
        return volume
    #--------------------------------------------Z,V PARA SRK-------------------------------------------------------------------------------------------
    def z_srk(self)->float: # retorna o coeficiente de compressibilidade para srk
        erro = 1
        it = 0
        i = 1 
        while abs(erro)>self.tol and it<self.nmx: 
            i = 1 + (self.betasrk()) - (self.qsrk()*self.betasrk())*((i-self.betasrk())/((i+(self.espsilonsrk()*self.betasrk()))*(i+(self.sigmasrk()*self.betasrk()))))
            erro = erro - i
            it = it + 1
            #print(it,i)
        return i 
    def v_srk(self)->float: # retorna o volume molar em cm^3/mol para srk
        volume = (self.z_srk()*self.r*self.t)/(self.p)
        return volume
    
class Auxiliar(vapor):
    def __init__(self, pressao, temperatura, pressao_critica, temperatura_critica, fator_acentrico, nmx, tol):
        super().__init__(pressao, temperatura, pressao_critica, temperatura_critica, fator_acentrico, nmx, tol)

    def dalpha_tr_pr_tr(self): # derivade de alpha de pr em função da temperatura reduzida 
        omega = 0.07780
        psi = 0.45724
        tr = self.t/self.temp_tc()
        # h = 0.000001
        # p1 = tr - 2*h
        # p2 = tr - h
        # p3 = tr + h
        # p4 = tr + 2*h
        # fp1 = (((1+(0.37464+1.54226*self.w_c()-0.26992*(self.w_c()**2))*(1-(p1**0.5)))**2)*psi)/(omega*p1)
        # fp2 = (((1+(0.37464+1.54226*self.w_c()-0.26992*(self.w_c()**2))*(1-(p2**0.5)))**2)*psi)/(omega*p2)
        # fp3 = (((1+(0.37464+1.54226*self.w_c()-0.26992*(self.w_c()**2))*(1-(p3**0.5)))**2)*psi)/(omega*p3)
        # fp4 = (((1+(0.37464+1.54226*self.w_c()-0.26992*(self.w_c()**2))*(1-(p4**0.5)))**2)*psi)/(omega*p4)
        # d = (fp1 - 8*fp2 + 8*fp3 - fp4)/(12*h)
        def dpr(tr):
            return (((1+(0.37464+1.54226*self.w_c()-0.26992*(self.w_c()**2))*(1-(tr**0.5)))**2)*psi)/(omega*tr)
        t0_pr = np.array([tr])
        dx_pr = 1e-6
        derivadapr = approx_fprime(t0_pr,dpr,dx_pr)
        return derivadapr[0]
        
    def dalpha_tr_srk_tr(self): # derivada de alpha de srk em função da temperatura reduzida 
        omega = 0.08664
        psi = 0.42748
        tr = self.t/self.temp_tc()
        h = 0.000001
        p1 = tr - 2*h
        p2 = tr - h
        p3 = tr + h
        p4 = tr + 2*h
        fp1 = (((1+(0.480+1.574*self.w_c()-0.176*(self.w_c()**2))*(1-(p1**0.5)))**2)*psi)/(omega*p1)
        fp2 = (((1+(0.480+1.574*self.w_c()-0.176*(self.w_c()**2))*(1-(p2**0.5)))**2)*psi)/(omega*p2)
        fp3 = (((1+(0.480+1.574*self.w_c()-0.176*(self.w_c()**2))*(1-(p3**0.5)))**2)*psi)/(omega*p3)
        fp4 = (((1+(0.480+1.574*self.w_c()-0.176*(self.w_c()**2))*(1-(p4**0.5)))**2)*psi)/(omega*p4)
        d = (fp1 - 8*fp2 + 8*fp3 - fp4)/(12*h)
        return d
    
    def dalpha_tr_rk_tr(self): # derivade de alpha de rk em função da temperatura reduzida
        omega = 0.08664
        psi = 0.42748 
        tr = self.t/self.temp_tc()
        h = 0.000001
        p1 = tr - 2*h
        p2 = tr - h
        p3 = tr + h
        p4 = tr + 2*h
        fp1 = ((p1**(-0.5))*psi)/(omega*p1)
        fp2 = ((p2**(-0.5))*psi)/(omega*p2)
        fp3 = ((p3**(-0.5))*psi)/(omega*p3)
        fp4 = ((p4**(-0.5))*psi)/(omega*p4)
        d = (fp1 - 8*fp2 + 8*fp3 - fp4)/(12*h)
        return d
    
    def dalpha_tr_vdw_tr(self): # derivade de alpha de vw em função da temperatura reduzida
        omega = 1/8
        psi = 27/64
        tr = self.t/self.temp_tc()
        h = 0.000001
        p1 = tr - 2*h
        p2 = tr - h
        p3 = tr + h
        p4 = tr + 2*h
        fp1 = (psi)/(omega*p1)
        fp2 = (psi)/(omega*p2)
        fp3 = (psi)/(omega*p3)
        fp4 = (psi)/(omega*p4)
        d = (fp1 - 8*fp2 + 8*fp3 - fp4)/(12*h)
        return d
    
    def i_pr(self): # retorna o parâmetro i de pr 
        ipr = (1/(self.sigmapr() - self.espsilonpr()))*(np.log((self.z_pr()+(self.sigmapr()*self.betapr()))/(self.z_pr()+(self.espsilonpr()*self.betapr()))))
        return ipr 
    def i_srk(self): # retorna o parâmetro i de srk
        isrk = (1/(self.sigmasrk() - self.espsilonsrk()))*(np.log((self.z_srk()+(self.sigmasrk()*self.betasrk()))/(self.z_srk()+(self.espsilonsrk()*self.betasrk()))))
        return isrk
    def i_rk(self): # retorna o parâmetro i de rk
        irk = (1/(self.sigmark() - self.espsilonrk()))*(np.log((self.z_rk()+(self.sigmark()*self.betark()))/(self.z_rk()+(self.espsilonrk()*self.betark()))))
        return irk
    def i_vdw(self): # retorna o parâmetro i de vdw
        ivdw = (self.betavw())/(self.z_vdw()+(self.espsilonvw()*self.betavw()))
        return ivdw

class Prop(Auxiliar):
    def __init__(self, pressao, temperatura, pressao_critica, temperatura_critica, fator_acentrico, nmx, tol):
        super().__init__(pressao, temperatura, pressao_critica, temperatura_critica, fator_acentrico, nmx, tol)
    #------------------------------------------entalpia_pr----------------------------------------------------------------------------------------------------
    def cp_gasid_formacao_298k(self,a,b,c,d): #h(298.15k)=0
        integral = (a*(self.t-298.15)) + ((b/2)*(((self.t**2)-(298.15**2)))) + ((c/3)*(((self.t**3)-(298.15**3)))) + ((d/4)*(((self.t**4)-(298.15**4))))
        return integral
    def entalpia_gas_ideal_formacao_models(self,a,b,c,d):
        h_id = ((self.cp_gasid_formacao_298k(a,b,c,d)))
        return h_id
    def entalpia_residual_pr(self):
        r = 8.314 # kJ/kmol.K
        tr = self.t/self.temp_tc()
        hr_pr = (self.z_pr() - 1 + (tr*self.dalpha_tr_pr_tr()*self.i_pr()))*r*self.t
        return hr_pr
    def entalpia_absolute_pr(self,a,b,c,d):
        h_absolute  = self.entalpia_gas_ideal_formacao_models(a,b,c,d)+self.entalpia_residual_pr()
        return h_absolute
    #------------------------------------------entropia_pr---------------------------------------------------------------------------------------------------
    def entropia_residual_pr(self):
        r = 8.314 # kJ/kmol.K
        tr = self.t/self.temp_tc()
        s_res_pr = ((np.log(self.z_pr() - self.betapr())) + ((self.qpr()+(tr*self.dalpha_tr_pr_tr()))*(self.i_pr())))*r
        return s_res_pr
    def entropia_gas_ideal_formacao(self,a,b,c,d): #s(298.15k,1bar)=0
        r = 8.314 # kJ/kmol.K
        sid = a*(np.log((self.t)/(298.15))) + b*(self.t - 298.15) + (c/2)*(((self.t)**2) - (298.15**2)) + (d/3)*(((self.t)**3) - (298.15**3)) - r*np.log(self.p)
        return sid
    def entropia_absolute_pr(self,a,b,c,d):
        s_absolute = self.entropia_gas_ideal_formacao(a,b,c,d)+self.entropia_residual_pr()
        return s_absolute
    #------------------------------------------gibbs_pr------------------------------------------------------------------------------------------------------
    def gibbs_absolute_pr(self,a,b,c,d):
        gibbs_absolute = self.entalpia_absolute_pr(a,b,c,d) - ((self.t)*self.entropia_absolute_pr(a,b,c,d))
        return gibbs_absolute
    #------------------------------------------interna_pr-------------------------------------------------------------------------------------------------
    def interna_absolute_pr(self,a,b,c,d):
        r = 8.314 # kJ/kmol.K
        interna_absolute = self.entalpia_absolute_pr(a,b,c,d) - (self.z_pr()*r*self.t)
        return interna_absolute
    #-----------------------------------------helmholtz_pr------------------------------------------------------------------------------------------------
    def helmholtz_absolute_pr(self,a,b,c,d):
        helmholtz_absolute = self.interna_absolute_pr(a,b,c,d) - (self.t*self.entropia_absolute_pr(a,b,c,d))
        return helmholtz_absolute
    #-----------------------------------------------entalpia_srk------------------------------------------------------------------------------------------
    def entalpia_residual_srk(self):
        r = 8.314 # kJ/kmol.K
        tr = self.t/self.temp_tc()
        hr_srk = (self.z_srk() - 1 + (tr*self.dalpha_tr_srk_tr()*self.i_srk()))*r*self.t
        return hr_srk
    def entalpia_absolute_srk(self,a,b,c,d):
        h_absolute_srk  = self.entalpia_gas_ideal_formacao_models(a,b,c,d)+self.entalpia_residual_srk()
        return h_absolute_srk
    #----------------------------------------------entropia_srk-------------------------------------------------------------------------------------------
    def entropia_residual_srk(self):
        r = 8.314 # kJ/kmol.K
        tr = self.t/self.temp_tc()
        s_res_srk = ((np.log(self.z_srk() - self.betasrk())) + ((self.qsrk()+(tr*self.dalpha_tr_srk_tr()))*(self.i_srk())))*r
        return s_res_srk
    def entropia_absolute_srk(self,a,b,c,d):
        s_absolute_srk = self.entropia_gas_ideal_formacao(a,b,c,d)+self.entropia_residual_srk()
        return s_absolute_srk
    #--------------------------------------------gibbs_srk------------------------------------------------------------------------------------------------
    def gibbs_absolute_srk(self,a,b,c,d):
        gibbs_absolutesrk = self.entalpia_absolute_srk(a,b,c,d) - ((self.t)*self.entropia_absolute_srk(a,b,c,d))
        return gibbs_absolutesrk
    #--------------------------------------------interna_srk----------------------------------------------------------------------------------------------
    def interna_absolute_srk(self,a,b,c,d):
        r = 8.314 # kJ/kmol.K
        interna_absolutesrk = self.entalpia_absolute_srk(a,b,c,d) - (self.z_srk()*r*self.t)
        return interna_absolutesrk
    #-------------------------------------------helmoltz_srk----------------------------------------------------------------------------------------------
    def helmholtz_absolute_srk(self,a,b,c,d):
        helmholtz_absolutesrk = self.interna_absolute_srk(a,b,c,d) - (self.t*self.entropia_absolute_srk(a,b,c,d))
        return helmholtz_absolutesrk
    #-------------------------------------------entalpia_rk-----------------------------------------------------------------------------------------------
    def entalpia_residual_rk(self):
        r = 8.314 # kJ/kmol.K
        tr = self.t/self.temp_tc()
        hr_rk = (self.z_rk() - 1 + (tr*self.dalpha_tr_rk_tr()*self.i_rk()))*r*self.t
        return hr_rk
    def entalpia_absolute_rk(self,a,b,c,d):
        h_absolute_rk  = self.entalpia_gas_ideal_formacao_models(a,b,c,d)+self.entalpia_residual_rk()
        return h_absolute_rk
    #------------------------------------------entropia_rk------------------------------------------------------------------------------------------------
    def entropia_residual_rk(self):
        r = 8.314 # kJ/kmol.K
        tr = self.t/self.temp_tc()
        s_res_rk = ((np.log(self.z_rk() - self.betark())) + ((self.qrk()+(tr*self.dalpha_tr_rk_tr()))*(self.i_rk())))*r
        return s_res_rk
    def entropia_absolute_rk(self,a,b,c,d):
        s_absolute_rk = self.entropia_gas_ideal_formacao(a,b,c,d)+self.entropia_residual_rk()
        return s_absolute_rk
    #-----------------------------------------gibbs_rk----------------------------------------------------------------------------------------------------
    def gibbs_absolute_rk(self,a,b,c,d):
        gibbs_absoluterk = self.entalpia_absolute_rk(a,b,c,d) - ((self.t)*self.entropia_absolute_rk(a,b,c,d))
        return gibbs_absoluterk
    #---------------------------------------interna_rk------------------------------------------------------------------------------------------------------
    def interna_absolute_rk(self,a,b,c,d):
        r = 8.314 # kJ/kmol.K
        interna_absoluterk = self.entalpia_absolute_rk(a,b,c,d) - (self.z_rk()*r*self.t)
        return interna_absoluterk
    #---------------------------------------helmoltz_rk-----------------------------------------------------------------------------------------------------
    def helmholtz_absolute_rk(self,a,b,c,d):
        helmholtz_absoluterk = self.interna_absolute_rk(a,b,c,d) - (self.t*self.entropia_absolute_rk(a,b,c,d))
        return helmholtz_absoluterk
    #--------------------------------------entalpiavw--------------------------------------------------------------------------------------------------------
    def entalpia_residual_vw(self):
        r = 8.314 # kJ/kmol.K
        tr = self.t/self.temp_tc()
        hr_vw = (self.z_vdw() - 1 + (tr*self.dalpha_tr_vdw_tr()*self.i_vdw()))*r*self.t
        return hr_vw
    def entalpia_absolute_vw(self,a,b,c,d):
        h_absolute_vw  = self.entalpia_gas_ideal_formacao_models(a,b,c,d)+self.entalpia_residual_vw()
        return h_absolute_vw
    #------------------------------------entropiavw---------------------------------------------------------------------------------------------------------
    def entropia_residual_vw(self):
        r = 8.314 # kJ/kmol.K
        tr = self.t/self.temp_tc()
        s_res_vw = ((np.log(self.z_vdw() - self.betavw())) + ((self.qvw()+(tr*self.dalpha_tr_vdw_tr()))*(self.i_vdw())))*r
        return s_res_vw
    def entropia_absolute_vw(self,a,b,c,d):
        s_absolute_vw = self.entropia_gas_ideal_formacao(a,b,c,d)+self.entropia_residual_vw()
        return s_absolute_vw
    #----------------------------------gibbsvw--------------------------------------------------------------------------------------------------------------
    def gibbs_absolute_vw(self,a,b,c,d):
        gibbs_absolutevw = self.entalpia_absolute_vw(a,b,c,d) - ((self.t)*self.entropia_absolute_vw(a,b,c,d))
        return gibbs_absolutevw
    #---------------------------------------internavw-------------------------------------------------------------------------------------------------------
    def interna_absolute_vw(self,a,b,c,d):
        r = 8.314 # kJ/kmol.K
        interna_absolutevw = self.entalpia_absolute_vw(a,b,c,d) - (self.z_vdw()*r*self.t)
        return interna_absolutevw
    #---------------------------------------helmoltzvw-------------------------------------------------------------------------------------------------------
    def helmholtz_absolute_vw(self,a,b,c,d):
        helmholtz_absolutevw = self.interna_absolute_vw(a,b,c,d) - (self.t*self.entropia_absolute_vw(a,b,c,d))
        return helmholtz_absolutevw
    #--------------------------------------------------fugacidades/coeficientes_de_fugacidade------------------------------------------------------------------
    #------------------------------------------------------------pr--------------------------------------------------------------------------------------------
    def coeficiente_fug_pr(self):
        lnc_pr = self.z_pr() - 1 - (np.log(self.z_pr() - self.betapr())) - (self.qpr()*self.i_pr())
        c_pr = np.exp(lnc_pr)
        f_pr = c_pr*self.p
        return c_pr,f_pr
    #------------------------------------------------------------srk-------------------------------------------------------------------------------------------
    def coeficiente_fug_srk(self):
        lnc_srk = self.z_srk() - 1 - (np.log(self.z_srk() - self.betasrk())) - (self.qsrk()*self.i_srk())
        c_srk = np.exp(lnc_srk)
        f_srk = c_srk*self.p
        return c_srk,f_srk
    #-------------------------------------------------------------rk-------------------------------------------------------------------------------------------
    def coeficiente_fug_rk(self):
        lnc_rk = self.z_rk() - 1 - (np.log(self.z_rk() - self.betark())) - (self.qrk()*self.i_rk())
        c_rk = np.exp(lnc_rk)
        f_rk = c_rk*self.p
        return c_rk,f_rk
    #--------------------------------------------------------vw------------------------------------------------------------------------------------------------
    def coeficiente_fug_vw(self):
        lnc_vw = self.z_vdw() - 1 - (np.log(self.z_vdw() - self.betavw())) - (self.qvw()*self.i_vdw())
        c_vw = np.exp(lnc_vw)
        f_vw = c_vw*self.p
        return c_vw,f_vw
class MisturaBinaria(Prop):
    def __init__(self,y1,y2,tm,pm,componente1,componente2): #componente = [pc,tc,w]
        self.componente1 = componente1
        self.componente2 = componente2
        self.y1 = y1
        self.y2 = y2
        self.tm = tm
        self.pm = pm
    #------------------------------------------------------------zm_pr(quando eu sei a temperatura e pressão da mistura)
    def z_m_pr(self,tol,nmx,kpr):
        componente_1pr = EqState(self.pm,self.tm,self.componente1[0],self.componente1[1],self.componente1[2])
        componente_2pr = EqState(self.pm,self.tm,self.componente2[0],self.componente2[1],self.componente2[2])
        beta1pr = componente_1pr.betapr()
        beta2pr = componente_2pr.betapr()
        q1pr = componente_1pr.qpr()
        q2pr = componente_2pr.qpr()
        betampr = beta1pr*self.y1 + beta2pr*self.y2
        qmpr = ((self.y1**2)*q1pr) + (2*self.y1*self.y2*((q1pr*q2pr)**(0.5))*(1-kpr)) + ((self.y2**2)*q2pr)
        erro = 1
        it = 0
        i = 1 
        while abs(erro)>tol and it<nmx: 
            i = 1 + (betampr) - (qmpr*betampr)*((i-betampr)/((i+(self.espsilonpr()*betampr))*(i+(self.sigmapr()*betampr))))
            erro = erro - i
            it = it + 1
            #print(it,i)
        return i 
    #-----------------------------------------------z_m_srk------------------(quando conheço a temperatura e pressão da mistura)
    def z_m_srk(self,tol,nmx,ksrk):
        componente_1srk = EqState(self.pm,self.tm,self.componente1[0],self.componente1[1],self.componente1[2])
        componente_2srk = EqState(self.pm,self.tm,self.componente2[0],self.componente2[1],self.componente2[2])
        beta1srk = componente_1srk.betasrk()
        beta2srk = componente_2srk.betasrk()
        q1srk = componente_1srk.qsrk()
        q2srk = componente_2srk.qsrk()
        betamsrk = beta1srk*self.y1 + beta2srk*self.y2
        qmsrk = ((self.y1**2)*q1srk) + (2*self.y1*self.y2*((q1srk*q2srk)**(0.5))*(1-ksrk)) + ((self.y2**2)*q2srk)
        erro = 1
        it = 0
        i = 1 
        while abs(erro)>tol and it<nmx: 
            i = 1 + (betamsrk) - (qmsrk*betamsrk)*((i-betamsrk)/((i+(self.espsilonsrk()*betamsrk))*(i+(self.sigmasrk()*betamsrk))))
            erro = erro - i
            it = it + 1
            #print(it,i)
        return i 
    #-------------------------------------------------------z_rk-----(quando sei a temperatura e pressão da mistura)
    def z_m_rk(self,tol,nmx,krk):
        componente_1rk = EqState(self.pm,self.tm,self.componente1[0],self.componente1[1],self.componente1[2])
        componente_2rk = EqState(self.pm,self.tm,self.componente2[0],self.componente2[1],self.componente2[2])
        beta1rk = componente_1rk.betark()
        beta2rk = componente_2rk.betark()
        q1rk = componente_1rk.qrk()
        q2rk = componente_2rk.qrk()
        betamrk = beta1rk*self.y1 + beta2rk*self.y2
        qmrk = ((self.y1**2)*q1rk) + (2*self.y1*self.y2*((q1rk*q2rk)**(0.5))*(1-krk)) + ((self.y2**2)*q2rk)
        erro = 1
        it = 0
        i = 1 
        while abs(erro)>tol and it<nmx: 
            i = 1 + (betamrk) - (qmrk*betamrk)*((i-betamrk)/((i+(self.espsilonrk()*betamrk))*(i+(self.sigmark()*betamrk))))
            erro = erro - i
            it = it + 1
            #print(it,i)
        return i 
    #------------------------------------------------z_m_vw----(quando sei a temperatura e pressão da mistura----------------)
    def z_m_vw(self,tol,nmx,kvw):
        componente_1vw = EqState(self.pm,self.tm,self.componente1[0],self.componente1[1],self.componente1[2])
        componente_2vw = EqState(self.pm,self.tm,self.componente2[0],self.componente2[1],self.componente2[2])
        beta1vw = componente_1vw.betavw()
        beta2vw = componente_2vw.betavw()
        q1vw = componente_1vw.qvw()
        q2vw = componente_2vw.qvw()
        betamvw = beta1vw*self.y1 + beta2vw*self.y2
        qmvw = ((self.y1**2)*q1vw) + (2*self.y1*self.y2*((q1vw*q2vw)**(0.5))*(1-kvw)) + ((self.y2**2)*q2vw)
        erro = 1
        it = 0
        i = 1 
        while abs(erro)>tol and it<nmx: 
            i = 1 + (betamvw) - (qmvw*betamvw)*((i-betamvw)/((i+(self.espsilonvw()*betamvw))*(i+(self.sigmavw()*betamvw))))
            erro = erro - i
            it = it + 1
            #print(it,i)
        return i 
class AuxiliarMist(MisturaBinaria):
    def __init__(self, y1, y2, tm, pm, componente1, componente2):
        super().__init__(y1, y2, tm, pm, componente1, componente2)
    #------------------------------------dalpha_trm_pr_trm(quando eu sei a temperatura e pressão da mistura) - pc,tc,w
    def dalpha_trm_pr_tr(self): # derivade de alpha de pr em função da temperatura reduzida (mistura)
        omega = 0.07780
        psi = 0.45724
        #pcmistura = self.y1*self.componente1[0] + self.y2*self.componente2[0]
        tcmistura = self.y1*self.componente1[1] + self.y2*self.componente2[1]
        wm = self.y1*self.componente1[2] + self.y2*self.componente2[2]
        trm = self.tm/tcmistura
        # h = 0.000001
        # p1 = tr - 2*h
        # p2 = tr - h
        # p3 = tr + h
        # p4 = tr + 2*h
        # fp1 = (((1+(0.37464+1.54226*self.w_c()-0.26992*(self.w_c()**2))*(1-(p1**0.5)))**2)*psi)/(omega*p1)
        # fp2 = (((1+(0.37464+1.54226*self.w_c()-0.26992*(self.w_c()**2))*(1-(p2**0.5)))**2)*psi)/(omega*p2)
        # fp3 = (((1+(0.37464+1.54226*self.w_c()-0.26992*(self.w_c()**2))*(1-(p3**0.5)))**2)*psi)/(omega*p3)
        # fp4 = (((1+(0.37464+1.54226*self.w_c()-0.26992*(self.w_c()**2))*(1-(p4**0.5)))**2)*psi)/(omega*p4)
        # d = (fp1 - 8*fp2 + 8*fp3 - fp4)/(12*h)
        def dpr(trm):
            return (((1+(0.37464+1.54226*wm-0.26992*(wm**2))*(1-(trm**0.5)))**2)*psi)/(omega*trm)
        t0_pr = np.array([trm])
        dx_pr = 1e-6
        derivadapr = approx_fprime(t0_pr,dpr,dx_pr)
        return derivadapr[0]
    #---------------------------------------------dalpha_trm_srk_trm--(quando eu sei a temperatura e pressão da mistura)
    def dalpha_trm_srk_tr(self): # derivada de alpha de srk em função da temperatura reduzida (mistura)
        omega = 0.08664
        psi = 0.42748
        tcmistura = self.y1*self.componente1[1] + self.y2*self.componente2[1]
        wm = self.y1*self.componente1[2] + self.y2*self.componente2[2]
        trm = self.tm/tcmistura
        h = 0.000001
        p1 = trm - 2*h
        p2 = trm - h
        p3 = trm + h
        p4 = trm + 2*h
        fp1 = (((1+(0.480+1.574*wm-0.176*(wm**2))*(1-(p1**0.5)))**2)*psi)/(omega*p1)
        fp2 = (((1+(0.480+1.574*wm-0.176*(wm**2))*(1-(p2**0.5)))**2)*psi)/(omega*p2)
        fp3 = (((1+(0.480+1.574*wm-0.176*(wm**2))*(1-(p3**0.5)))**2)*psi)/(omega*p3)
        fp4 = (((1+(0.480+1.574*wm-0.176*(wm**2))*(1-(p4**0.5)))**2)*psi)/(omega*p4)
        d = (fp1 - 8*fp2 + 8*fp3 - fp4)/(12*h)
        return d
    #----------------------------dalpha_trm_rk_trm----(quando eu sei a temperatura e pressão da mistura)
    def dalpha_trm_rk_tr(self): # derivade de alpha de rk em função da temperatura reduzida (mistura)
        omega = 0.08664
        psi = 0.42748 
        tcmistura = self.y1*self.componente1[1] + self.y2*self.componente2[1]
        #wm = self.y1*self.componente1[2] + self.y2*self.componente2[2]
        trm = self.tm/tcmistura
        h = 0.000001
        p1 = trm - 2*h
        p2 = trm - h
        p3 = trm + h
        p4 = trm + 2*h
        fp1 = ((p1**(-0.5))*psi)/(omega*p1)
        fp2 = ((p2**(-0.5))*psi)/(omega*p2)
        fp3 = ((p3**(-0.5))*psi)/(omega*p3)
        fp4 = ((p4**(-0.5))*psi)/(omega*p4)
        d = (fp1 - 8*fp2 + 8*fp3 - fp4)/(12*h)
        return d
    #----------------------------------dalpha_trm_vw_trm---(quando eu sei a temperatura e pressão da mistura)
    def dalpha_trm_vdw_tr(self): # derivade de alpha de vw em função da temperatura reduzida (mistura)
        omega = 1/8
        psi = 27/64
        tcmistura = self.y1*self.componente1[1] + self.y2*self.componente2[1]
        #wm = self.y1*self.componente1[2] + self.y2*self.componente2[2]
        trm = self.tm/tcmistura
        h = 0.000001
        p1 = trm - 2*h
        p2 = trm - h
        p3 = trm + h
        p4 = trm + 2*h
        fp1 = (psi)/(omega*p1)
        fp2 = (psi)/(omega*p2)
        fp3 = (psi)/(omega*p3)
        fp4 = (psi)/(omega*p4)
        d = (fp1 - 8*fp2 + 8*fp3 - fp4)/(12*h)
        return d
    def i_prm(self,tol,nmx,kpr): # retorna o parâmetro i de pr (mistura) - conheço t e p da mistura 
        componente_1pr = EqState(self.pm,self.tm,self.componente1[0],self.componente1[1],self.componente1[2])
        componente_2pr = EqState(self.pm,self.tm,self.componente2[0],self.componente2[1],self.componente2[2])
        beta1pr = componente_1pr.betapr()
        beta2pr = componente_2pr.betapr()
        betampr = beta1pr*self.y1 + beta2pr*self.y2
        ipr = (1/(self.sigmapr() - self.espsilonpr()))*(np.log((self.z_m_pr(tol,nmx,kpr)+(self.sigmapr()*betampr))/(self.z_m_pr(tol,nmx,kpr)+(self.espsilonpr()*betampr))))
        return ipr 
    def i_srkm(self,tol,nmx,ksrk): # retorna o parâmetro i de srk (mistura) - conheço t e p da mistura
        componente_1srk = EqState(self.pm,self.tm,self.componente1[0],self.componente1[1],self.componente1[2])
        componente_2srk = EqState(self.pm,self.tm,self.componente2[0],self.componente2[1],self.componente2[2])
        beta1srk = componente_1srk.betasrk()
        beta2srk = componente_2srk.betasrk()
        betamsrk = beta1srk*self.y1 + beta2srk*self.y2
        isrk = (1/(self.sigmasrk() - self.espsilonsrk()))*(np.log((self.z_m_srk(tol,nmx,ksrk)+(self.sigmasrk()*betamsrk))/(self.z_m_srk(tol,nmx,ksrk)+(self.espsilonsrk()*betamsrk))))
        return isrk
    def i_rkm(self,tol,nmx,krk): # retorna o parâmetro i de rk (mistura) - conheço t e p da mistura
        componente_1rk = EqState(self.pm,self.tm,self.componente1[0],self.componente1[1],self.componente1[2])
        componente_2rk = EqState(self.pm,self.tm,self.componente2[0],self.componente2[1],self.componente2[2])
        beta1rk = componente_1rk.betark()
        beta2rk = componente_2rk.betark()
        betamrk = beta1rk*self.y1 + beta2rk*self.y2
        irk = (1/(self.sigmark() - self.espsilonrk()))*(np.log((self.z_m_rk(tol,nmx,krk)+(self.sigmark()*betamrk))/(self.z_m_rk(tol,nmx,krk)+(self.espsilonrk()*betamrk))))
        return irk
    def i_vdwm(self,tol,nmx,kvw): # retorna o parâmetro i de vdw (mistura) - conheço t e p da mistura
        componente_1vw = EqState(self.pm,self.tm,self.componente1[0],self.componente1[1],self.componente1[2])
        componente_2vw = EqState(self.pm,self.tm,self.componente2[0],self.componente2[1],self.componente2[2])
        beta1vw = componente_1vw.betavw()
        beta2vw = componente_2vw.betavw()
        betamvw = beta1vw*self.y1 + beta2vw*self.y2
        ivdw = (betamvw)/(self.z_m_vw(tol,nmx,kvw)+(self.espsilonvw()*betamvw))
        return ivdw
class PropMist(AuxiliarMist):
    def __init__(self, y1, y2, tm, pm, componente1, componente2):
        super().__init__(y1, y2, tm, pm, componente1, componente2)
    def cpm_gasid_formacao_298k(self,cp1,cp2): #h(298.15k)=0 -- cp = [a,b,c,d]
        am = self.y1*cp1[0] + self.y2*cp2[0]
        bm = self.y1*cp1[1] + self.y2*cp2[1]
        cm = self.y1*cp1[2] + self.y2*cp2[2]
        dm = self.y1*cp1[3] + self.y2*cp2[3]
        cpm = self.cp_gasid_formacao_298k(am,bm,cm,dm)
        return cpm
    def entalpiam_gas_ideal_298k(self,cp1,cp2):
        hmid = self.cpm_gasid_formacao_298k(cp1,cp2)
        return hmid
#propriedades da mistura usando o método de PR(PENG-ROBINSON)
    def entalpia_residualm_pr(self,tol,nmx,kpr):
        r = 8.314 # kJ/kmol.K
        tcmistura = self.y1*self.componente1[1] + self.y2*self.componente2[1]
        tr = self.tm/tcmistura
        hr_pr = (self.z_m_pr(tol,nmx,kpr) - 1 + (tr*self.dalpha_trm_pr_tr()*self.i_prm(tol,nmx,kpr)))*r*self.tm
        return hr_pr
    def entalpia_absolutam_pr(self,cp1,cp2,tol,nmx,kpr):
        h_absolute  = self.entalpiam_gas_ideal_298k(cp1,cp2)+self.entalpia_residualm_pr(tol,nmx,kpr)
        return h_absolute
    def entropia_residualm_pr(self,tol,nmx,kpr):
        r = 8.314 # kJ/kmol.K
        tcmistura = self.y1*self.componente1[1] + self.y2*self.componente2[1]
        tr = self.tm/tcmistura
        componente_1pr = EqState(self.pm,self.tm,self.componente1[0],self.componente1[1],self.componente1[2])
        componente_2pr = EqState(self.pm,self.tm,self.componente2[0],self.componente2[1],self.componente2[2])
        beta1pr = componente_1pr.betapr()
        beta2pr = componente_2pr.betapr()
        betampr = beta1pr*self.y1 + beta2pr*self.y2
        q1pr = componente_1pr.qpr()
        q2pr = componente_2pr.qpr()
        qmpr = ((self.y1**2)*q1pr) + (2*self.y1*self.y2*((q1pr*q2pr)**(0.5))*(1-kpr)) + ((self.y2**2)*q2pr)
        s_res_pr = ((np.log(self.z_m_pr(tol,nmx,kpr) - betampr)) + ((qmpr+(tr*self.dalpha_trm_pr_tr()))*(self.i_prm(tol,nmx,kpr))))*r
        return s_res_pr
    def entropia_gasidealm(self,cp1,cp2):
        r = 8.314 # kJ/kmol.K
        am = self.y1*cp1[0] + self.y2*cp2[0]
        bm = self.y1*cp1[1] + self.y2*cp2[1]
        cm = self.y1*cp1[2] + self.y2*cp2[2]
        dm = self.y1*cp1[3] + self.y2*cp2[3]
        sid = am*(np.log((self.tm)/(298.15))) + bm*(self.tm - 298.15) + (cm/2)*(((self.tm)**2) - (298.15**2)) + (dm/3)*(((self.tm)**3) - (298.15**3)) - r*np.log(self.pm)
        return sid
    def entropia_absolutam_pr(self,cp1,cp2,tol,nmx,kpr):
        s_absolute = self.entropia_gasidealm(cp1,cp2)+self.entropia_residualm_pr(tol,nmx,kpr)
        return s_absolute
    def gibbs_absolutem_pr(self,cp1,cp2,tol,nmx,kpr):
        gibbs_absolute = self.entalpia_absolutam_pr(cp1,cp2,tol,nmx,kpr) - ((self.tm)*self.entropia_absolutam_pr(cp1,cp2,tol,nmx,kpr))
        return gibbs_absolute
    #------------------------------------------interna_pr-------------------------------------------------------------------------------------------------
    def interna_absolutem_pr(self,cp1,cp2,tol,nmx,kpr):
        r = 8.314 # kJ/kmol.K
        interna_absolute = self.entalpia_absolutam_pr(cp1,cp2,tol,nmx,kpr) - (self.z_m_pr(tol,nmx,kpr)*r*self.tm)
        return interna_absolute
    #-----------------------------------------helmholtz_pr------------------------------------------------------------------------------------------------
    def helmholtz_absolutem_pr(self,cp1,cp2,tol,nmx,kpr):
        helmholtz_absolute = self.interna_absolutem_pr(cp1,cp2,tol,nmx,kpr) - (self.tm*self.entropia_absolutam_pr(cp1,cp2,tol,nmx,kpr))
        return helmholtz_absolute
    #propriedades da mistura usando o método SRK-----------------------------------------------------------------------------------------------------------------------
    def entalpia_residualm_srk(self,tol,nmx,ksrk):
        r = 8.314 # kJ/kmol.K
        tcmistura = self.y1*self.componente1[1] + self.y2*self.componente2[1]
        tr = self.tm/tcmistura
        hr_pr = (self.z_m_srk(tol,nmx,ksrk) - 1 + (tr*self.dalpha_trm_srk_tr()*self.i_srkm(tol,nmx,ksrk)))*r*self.tm
        return hr_pr
    def entalpia_absolutam_srk(self,cp1,cp2,tol,nmx,ksrk):
        h_absolute  = self.entalpiam_gas_ideal_298k(cp1,cp2)+self.entalpia_residualm_srk(tol,nmx,ksrk)
        return h_absolute
    def entropia_residualm_srk(self,tol,nmx,ksrk):
        r = 8.314 # kJ/kmol.K
        tcmistura = self.y1*self.componente1[1] + self.y2*self.componente2[1]
        tr = self.tm/tcmistura
        componente_1srk = EqState(self.pm,self.tm,self.componente1[0],self.componente1[1],self.componente1[2])
        componente_2srk = EqState(self.pm,self.tm,self.componente2[0],self.componente2[1],self.componente2[2])
        beta1srk = componente_1srk.betasrk()
        beta2srk = componente_2srk.betasrk()
        betamsrk = beta1srk*self.y1 + beta2srk*self.y2
        q1srk = componente_1srk.qsrk()
        q2srk = componente_2srk.qsrk()
        qmsrk = ((self.y1**2)*q1srk) + (2*self.y1*self.y2*((q1srk*q2srk)**(0.5))*(1-ksrk)) + ((self.y2**2)*q2srk)
        s_res_srk = ((np.log(self.z_m_srk(tol,nmx,ksrk) - betamsrk)) + ((qmsrk+(tr*self.dalpha_trm_srk_tr()))*(self.i_srkm(tol,nmx,ksrk))))*r
        return s_res_srk
    def entropia_absolutam_srk(self,cp1,cp2,tol,nmx,ksrk):
        s_absolute = self.entropia_gasidealm(cp1,cp2)+self.entropia_residualm_srk(tol,nmx,ksrk)
        return s_absolute
    def gibbs_absolutem_srk(self,cp1,cp2,tol,nmx,ksrk):
        gibbs_absolute = self.entalpia_absolutam_srk(cp1,cp2,tol,nmx,ksrk) - ((self.tm)*self.entropia_absolutam_srk(cp1,cp2,tol,nmx,ksrk))
        return gibbs_absolute
    #------------------------------------------interna_srk-------------------------------------------------------------------------------------------------
    def interna_absolutem_srk(self,cp1,cp2,tol,nmx,ksrk):
        r = 8.314 # kJ/kmol.K
        interna_absolute = self.entalpia_absolutam_srk(cp1,cp2,tol,nmx,ksrk) - (self.z_m_srk(tol,nmx,ksrk)*r*self.tm)
        return interna_absolute
    #-----------------------------------------helmholtz_srk------------------------------------------------------------------------------------------------
    def helmholtz_absolutem_srk(self,cp1,cp2,tol,nmx,ksrk):
        helmholtz_absolute = self.interna_absolutem_srk(cp1,cp2,tol,nmx,ksrk) - (self.tm*self.entropia_absolutam_srk(cp1,cp2,tol,nmx,ksrk))
        return helmholtz_absolute
    #propriedades de mistura usando o rk 
    def entalpia_residualm_rk(self,tol,nmx,krk):
        r = 8.314 # kJ/kmol.K
        tcmistura = self.y1*self.componente1[1] + self.y2*self.componente2[1]
        tr = self.tm/tcmistura
        hr_rk = (self.z_m_rk(tol,nmx,krk) - 1 + (tr*self.dalpha_trm_rk_tr()*self.i_rkm(tol,nmx,krk)))*r*self.tm
        return hr_rk
    def entalpia_absolutam_rk(self,cp1,cp2,tol,nmx,krk):
        h_absolute  = self.entalpiam_gas_ideal_298k(cp1,cp2)+self.entalpia_residualm_rk(tol,nmx,krk)
        return h_absolute
    def entropia_residualm_rk(self,tol,nmx,krk):
        r = 8.314 # kJ/kmol.K
        tcmistura = self.y1*self.componente1[1] + self.y2*self.componente2[1]
        tr = self.tm/tcmistura
        componente_1rk = EqState(self.pm,self.tm,self.componente1[0],self.componente1[1],self.componente1[2])
        componente_2rk = EqState(self.pm,self.tm,self.componente2[0],self.componente2[1],self.componente2[2])
        beta1rk = componente_1rk.betark()
        beta2rk = componente_2rk.betark()
        betamrk = beta1rk*self.y1 + beta2rk*self.y2
        q1rk = componente_1rk.qrk()
        q2rk = componente_2rk.qrk()
        qmrk = ((self.y1**2)*q1rk) + (2*self.y1*self.y2*((q1rk*q2rk)**(0.5))*(1-krk)) + ((self.y2**2)*q2rk)
        s_res_rk = ((np.log(self.z_m_rk(tol,nmx,krk) - betamrk)) + ((qmrk+(tr*self.dalpha_trm_rk_tr()))*(self.i_rkm(tol,nmx,krk))))*r
        return s_res_rk
    def entropia_absolutam_rk(self,cp1,cp2,tol,nmx,krk):
        s_absolute = self.entropia_gasidealm(cp1,cp2)+self.entropia_residualm_rk(tol,nmx,krk)
        return s_absolute
    def gibbs_absolutem_rk(self,cp1,cp2,tol,nmx,krk):
        gibbs_absolute = self.entalpia_absolutam_rk(cp1,cp2,tol,nmx,krk) - ((self.tm)*self.entropia_absolutam_rk(cp1,cp2,tol,nmx,krk))
        return gibbs_absolute
    #------------------------------------------interna_rk-------------------------------------------------------------------------------------------------
    def interna_absolutem_rk(self,cp1,cp2,tol,nmx,krk):
        r = 8.314 # kJ/kmol.K
        interna_absolute = self.entalpia_absolutam_rk(cp1,cp2,tol,nmx,krk) - (self.z_m_rk(tol,nmx,krk)*r*self.tm)
        return interna_absolute
    #-----------------------------------------helmholtz_rk------------------------------------------------------------------------------------------------
    def helmholtz_absolutem_rk(self,cp1,cp2,tol,nmx,krk):
        helmholtz_absolute = self.interna_absolutem_rk(cp1,cp2,tol,nmx,krk) - (self.tm*self.entropia_absolutam_rk(cp1,cp2,tol,nmx,krk))
        return helmholtz_absolute
    #propriedades de mistura usando vw---------------------------------------------------------------------------------------------------------------------
    def entalpia_residualm_vw(self,tol,nmx,kvw):
        r = 8.314 # kJ/kmol.K
        tcmistura = self.y1*self.componente1[1] + self.y2*self.componente2[1]
        tr = self.tm/tcmistura
        hr_vw = (self.z_m_vw(tol,nmx,kvw) - 1 + (tr*self.dalpha_trm_vw_tr()*self.i_vdwm(tol,nmx,kvw)))*r*self.tm
        return hr_vw
    def entalpia_absolutam_vw(self,cp1,cp2,tol,nmx,kvw):
        h_absolute  = self.entalpiam_gas_ideal_298k(cp1,cp2)+self.entalpia_residualm_vw(tol,nmx,kvw)
        return h_absolute
    def entropia_residualm_vw(self,tol,nmx,kvw):
        r = 8.314 # kJ/kmol.K
        tcmistura = self.y1*self.componente1[1] + self.y2*self.componente2[1]
        tr = self.tm/tcmistura
        componente_1vw = EqState(self.pm,self.tm,self.componente1[0],self.componente1[1],self.componente1[2])
        componente_2vw = EqState(self.pm,self.tm,self.componente2[0],self.componente2[1],self.componente2[2])
        beta1vw = componente_1vw.betavw()
        beta2vw = componente_2vw.betavw()
        betamvw = beta1vw*self.y1 + beta2vw*self.y2
        q1vw = componente_1vw.qvw()
        q2vw = componente_2vw.qvw()
        qmvw = ((self.y1**2)*q1vw) + (2*self.y1*self.y2*((q1vw*q2vw)**(0.5))*(1-kvw)) + ((self.y2**2)*q2vw)
        s_res_vw = ((np.log(self.z_m_vw(tol,nmx,kvw) - betamvw)) + ((qmvw+(tr*self.dalpha_trm_vdw_tr()))*(self.i_vdwm(tol,nmx,kvw))))*r
        return s_res_vw
    def entropia_absolutam_vw(self,cp1,cp2,tol,nmx,kvw):
        s_absolute = self.entropia_gasidealm(cp1,cp2)+self.entropia_residualm_vw(tol,nmx,kvw)
        return s_absolute
    def gibbs_absolutem_vw(self,cp1,cp2,tol,nmx,kvw):
        gibbs_absolute = self.entalpia_absolutam_vw(cp1,cp2,tol,nmx,kvw) - ((self.tm)*self.entropia_absolutam_vw(cp1,cp2,tol,nmx,kvw))
        return gibbs_absolute
    #------------------------------------------interna_vdw-------------------------------------------------------------------------------------------------
    def interna_absolutem_rk(self,cp1,cp2,tol,nmx,krk):
        r = 8.314 # kJ/kmol.K
        interna_absolute = self.entalpia_absolutam_rk(cp1,cp2,tol,nmx,krk) - (self.z_m_rk(tol,nmx,krk)*r*self.tm)
        return interna_absolute
    #-----------------------------------------helmholtz_vdw------------------------------------------------------------------------------------------------
    def helmholtz_absolutem_rk(self,cp1,cp2,tol,nmx,krk):
        helmholtz_absolute = self.interna_absolutem_rk(cp1,cp2,tol,nmx,krk) - (self.tm*self.entropia_absolutam_rk(cp1,cp2,tol,nmx,krk))
        return helmholtz_absolute



    
'''
n_butano  = EqState(9.4573,350,37.96,425.1,0.2)
n_butano_vap = vapor(9.4573,350,37.96,425.1,0.2,300,0.0000001)
nbutanoprop = Prop(50,500,37.96,425.1,0.2,300,0.0000001)
zpr = n_butano_vap.z_pr()
vpr = n_butano_vap.v_pr()
zvw = n_butano_vap.z_vdw()
vvdw  = n_butano_vap.v_vdw()
zrk = n_butano_vap.z_rk()
vrk  = n_butano_vap.v_rk()
zsrk = n_butano_vap.z_srk()
vsrk  = n_butano_vap.v_srk()
hr = nbutanoprop.entalpia_residual_pr()
hid = nbutanoprop.entalpia_gas_ideal_formacao_models(3.96,0.3715,-0.0001834,0.000000035)
cp = nbutanoprop.cp_gasid_formacao_298k(3.96,0.3715,-0.0001834,0.000000035)
habsulute = nbutanoprop.entalpia_absolute_pr(3.96,0.3715,-0.0001834,0.000000035)
sr = nbutanoprop.entropia_residual_pr()
entropiaid = nbutanoprop.entropia_gas_ideal_formacao(3.96,0.3715,-0.0001834,0.000000035)
entropiaabsolute = nbutanoprop.entropia_absolute_pr(3.96,0.3715,-0.0001834,0.000000035)
gibbsabs = nbutanoprop.gibbs_absolute_pr(3.96,0.3715,-0.0001834,0.000000035)
internaabs = nbutanoprop.interna_absolute_pr(3.96,0.3715,-0.0001834,0.000000035)
helmabs = nbutanoprop.helmholtz_absolute_pr(3.96,0.3715,-0.0001834,0.000000035)
print(zpr,vpr)
print(zvw,vvdw)
print(zrk,vrk)
print(zsrk,vsrk)
print(hr)
print(hid)
print(cp)
print(habsulute)
print(sr)
print(entropiaid)
print(entropiaabsolute)
print(gibbsabs)
print(internaabs)
print(helmabs)
print()
entalpiasrk = nbutanoprop.entalpia_absolute_srk(3.96,0.3715,-0.0001834,0.000000035)
entropiasrk = nbutanoprop.entropia_absolute_srk(3.96,0.3715,-0.0001834,0.000000035)
gibbssrk = nbutanoprop.gibbs_absolute_srk(3.96,0.3715,-0.0001834,0.000000035)
internasrk = nbutanoprop.interna_absolute_srk(3.96,0.3715,-0.0001834,0.000000035)
helmoltzsrk = nbutanoprop.helmholtz_absolute_srk(3.96,0.3715,-0.0001834,0.000000035)
print(entalpiasrk)
print(entropiasrk)
print(gibbssrk)
print(internasrk)
print(helmoltzsrk)
print()
entalpiark = nbutanoprop.entalpia_absolute_rk(3.96,0.3715,-0.0001834,0.000000035)
enttropiark = nbutanoprop.entropia_absolute_rk(3.96,0.3715,-0.0001834,0.000000035)
gibbsrk = nbutanoprop.gibbs_absolute_rk(3.96,0.3715,-0.0001834,0.000000035)
internark = nbutanoprop.interna_absolute_rk(3.96,0.3715,-0.0001834,0.000000035)
helmoltzrk = nbutanoprop.helmholtz_absolute_rk(3.96,0.3715,-0.0001834,0.000000035)
print(entalpiark)
print(enttropiark)
print(gibbsrk)
print(internark)
print(helmoltzrk)
print()
entalpiavw = nbutanoprop.entalpia_absolute_vw(3.96,0.3715,-0.0001834,0.000000035)
enttropiavw = nbutanoprop.entropia_absolute_vw(3.96,0.3715,-0.0001834,0.000000035)
gibbsvw = nbutanoprop.gibbs_absolute_vw(3.96,0.3715,-0.0001834,0.000000035)
internavw = nbutanoprop.interna_absolute_vw(3.96,0.3715,-0.0001834,0.000000035)
helmoltzvw = nbutanoprop.helmholtz_absolute_vw(3.96,0.3715,-0.0001834,0.000000035)
print(entalpiavw)
print(enttropiavw)
print(gibbsvw)
print(internavw)
print(helmoltzvw)
print()
fugpr = nbutanoprop.coeficiente_fug_pr()
fugsrk = nbutanoprop.coeficiente_fug_srk()
fugrk = nbutanoprop.coeficiente_fug_rk()
fugvw = nbutanoprop.coeficiente_fug_vw()
print(fugpr[0],fugpr[1])
print(fugsrk[0],fugsrk[1])
print(fugrk[0],fugrk[1])
print(fugvw[0],fugvw[1])
print()
butano = [37.96,425.1,0.2]
propano = [42.49,369.82,0.152]
misturabp = MisturaBinaria(0.6,0.4,550,55,butano,propano)
zmpr = misturabp.z_m_pr(0.00001,200,0)
zmsrk = misturabp.z_m_srk(0.00001,200,0)
zmrk = misturabp.z_m_rk(0.00001,200,0)
zmvw = misturabp.z_m_vw(0.00001,200,0)
print(zmpr)
print(zmsrk)
print(zmrk)
print(zmvw)
print(nbutanoprop.z_pr())
'''


    
    
        