#load modules
import datetime as dt
import numpy as np
import SEB_functions as SEBf
import matplotlib.pyplot as plt

#load data
station = "AWS17"   #AWS15 and AWS17 work best
dtime = "daily"     #'daily' or 'hourly'
# Options Greenland:
# S5 27-08-2003 - 24-08-2023
# S6 29-08-2003 - 21-08-2023
# S9 23-08-2000 - 21-08-2023
# S10 17-08-2010 - 17-04-2026
# Options Antarctica:
# AWS14 21-01-2009 - 01-01-2023
# AWS15 21-09-2009 - 24-06-2014
# AWS17 19-02-2011 - 10-03-2016
# AWS18 25-11-2014 - 02-12-2022
if station == "S5" or station == "S6" or station == "S9" or station == "S10":
    SEBdata = SEBf.SEB_data(FileName="Pangaea-Data/GRL_"+ station +"_"+dtime+"_all.csv")
elif station == "AWS14" or station == "AWS15" or station == "AWS17" or station == "AWS18":     
    SEBdata = SEBf.SEB_data(FileName="Pangaea-Data/ANT_"+ station +"_"+dtime+"_all.csv")

# Extracting the variables from the data
SWdown = SEBdata.Extract_Variable("SWd") # from observations
SWup   = SEBdata.Extract_Variable("SWu") # from observations
SWnet  = SWdown - SWup
LWdown = SEBdata.Extract_Variable("LWd") # from observations
LWup   = SEBdata.Extract_Variable("LWu") # from observations
LWup_mod = SEBdata.Extract_Variable("LWu_mod") # from modeled TS
LWnet  = LWdown - LWup_mod
SHF    = SEBdata.Extract_Variable("SHFdown_mod") # calculated with MO using modeled TS
LHF    = SEBdata.Extract_Variable("LHFdown_mod") # calculated with MO using modeled TS
Gs     = SEBdata.Extract_Variable("GHFup_mod") # calculated using modelled TS
MeltS  = SEBdata.Extract_Variable("meltE") # calculated using modelled TS
MeltT  = MeltS  #SEBdata.Extract_Variable("totm_nrg") #not available is melt including subsurface when implementing radiation penetration
Residu = SWnet + LWnet + SHF + LHF + Gs - MeltS # residu combines melt with non closure of energy terms
T2m    = SEBdata.Extract_Variable("t2m")
pressure=SEBdata.Extract_Variable("p")
ff10m   = SEBdata.Extract_Variable("ff10m")
cloud_cover = SEBdata.Extract_Variable("cloud_cover")
Ts_mod = SEBdata.Extract_Variable("Ts_mod")
c_s=SHF/(ff10m*(T2m-Ts_mod)+1e-6) #correction factor of 1e-6 to avoid division by zero

########################################################################################
#adapt conditions for our scenario: Always cloudy
#variable cloud_cover: 0=clear sky, 1=completely cloudy

#use values from other script: Check if we can set some of them to zero
corr_SWdown = -0.19
corr_SWup = -0.15
corr_LWdown = 0.43
corr_LWup = -0.03
corr_SHF = -0.50
corr_LHF = -0.11
corr_Gs = -0.08

#correct the values for cloud cover: Multiply correlation with increase in cloud cover and add to original value
cloud_cover=np.nan_to_num(cloud_cover,nan=1.0)
#Explanation: if there is no measurement for cloud cover for the day: Set to 1, so that other variables unchanged
SWdown=SWdown*(1+(1-cloud_cover)*corr_SWdown)
SWup=SWup*(1+(1-cloud_cover)*corr_SWup)
LWdown=LWdown*(1+(1-cloud_cover)*corr_LWdown)
LWup=LWup*(1+(1-cloud_cover)*corr_LWup)
SHF=SHF*(1+(1-cloud_cover)*corr_SHF)
LHF=LHF*(1+(1-cloud_cover)*corr_LHF)
Gs=Gs*(1+(1-cloud_cover)*corr_Gs)



########################################################################################
#from here: Solve equation

#variables to track success/failures of the root finding algorithm
fail_bracket=0
fail_maxiter=0
fail_invalid_inputdata=0
success=0

#SEB equation. Can include or leave out melt and residual terms
def SEB_equation(Ts_mod,SWdown,SWup,LWdown,ff10m,T2m,LHF,Gs,MeltS,Residu,c_s,incl_melt=False,incl_residual=False):
    sigma = 5.67*10**-8
    LHS=0
    if incl_melt:
        LHS+=MeltS
    if incl_residual:
        LHS+=Residu
    RHS=SWdown-SWup+LWdown-sigma*(Ts_mod+273.16)**4+c_s*ff10m*(T2m-Ts_mod)+LHF+Gs
    return RHS-LHS

#bisection method for root finding
def find_root(func,x1,x2,args,tol=1e-6,N_max=100):
    global fail_bracket, fail_maxiter, success
    #check if the initial values are correct
    if (func(x1,*args)*func(x2,*args))>0:
        #print("Root not found: Initial values do not bracket the root")
        fail_bracket+=1
        return None
    #iterate
    i=0
    while i<N_max:
        x=(x1+x2)/2.0
        if func(x,*args)==0 or np.abs((x2-x1)/2.0)<tol:
            #print("Root found after",i,"iterations")
            success+=1
            return x
        if np.sign(func(x,*args))==np.sign(func(x1,*args)):
            x1=x
        if np.sign(func(x,*args))==np.sign(func(x2,*args)):
            x2=x
        i+=1
    print("Root not found: Maximum number of iterations reached")
    fail_maxiter+=1
    return None

#function to solve the SEB equation and retrun the the new T and M
def close_equation(func,T_init_1,T_init_2,args):
    new_temp=find_root(func,T_init_1,T_init_2,args)
    if new_temp is None:
        #print("No solution found")
        return np.nan, np.nan
    if new_temp>0:
        #print("Temperature is above zero degrees Celsius")
        new_temp=0
        melt=func(new_temp,*args)
        #print('T = 0°C')
        #print(f'M = {melt:.3f} W/m^2')
    else:
        melt=0
        #print(f'T = {new_temp:.3f}°C')
        #print('M = 0 W/m^2')
    return new_temp, melt

#solve the SEB equation for each given time step
T_solve=[]
M_solve=[]
for i in range(len(SWdown)):
    if (np.isnan(SWdown[i]) or np.isnan(SWup[i]) or np.isnan(LWdown[i]) or
        np.isnan(ff10m[i]) or np.isnan(T2m[i]) or np.isnan(LHF[i]) or
        np.isnan(Gs[i]) or np.isnan(MeltS[i]) or np.isnan(Residu[i]) or
        np.isnan(c_s[i])):
        fail_invalid_inputdata+=1
        T_solve.append(np.nan)
        M_solve.append(np.nan)
        continue    #skip iteration if any of the values is NaN
    args=[SWdown[i],SWup[i],LWdown[i],ff10m[i],T2m[i],LHF[i],Gs[i],MeltS[i],Residu[i],c_s[i]]
    T1_init,T2_init=Ts_mod[i]+50,Ts_mod[i]-50
    T,M=close_equation(SEB_equation,T1_init,T2_init,args)
    T_solve.append(T)
    M_solve.append(M)
T_solve=np.array(T_solve)
M_solve=np.array(M_solve)
print(f"Bracket failure: {fail_bracket}")
print(f"Maxiter failure: {fail_maxiter}")
print(f"Invalid input data: {fail_invalid_inputdata}")
print(f"Success: {success}")



##########################################################################
#from here: Plot the results

#change data arrays to daily/monthly/avgmonthly for plotting. Also, create the according time data for x-axis
timestep='avgmonthly' # 'daily', 'monthly', 'avgmonthly'
if timestep=='daily':
    f=SEBf.get_daily_average
    Time_data  = f(SHF, SEBdata.DateTime, GiveDatesBack=True)[1]
    Time_range  = [dt.datetime.fromisoformat("2012-01-01"), dt.datetime.fromisoformat("2013-01-01")]
    Time_Label  = "Date"
elif timestep=='monthly':
    f=SEBf.get_monthly_average
    Time_data  = f(SHF, SEBdata.DateTime, GiveDatesBack=True)[1]
    Time_range  = [SEBdata.DateTime[0], SEBdata.DateTime[-1]]
    Time_Label  = "Year"
elif timestep=='avgmonthly':
    f=SEBf.get_avg_monthly_value
    Time_data  = np.arange(14)
    Time_range  = [0, 12.5]
    Time_Label  = "Month"
SWdown=f(SWdown,SEBdata.DateTime)
SWup=f(SWup,SEBdata.DateTime)
LWdown=f(LWdown,SEBdata.DateTime)
ff10m=f(ff10m,SEBdata.DateTime)
T2m=f(T2m,SEBdata.DateTime)
LHF=f(LHF,SEBdata.DateTime)
Gs=f(Gs,SEBdata.DateTime)
MeltS=f(MeltS,SEBdata.DateTime)
Residu=f(Residu,SEBdata.DateTime)
c_s=f(c_s,SEBdata.DateTime)
Ts_mod=f(Ts_mod,SEBdata.DateTime)
cloud_cover=f(cloud_cover,SEBdata.DateTime)
T_solve=f(T_solve,SEBdata.DateTime)
M_solve=f(M_solve,SEBdata.DateTime)

#find NaN values and remove them
indices_to_remove = ~np.isnan(T_solve)
SWdown=SWdown[indices_to_remove]
SWup=SWup[indices_to_remove]
LWdown=LWdown[indices_to_remove]
ff10m=ff10m[indices_to_remove]
T2m=T2m[indices_to_remove]
LHF=LHF[indices_to_remove]
Gs=Gs[indices_to_remove]
MeltS=MeltS[indices_to_remove]
Residu=Residu[indices_to_remove]
c_s=c_s[indices_to_remove]
Ts_mod=Ts_mod[indices_to_remove]
T_solve=T_solve[indices_to_remove]
M_solve=M_solve[indices_to_remove]
Time_data=Time_data[indices_to_remove]

#print some results
print('Change in total melt:',np.sum(M_solve)-np.sum(MeltS),'W/m^2')
print('Average increase in surface Temperature:',np.mean(T_solve-Ts_mod),'°C')


#compare visually the solved M,T with the original values from the data
fig, axs = plt.subplots(2,sharex=True)
fig.suptitle(f'Station: {station}, timestep:{timestep} ')
axs[0].set_ylabel('Temperature [°C]')
axs[0].plot(Time_data,Ts_mod,label='Ts_mod')
axs[0].plot(Time_data,T_solve,label='T_solve')
axs[0].legend()
axs[0].grid(True)
axs[1].set_ylabel(r'Melt [$\frac{W}{m^2}$]')
axs[1].plot(Time_data,MeltS,label='MeltS')
axs[1].plot(Time_data,M_solve,label='M_solve')
axs[1].legend()
axs[1].grid(True)
axs[1].set_xlabel(Time_Label)
#axs[1].set_xlim(Time_range)
plt.show()



