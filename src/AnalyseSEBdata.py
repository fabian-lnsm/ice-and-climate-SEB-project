# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 16:37:18 2023

@author: berg
"""

import datetime as dt
import numpy as np
import SEB_functions as SEBf
import matplotlib.pyplot as plt

station = "AWS14"
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


#%%
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

#calculate the coefficient c_s
SHF_avg_monthly=SEBf.get_avg_monthly_value(SHF, SEBdata.DateTime)
T2m_avg_monthly=SEBf.get_avg_monthly_value(T2m, SEBdata.DateTime)
Tsurface_avg_monthly=SEBf.get_avg_monthly_value(Ts_mod, SEBdata.DateTime)
U10m_avg_monthly=SEBf.get_avg_monthly_value(ff10m, SEBdata.DateTime)
c_s_avg_monthly=SHF_avg_monthly/(U10m_avg_monthly*(T2m_avg_monthly-Tsurface_avg_monthly))
c_s_mean=np.mean(c_s_avg_monthly)

""" plt.plot(np.linspace(0,1,14),c_s_avg_monthly)
plt.hlines(c_s_mean,0,1,colors='r',linestyles='dashed')
plt.xlabel("Month")
plt.ylabel("c_s")
plt.title("c_s for "+station)
plt.show() """

""" SHF_monthly=SEBf.get_monthly_average(SHF, SEBdata.DateTime)
T2m_monthly=SEBf.get_monthly_average(T2m, SEBdata.DateTime)
Tsurface_monthly=SEBf.get_monthly_average(Tsurface, SEBdata.DateTime)
U10m_monthly=SEBf.get_monthly_average(U10m, SEBdata.DateTime)
c_s_monthly=SHF_monthly/(U10m_monthly*(T2m_monthly-Tsurface_monthly))
plt.plot(np.linspace(0,1,c_s_monthly.shape[0]),c_s_monthly) """


#time,t,t2m,t_uncorr,q,q2m,rh,rh2m,ff,ff10m,ffmax,ffdir,p,
#SWd,SWd_notiltcorr,SWd_uncorr,SWu,SWu_uncorr,LWd,LWd_uncorr,LWu,LWu_uncorr,
#Tsub1a,Tsub2a,Tsub3a,Tsub4a,Tsub5a,Tsub1b,Tsub2b,Tsub3b,Tsub4b,Tsub5b,
#z_surf,z_surf_filtered,cum_surface_height_zboom,cum_surface_height,z_u,Vbat,TCNR,Tlogger,
#LAT,LON,alb,sza,Ts_obs,
#SHFdown,LHFdown,GHFup,cloud_cover,LWu_mod,SHFdown_mod,LHFdown_mod,GHFup_mod,Ts_mod,
#meltE,melt,cumulative_melt,sublimation,accumulation,
#LWC,LWC_10cm,z_boom,z_stakes,z_adw,z_adw_geus,z_pt_cor,dens,
#SHFdown_obs,LHFdown_obs,SHFdown_c,ff_EC,ffdir_EC,ustar_EC,Tsonic,Tcouple,z_EC


#%%
# Example plots of the SEB
# OPTIONS HERE ARE
# Typical yearly cycle -> AvgMonth
# Monthly averages     -> Monthly
# Daily averages       -> Daily

print("some example plots")

PlotType = "Monthly"

if PlotType == "AvgMonth":
    # Example of plotting the typical yearly cycle
    MyFunc = SEBf.get_avg_monthly_value
    Range  = [0, 12.5]
    Label  = "Month"
    Xdata  = np.arange(14)
elif PlotType == "Monthly":
    # Example of montly averages
    MyFunc = SEBf.get_monthly_average
    Range  = [SEBdata.DateTime[0], SEBdata.DateTime[-1]]
    Label  = "Year"
    # do dummy call to function to get the date-times for the x-axis
    Xdata  = MyFunc(SHF, SEBdata.DateTime, GiveDatesBack=True)[1]
elif PlotType == "Daily":
    # Example of montly averages
    MyFunc = SEBf.get_daily_average
    # restrain ourselves to one year, I take here 2016
    Range  = [dt.datetime.fromisoformat("2012-01-01"), dt.datetime.fromisoformat("2013-01-01")]
    Label  = "Date"
    # do dummy call to function to get the date-times for the x-axis
    Xdata  = MyFunc(SHF, SEBdata.DateTime, GiveDatesBack=True)[1]



# however, it is easier to do the call to the SEBf-function in the plot call...


#fig, axs = plt.subplots(2, sharex=True)
#fig.suptitle(station)

# upper figure are the radiative fluxes
#axs[0].plot(Xdata, MyFunc(SWdown, SEBdata.DateTime), 'b', linewidth=1.5, label="$SW_{down}$")
#axs[0].plot(Xdata, MyFunc(-SWup, SEBdata.DateTime), 'b:', linewidth=1.5, label="$SW_{up}$")
#axs[0].plot(Xdata, MyFunc(LWdown, SEBdata.DateTime), 'r', linewidth=0.5, label="$LW_{down}$")
#axs[0].plot(Xdata, MyFunc(-LWup, SEBdata.DateTime), 'r:', linewidth=0.5, label="$LW_{up}$")
#axs[0].plot(Xdata, MyFunc(SWnet+LWnet, SEBdata.DateTime), 'k', linewidth=0.5, label="$R_{net}$")
#axs[0].set_ylabel("Energy flux (W/m2)")
#axs[0].legend(loc='lower left') # well, no spot is nice
#axs[0].grid(True)

""" axs[1].plot(Xdata, MyFunc(-cloud_cover, SEBdata.DateTime), 'purple', linewidth=1.5, label="Cloud cover")
axs[1].plot(Xdata, MyFunc(Ts_mod, SEBdata.DateTime), 'r', linewidth=0.5, label="Ts_mod")
axs[1].legend(loc='lower left') # well, no spot is nice
axs[1].grid(True)
 """

""" axs[1].plot(Xdata, MyFunc(SWnet, SEBdata.DateTime), 'b', linewidth=0.5, label="$SW_{net}$")
axs[1].plot(Xdata, MyFunc(LWnet, SEBdata.DateTime), 'r', linewidth=0.5, label="$LW_{net}$")
axs[1].plot(Xdata, MyFunc(SHF, SEBdata.DateTime), 'seagreen', linewidth=0.5, label="$SHF$")
axs[1].plot(Xdata, MyFunc(LHF, SEBdata.DateTime), 'orange', linewidth=0.5, label="$LHF$")
axs[1].plot(Xdata, MyFunc(Gs, SEBdata.DateTime), 'grey', linewidth=0.5, label="$Gs$")
axs[1].plot(Xdata, MyFunc(MeltS, SEBdata.DateTime), 'purple', linewidth=0.5, label="$M$")
axs[1].plot(Xdata, MyFunc(Residu, SEBdata.DateTime), 'k', linewidth=0.5, label="Residue SEB model")
axs[1].set_ylabel("Energy flux (W/m2)")
axs[1].legend(loc='lower left') # Again, no spot is nice
axs[1].set_xlabel(Label)
axs[1].set_xlim(Range)
axs[1].grid(True) """

# Example, a simple adjustment: Increase LW down by making the apparent atmospheric temperature 1 K warmer

dtemp = 1.
Tsurf = SEBdata.Extract_Variable("Ts_mod")
LWmod = SEBf.convert_T_in_LWout(Tsurf) # calculates LWout from modelled Ts
LWout = LWnet - LWdown # original LWout, with correct sign convention, thus minus for outgoing, same as date LWup
Tatm   = SEBf.convert_LWout_in_T(-LWdown, Celcius=False) # calculate Ta from observed Ldown, -sign necessary to compensate for sign in the function
LWinDT = -SEBf.convert_T_in_LWout(Tatm+dtemp, Celcius=False) # calculate new LWin based on increased air temperature, including correct sign
epsilon = SEBf.get_epsilon_from_LWin_Temp(LWdown,T2m)
LWoutDTuc = LWnet - LWinDT # new LWout asuming LWnet remains the same, with correct sign convention
TskinDT = SEBf.convert_LWout_in_T(LWoutDTuc) # calculate new surface temp from new LWout
TskinDT = np.where(TskinDT>0., 0., TskinDT) # in case now Tskin > 0 set to 0degC
LWoutDT = SEBf.convert_T_in_LWout(TskinDT) # recalculate new LWout based on new surface temp
MeltTDT = MeltT + LWoutDT - LWoutDTuc # recalculate corrected amount of melt 
fig2, ax2 = plt.subplots()
ax2.plot(SEBdata.DateTime, SEBf.get_running_melt_sum(MeltTDT, SEBdata.TimeStep,ResetAtNan=False), 'r', label='Melt for +1K')
ax2.plot(SEBdata.DateTime, SEBf.get_running_melt_sum(MeltT,   SEBdata.TimeStep,ResetAtNan=False), 'g', label='Observed Melt')
ax2.legend(loc='upper left')
ax2.set_ylabel('Accumulated melt (m w.e.)')
ax2.set_xlabel('Year')
ax2.set_xlim([SEBdata.DateTime[0], SEBdata.DateTime[-1]])
ax2.set_ylim([0, 19])
ax2.text(0.85, 0.85, station, transform = ax2.transAxes) 

############################################################################################################

#test

def calc_melt(Tsurf,SW_down, SW_up, LW_down,Gs, Tsurf_old,SHF_old, LHF_old,c_s,U10m,sigma=5.67e-8):
    LW_up=sigma*Tsurf**4
    SHF=SHF_old+c_s*U10m*(Tsurf_old-Tsurf)
    LHF=LHF_old #maybe change later
    melt = SW_down - SW_up + LW_down - LW_up + SHF + LHF + Gs
    return melt


#bisection method for root finding
def find_root(func,x1,x2,args,tol=1e-10,N_max=100):
    #check if the initial values are correct
    if func(x1,*args)*func(x2,*args)>0:
        print("Root not found: Initial values do not bracket the root")
        return None
    #iterate
    i=0
    while i<N_max:
        x=(x1+x2)/2
        if func(x,*args)==0 or (x2-x1)/2<tol:
            print("Root found after",i,"iterations")
            return x
        if np.sign(func(x,*args))==np.sign(func(x1,*args)):
            x1=x
        if np.sign(func(x,*args))==np.sign(func(x2,*args)):
            x2=x
        i+=1
    print("Root not found: Maximum number of iterations reached")
    return None

    
def close_equation(T_init_min,T_init_max,args):
    new_temp=find_root(calc_melt,T_init_min,T_init_max,args)
    if new_temp is None:
        print("No solution found")
        return None, None
    if new_temp>0:
        #print("Temperature is above zero degrees Celsius")
        new_temp=0
        melt=calc_melt(new_temp,*args)
        print('T = 0°C')
        print(f'M = {melt:.3f} W/m^2')
    else:
        melt=0
        print(f'T = {new_temp:.3f}°C')
        print('M = 0 W/m^2')
    return new_temp, melt
    
#old values, paramaters needed, insert values here
func=SEBf.get_monthly_average
SW_down_old=func(SWdown, SEBdata.DateTime)
SW_up_old=func(SWup, SEBdata.DateTime)
LW_down_old=func(LWdown, SEBdata.DateTime)
Gs_old=func(Gs, SEBdata.DateTime)
Tsurf_old=func(Ts_mod, SEBdata.DateTime)
SHF_old=func(SHF, SEBdata.DateTime)
LHF_old=func(LHF, SEBdata.DateTime)
U10m = func(ff10m, SEBdata.DateTime)
T2m = func(T2m, SEBdata.DateTime)
c_s=SHF_old/(U10m*(T2m-Tsurf_old))


#sensitivity analysis
SW_down=SW_down_old*1 #adapt this
SW_up=SW_up_old*1 #adapt this
LW_down=LW_down_old*1 #adapt this
Gs=Gs_old #adapt this

#determine new surface temperature and melt
T_init_min=-100.
T_init_max=100.
N=U10m.shape[0]

stats_newTemp=np.zeros(N)
stats_melt=np.zeros(N)
for i in range(0,N):
    args=[SW_down[i], SW_up[i], LW_down[i],Gs[i], Tsurf_old[i],SHF_old[i], LHF_old[i],c_s[i],U10m[i]]
    stats_newTemp[i], stats_melt[i]=close_equation(T_init_min,T_init_max,args)

fig3, ax3 = plt.subplots(2)
ax3[0].plot(SEBdata.DateTime[:N], stats_newTemp, 'r', label='New Surface Temperature')
ax3[0].plot(SEBdata.DateTime[:N], Tsurf[:N], 'b', label='Original Surface Temperature')
ax3[0].legend(loc='upper left')
ax3[1].plot(SEBdata.DateTime[:N], stats_melt, 'g', label='Melt')
ax3[1].plot(SEBdata.DateTime[:N], MeltS[:N], 'b', label='Original Melt')
ax3[1].legend(loc='upper left')
plt.show()


