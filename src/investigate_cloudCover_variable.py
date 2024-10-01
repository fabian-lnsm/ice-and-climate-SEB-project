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
c_s=SHF/(ff10m*(T2m-Ts_mod)+1e-6) #correction factor to avoid division by zero

#change data to daily/monthly/avgmonthly. Also, create the according time data for plotting later
timestep='daily' # 'daily', 'monthly', 'avgmonthly'
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
LWup=f(LWup,SEBdata.DateTime)
ff10m=f(ff10m,SEBdata.DateTime)
T2m=f(T2m,SEBdata.DateTime)
SHF=f(SHF,SEBdata.DateTime)
LHF=f(LHF,SEBdata.DateTime)
Gs=f(Gs,SEBdata.DateTime)
MeltS=f(MeltS,SEBdata.DateTime)
Residu=f(Residu,SEBdata.DateTime)
c_s=f(c_s,SEBdata.DateTime)
Ts_mod=f(Ts_mod,SEBdata.DateTime)
cloud_cover=f(cloud_cover,SEBdata.DateTime)


############################################################################

#remove nan values
indices_to_remove = ~np.isnan(cloud_cover)
Time_data=Time_data[indices_to_remove]
cloud_cover=cloud_cover[indices_to_remove]
SWdown=SWdown[indices_to_remove]
SWup=SWup[indices_to_remove]
LWdown=LWdown[indices_to_remove]
LWup=LWup[indices_to_remove]
SHF=SHF[indices_to_remove]
LHF=LHF[indices_to_remove]
Gs=Gs[indices_to_remove]
Ts_mod=Ts_mod[indices_to_remove]


#print some values and correlations
print(f'Correlation: Cloud cover & SWdown: {np.corrcoef(cloud_cover,SWdown)[0,1]:.2f}')
print(f'Correlation: Cloud cover & SWup: {np.corrcoef(cloud_cover,SWup)[0,1]:.2f}')
print(f'Correlation: Cloud cover & LWdown: {np.corrcoef(cloud_cover,LWdown)[0,1]:.2f}')
print(f'Correlation: Cloud cover & LWup: {np.corrcoef(cloud_cover,LWup)[0,1]:.2f}')
print(f'Correlation: Cloud cover & SHF: {np.corrcoef(cloud_cover,SHF)[0,1]:.2f}')
print(f'Correlation: Cloud cover & LHF: {np.corrcoef(cloud_cover,LHF)[0,1]:.2f}')
print(f'Correlation: Cloud cover & Gs: {np.corrcoef(cloud_cover,Gs)[0,1]:.2f}')


###########################################################################

fig, axs = plt.subplots(2,sharex=True)
fig.suptitle(f'Station: {station}, timestep:{timestep} ')
axs[0].set_ylabel('')
axs[0].plot(Time_data,cloud_cover,label='cloud_cover')
axs[0].legend()
axs[0].grid(True)
#axs[1].set_ylabel('')
#axs[1].legend()
#axs[1].grid(True)
#axs[1].set_xlabel(Time_Label)
#plt.show()

