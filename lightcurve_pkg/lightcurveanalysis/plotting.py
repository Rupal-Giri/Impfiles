import matplotlib.pyplot as plt
import pandas as pd

def band_sep_FG(data, savepath):
    """
    Plot flux vs MJD for different frequency bands.
    """
    data_2 = data[(data['Freq'] >= 2) & (data['Freq'] <= 3)]
    data_5 = data[(data['Freq'] >= 4) & (data['Freq'] <= 5)]
    data_8 = data[(data['Freq'] >= 7) & (data['Freq'] <= 9)]
    data_10 = data[(data['Freq'] >= 10) & (data['Freq'] <= 11)]
    data_14 = data[(data['Freq'] >= 14) & (data['Freq'] <= 15)]
    data_23 = data[(data['Freq'] >= 23) & (data['Freq'] <= 24)]
    data_32 = data[(data['Freq'] >= 31) & (data['Freq'] <= 33)]
    data_43 = data[(data['Freq'] >= 42) & (data['Freq'] <= 44)]

    plt.figure(figsize=(10, 6))
    plt.errorbar(data_2['MJD']-2400000.5, data_2['Flux'], yerr=data_2['Flux err'], fmt='o', label='2GHz',color='crimson')
    plt.errorbar(data_5['MJD']-2400000.5, data_5['Flux'], yerr=data_5['Flux err'], fmt='o', label='5GHz',color='orange')
    plt.errorbar(data_8['MJD']-2400000.5, data_8['Flux'], yerr=data_8['Flux err'], fmt='o', label='8GHz',color='gold')
    plt.errorbar(data_10['MJD']-2400000.5, data_10['Flux'], yerr=data_10['Flux err'], fmt='o', label='10GHz',color='forestgreen')
    plt.errorbar(data_14['MJD']-2400000.5, data_14['Flux'], yerr=data_14['Flux err'], fmt='o', label='14GHz',color='mediumturquoise')
    plt.errorbar(data_23['MJD']-2400000.5, data_23['Flux'], yerr=data_23['Flux err'], fmt='o', label='23GHz',color='dodgerblue')
    plt.errorbar(data_32['MJD']-2400000.5, data_32['Flux'], yerr=data_32['Flux err'], fmt='o', label='32GHz',color='rebeccapurple')
    plt.errorbar(data_43['MJD']-2400000.5, data_43['Flux'], yerr=data_43['Flux err'], fmt='o', label='43GHz',color='black')
    
    plt.plot(data_2['MJD']-2400000.5, data_2['Flux'], color='crimson')
    plt.plot(data_5['MJD']-2400000.5, data_5['Flux'], color='orange')
    plt.plot(data_8['MJD']-2400000.5, data_8['Flux'], color='gold')
    plt.plot(data_10['MJD']-2400000.5, data_10['Flux'], color='forestgreen')
    plt.plot(data_14['MJD']-2400000.5, data_14['Flux'], color='mediumturquoise')
    plt.plot(data_23['MJD']-2400000.5, data_23['Flux'], color='dodgerblue')
    plt.plot(data_32['MJD']-2400000.5, data_32['Flux'], color='rebeccapurple')
    plt.plot(data_43['MJD']-2400000.5, data_43['Flux'], color='black')
    
    plt.grid(True)
    plt.xlabel('MJD')
    plt.ylabel('Flux (Jy)')
    plt.title('Flux vs. MJD')
    plt.legend()
    plt.savefig(savepath, dpi=300)
    
    return plt.gcf()

def band_sep_med(data):
    data_5 =  data[(data['Freq (MHz)'] >= 4000) & (data['Freq (MHz)'] <= 7000)]
    data_8 =  data[(data['Freq (MHz)'] >= 7000) & (data['Freq (MHz)'] <= 9000)]
    data_20 =  data[(data['Freq (MHz)'] >= 15000) & (data['Freq (MHz)'] <= 30000)]

    plt.figure(figsize=(10, 6))
    plt.scatter(data_5['# MJD'], data_5['Flux (Jy)'], label='5GHz',color='orange')
    plt.scatter(data_8['# MJD'], data_8['Flux (Jy)'], label='8GHz',color='forestgreen')
    plt.scatter(data_20['# MJD'], data_20['Flux (Jy)'], label='20GHz',color='dodgerblue')
    plt.plot(data_5['# MJD'], data_5['Flux (Jy)'],color='orange')
    plt.plot(data_8['# MJD'], data_8['Flux (Jy)'],color='forestgreen')
    plt.plot(data_20['# MJD'], data_20['Flux (Jy)'],color='dodgerblue')
    plt.grid(True)
    plt.xlabel('MJD')
    plt.ylabel('Flux (Jy)')
    plt.title('Flux vs. MJD ')
    plt.legend()
    plt.savefig('Impfiles/images/flux_vs_mjd_.png',dpi=300)

    return plt.gcf()

def band_sep_FG_and_med(data_FG, data_med,savepath):
    #Define colors
    colors_med=['orange','forestgreen','dodgerblue']
    colors_FG=['crimson','orange','gold','forestgreen','mediumturquoise','dodgerblue','rebeccapurple','black']

    # Define frequency bands (GHz)
    freq_bands = {
        '2GHz': (1, 3), '5GHz': (4, 6), '8GHz': (7, 9), '10GHz': (9, 11),
        '14GHz': (13, 15), '23GHz': (22, 25), '32GHz': (31, 33), '43GHz': (42, 44)
    }

    plt.figure(figsize=(10, 6))

    # Plot data from band_sep_FG (source_data)
    for label, (low, high) in freq_bands.items():
        band_data_FG = data_FG[(data_FG['Freq'] >= low) & (data_FG['Freq'] <= high)]
        if not band_data_FG.empty:
            plt.errorbar(
                band_data_FG['MJD'] - 2400000.5, band_data_FG['Flux'], yerr=band_data_FG['Flux err'], 
                fmt='o', label=f"{label} (FG)", alpha=0.5, capsize=3,color=colors_FG[list(freq_bands.keys()).index(label)]
            )

    # Convert MHz to GHz in data_med if necessary
    if 'Freq (MHz)' in data_med.columns:
        data_med['Freq'] = data_med['Freq (MHz)'] / 1000  # Convert MHz to GHz

    # Plot data from band_sep_med (f_data)
    for label, (low, high) in freq_bands.items():
        band_data_med = data_med[(data_med['Freq'] >= low) & (data_med['Freq'] <= high)]
        if not band_data_med.empty:
            plt.errorbar(
                band_data_med['# MJD'], band_data_med['Flux (Jy)'], yerr=band_data_med['e_flux (Jy)'], 
                fmt='x', label=f"{label} (Med)", alpha=0.7, capsize=3,color=colors_FG[list(freq_bands.keys()).index(label)]
            )

    # Plot formatting
    plt.grid(True)
    plt.xlabel('MJD')
    plt.ylabel('Flux (Jy)')
    plt.title('Flux vs. MJD (FG & Med Data with Errors)')
    plt.legend()
    plt.savefig(savepath,dpi=300)
    plt.show()
   
    return plt.gcf()

def band_sep_med_bach(f_data, source_data_bach, filepath):

    # Filtering data for med
    f_data_5 = f_data[(f_data['Freq (MHz)'] >= 4000) & (f_data['Freq (MHz)'] <= 7000)]
    f_data_8 = f_data[(f_data['Freq (MHz)'] >= 7000) & (f_data['Freq (MHz)'] <= 9000)]
    f_data_20 = f_data[(f_data['Freq (MHz)'] >= 15000) & (f_data['Freq (MHz)'] <= 30000)]

    # Filtering data for bach
    f_datab_5 = source_data_bach[(source_data_bach['Freq'] >= 4000) & (source_data_bach['Freq'] <= 7000)]
    f_datab_8 = source_data_bach[(source_data_bach['Freq'] >= 7000) & (source_data_bach['Freq'] <= 9000)]
    f_datab_20 = source_data_bach[(source_data_bach['Freq'] >= 15000) & (source_data_bach['Freq'] <= 30000)]
    f_datab_43 = source_data_bach[(source_data_bach['Freq'] >= 42000) & (source_data_bach['Freq'] <= 44000)]

    # Plot
    plt.figure(figsize=(12, 7))

    # 0235+164 (Circles)
    plt.errorbar(f_data_5['# MJD'], f_data_5['Flux (Jy)'], yerr=f_data_5['e_flux (Jy)'], fmt='o', 
                label='5GHz', color='orange', alpha=0.7)
    plt.errorbar(f_data_8['# MJD'], f_data_8['Flux (Jy)'], yerr=f_data_8['e_flux (Jy)'], fmt='o', 
                label='8GHz', color='forestgreen', alpha=0.7)
    plt.errorbar(f_data_20['# MJD'], f_data_20['Flux (Jy)'], yerr=f_data_20['e_flux (Jy)'], fmt='o', 
                label='20GHz', color='dodgerblue', alpha=0.7)
    plt.plot(f_data_5['# MJD'], f_data_5['Flux (Jy)'], color='orange')
    plt.plot(f_data_8['# MJD'], f_data_8['Flux (Jy)'], color='forestgreen')
    plt.plot(f_data_20['# MJD'], f_data_20['Flux (Jy)'], color='dodgerblue')
    # 3C66A (Squares)
    plt.errorbar(f_datab_5['MJD'], f_datab_5['Flux'], yerr=f_datab_5['Flux err'], fmt='s', 
                label='5GHz', color='orange', alpha=0.7)
    plt.errorbar(f_datab_8['MJD'], f_datab_8['Flux'], yerr=f_datab_8['Flux err'], fmt='s', 
                label='8GHz', color='forestgreen', alpha=0.7)
    plt.errorbar(f_datab_20['MJD'], f_datab_20['Flux'], yerr=f_datab_20['Flux err'], fmt='s', 
                label='20GHz', color='dodgerblue', alpha=0.7)
    plt.errorbar(f_datab_43['MJD'], f_datab_43['Flux'], yerr=f_datab_43['Flux err'], fmt='s', 
                label='43GHz', color='black', alpha=0.7)
    plt.plot(f_datab_5['MJD'], f_datab_5['Flux'], color='orange')
    plt.plot(f_datab_8['MJD'], f_datab_8['Flux'], color='forestgreen')
    plt.plot(f_datab_20['MJD'], f_datab_20['Flux'], color='dodgerblue')
    plt.plot(f_datab_43['MJD'], f_datab_43['Flux'], color='black')
    # Labels and Title
    plt.xlabel('MJD', fontsize=14)
    plt.ylabel('Flux (Jy)', fontsize=14)
    #plt.title('Flux vs. MJD for 0235+164 and 3C66A', fontsize=16)
    plt.legend(loc='best', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.savefig(filepath, dpi=300)
    # Show plot
    plt.show()
    return plt.gcf()

