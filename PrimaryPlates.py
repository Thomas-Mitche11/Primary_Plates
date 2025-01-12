# These are the imports needed for the code to work

import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

from astropy.coordinates import (SkyCoord, Distance, Galactic, 
                                 EarthLocation, AltAz)
import astropy.coordinates as coord

from astropy import units as u
import os

import SciServer.CasJobs as CasJobs
import SciServer.SkyServer as SkyServer   # show individual objects and generate thumbnail images through SkyServer

import ipywidgets as widgets
from IPython.display import display



df = pd.read_csv('APOGEE_Plates.txt')
df1 = pd.read_csv('MaNGA_Plates.txt') #Reads in the plates numbers that I grabbed from HTML

# plate = 8603 #default plate

def get_plate(df):
    #This function cleans up the HTML and gives the Plate and MJD for a given data file
    Plate = []
    MJD = []
    for i in np.array(df):
        Plate.append(i[0].split(">", 1)[1].split("<")[0].split("/")[0]) #
        MJD.append(i[0].split(">", 1)[1].split("<")[0].split("/")[1])

    Plate = np.array(Plate, dtype=int)
    MJD = np.array(MJD, dtype=int)
    return Plate, MJD

def APOGEE_MaNGA_Plate(df, df1):
    #This function outputs an array of plates that are both MaNGA and APOGEE
    A_Plate, A_MJD = get_plate(df)
    M_Plate, M_MJD = get_plate(df1)
    
    AM_Ar = []
    for i in M_Plate:
        P, M = A_Plate[np.where(A_Plate == i)[0]], A_MJD[np.where(A_Plate == i)[0]]
        if len(P) > 0:
            for j in np.arange(len(P)):
                AM_Ar.append(np.array([P[j],M[j]]))
                
    return np.array(AM_Ar).T

def set_plate(No_plate = 'No', Plate_num_val = 8603, Plates = APOGEE_MaNGA_Plate(df, df1)[0]):
    global plate
    if No_plate == 'No':
        plate = np.random.choice(Plates)
        print('Plate number '+str(plate)+' was selected.')
        print('Move along to the next code cell')
    if No_plate == 'Yes':
        a = Plate_num_val
        if not a.isnumeric():
            print('')
            print('You need put in a number. No spaces or other characters')
            print('')
            return False
        if a.isnumeric():
            a = int(a)
            if a not in Plates:
                if a in get_plate(df)[0]:
                    print('')
                    print('This is just an APOGEE plate')
                    print('You will not be able to do the entirety of this activity')
                    print('Try again or select NO to recieve a random plate') 
                    print('')
                    return False

                else:
                    if a in get_plate(df1)[0]:
                        print('')
                        print('Just MaNGA plate')
                        print('You will not be able to do the entirety of this activity')
                        print('Try again or select NO to recieve a random plate')
                        print('')
                        return False

                    else:
                        print('')
                        print('This plate does not appear in our database.')
                        print('Try again or select NO to recieve a random plate') 
                        print('')
                        return False

            else:
                plate = a

    return plate            

def clean_stars(stars):
    
    stars = stars[(abs(stars['ra'] - stars['racen']) < stars['radius'])]
    
    stars_nr = stars.drop_duplicates('apogee_id')
    stars_nr

    dupindex_ar = []
    for aid in stars_nr['apogee_id']:
    #     print(aid)
        if all(stars[(stars['apogee_id'] == aid)]['age'] < 0):
            mindex = pd.Series.idxmin(stars[(stars['apogee_id'] == aid)]['mass_err'])
            dupindex = stars[(stars['apogee_id'] == aid)].index[stars[(stars['apogee_id'] == aid)].index != mindex]
        else:
            mindex = pd.Series.idxmin(stars[(stars['apogee_id'] == aid)&(stars['age'] > 0)]['mass_err'])
            dupindex = stars[(stars['apogee_id'] == aid)].index[stars[(stars['apogee_id'] == aid)].index != mindex]
        dupindex_ar.append(dupindex)
    c_stars = stars
    for i in dupindex_ar:
        c_stars = c_stars.drop(i)
#         c_stars = c_stars.reset_index()
    return c_stars

def Assign_Plate(Plates):
    AM_Button = widgets.widgets.ToggleButtons(
        options=['Yes', 'No'],
        value=None,
    #     description=' Do you have an APOGEE/MaNGA Plate? ',
        disabled=False,
        button_style='info', # 'success', 'info', 'warning', 'danger' or ''
        tooltips=['Enter the plate number in the cell below', 
                  'Just carry on with the activity and a random plate will be selected for you'],
    #     icons=['check'] * 3
    )
    output1 = widgets.Output()
    print('Do you have an APOGEE/MaNGA Plate?')
    display(AM_Button, output1)
    
    Plate_num = widgets.Combobox(
        # value='John',
        placeholder='APOGEE or MaNGA plate #',
        options= list(np.unique(Plates).astype(str)),
        description='Plate #:',
        ensure_option=False,
        disabled=False
    )

    select_button = widgets.Button(
        description='Set Plate',
        disabled=False,
        button_style='info', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click to set the plate',
    #     icon='check' # (FontAwesome names without the `fa-` prefix)
    )
    
    # tr = 0
    plate = 'test'
    def AM_Button_Click(b):
        global plate
        with output1:
    #         print('test')
    #         print(AM_Button.value)
            if AM_Button.value == 'No':
                plate = set_plate(AM_Button.value, Plate_num.value, Plates)
            if AM_Button.value == 'Yes':
                display(Plate_num)
                select_button.button_style = 'info'
                display(select_button)
                def push(g):
                    plate = set_plate(AM_Button.value, Plate_num.value, Plates)
                    if plate:
                        select_button.button_style = 'success'
                        print('Plate number '+str(plate)+' was selected.')
                        print('Move along to the next code cell')
                select_button.on_click(push)


    AM_Button.observe(AM_Button_Click, 'value')

def plate_set(p):
    global plate
    plate = p
    return plate
                
# plate = Plates[0]
def APOGEE_q(plate, query = False):
    if query:
        return query
    else:
        query_s = """
    select distinct star.ra, star.dec,
      217.7358*(star.ra - plate.racen)*cos(plate.deccen*PI()/180.0) as x,
    -217.7358*(star.dec - plate.deccen) as y,
    obj.j as mag_j, obj.k as mag_k,
    1./(star.gaiaedr3_parallax/1000.) as dist,
    plate.racen, plate.deccen, plate.plate, plate.radius, gaiaedr3_pmra as pmra,
    distmass.mass as mass, distmass.mass_err, distmass.age, star.apogee_id
    from apogeestar as star
      --join apogeeObject as obj on obj.target_id=star.target_id
      join apogeeObject as obj on obj.apogee_id=star.apogee_id
      join apogeeDistMass as distmass on distmass.apogee_id=star.apogee_id
      join apogeePlate as plate on plate.location_id=star.location_id
    where
      plate.plate = """+str(plate)+'\norder by mag_j\n'
        return(query_s)

def MaNGA_q(plate, query = False):
    if query:
        return query
    else:
        query_g = """
    SELECT distinct drp.mangaid as id, drp.ifudsgn, drp.ifura, drp.ifudec, 
    drp.nsa_elpetro_mass, drp.plate,
    nsa_zdist*300000./67. as dist,
    217.7358*(drp.ifura - drp.cenra)*cos(drp.cendec*PI()/180.0) as x,
    -217.7358*(drp.ifudec - drp.cendec) as y,
    22.5 - 2.5*LOG(drp.nsa_elpetro_flux_r, 10) as mag_r, drp.nsa_elpetro_absmag_r as absmag_r,
    22.5 - 2.5*LOG(drp.nsa_elpetro_flux_g, 10) as mag_g, drp.nsa_elpetro_absmag_g as absmag_g,
    22.5 - 2.5*LOG(drp.nsa_elpetro_flux_u, 10) as mag_u, drp.nsa_elpetro_absmag_u as absmag_u
    FROM mangadrpall as drp
    WHERE
    drp.plate ="""+str(plate)+'\norder by mag_r\n'
        return query_g



def Get_Star_Gal(plate, query = False):
    query_test = True
    while query_test == True:
        try:
            unfiltered_stars = CasJobs.executeQuery(APOGEE_q(plate, query = query), "dr17")
            if query:
                stars = unfiltered_stars
            else:
                stars = clean_stars(unfiltered_stars)
                stars = stars.reset_index()
                stars = stars.drop(columns='index')
            
            gals = CasJobs.executeQuery(MaNGA_q(plate, query = query), "dr17")
            query_test = False
        except Exception:
            print('There appears to be and error with recieving data with plate '+str(plate))
            print('a new plate will be found for you')
            plate = plate_set(set_plate('No'))
            
    return stars, gals

# blue = '#307efb'
# green = '#5acc77'

def star_sky_plot(stars, gals, blue = '#307efb', green = '#5acc77', flip = True):


    #make a circle
    xc = np.arange(-400, 401)  # don't understand why this is not 400
    yc = np.sqrt(400**2-xc**2)

    #make a tab
    sztab = 50.
    xtab = np.array([-sztab, -sztab, sztab, sztab])
    ytabval = np.sqrt(400**2-xtab**2)
    ytab = ytabval + np.array([0, 40, 40, 0])


    def plate_plot(x,y, blue = '#307efb', green = '#5acc77', gal_tf = False):
        plt.figure(figsize=(15,15))

        # plt.plot(x, y, marker='$o$', ls='None', c='k')
        plt.plot(x, y, 'o', c=blue)
        if gal_tf:
            plt.plot(gals['x'], gals['y'], 'o', c=green)
        plt.plot(0,0,'ko') 

        plt.plot(xc, yc, c='k')
        plt.plot(xc, -1.*yc, c='k')


        plt.plot(xtab, ytab, c='k')



        ax = plt.gca()
        ax.set_aspect(aspect=1)
        ax.text(xtab[0]+20, ytab[0]+10, plate, fontsize=17)
        #ax.get_xaxis().set_visible(False)
        #ax.get_yaxis().set_visible(False)
        plt.axis('off')  # also removes border


        plt.show()
    #------------------------------------------------------------------------------------------------------------------
    # set thumbnail parameters
    width=900           # image width
    height= width          # height
    pixelsize=1    # image scale

    scale=2.785*(stars['radius'][0] * u.deg).to(u.arcsec).value/pixelsize/width #Needed to scale the image (works with 8625)
    # scale=2*(stars['radius'][0] * u.arcsec).value/pixelsize/width
    img = SkyServer.getJpegImgCutout(ra=stars['racen'][0], dec=stars['deccen'][0], width=width, height=height,
                                     scale=scale)


    #---------------------------------------------------------------------------------------------------------------


    sky = widgets.ToggleButtons(
        options=['Plate','Sky','Both'],
    #     description=' Do you have an APOGEE/MaNGA Plate? ',
        disabled=False,
        button_style='info', # 'success', 'info', 'warning', 'danger' or ''
        tooltips=['A plot of the APOGEE/MaNGA plate that matches the number you submitted', 
                  'A look at the night sky in the area of the plate',
                 'A side-by side of the plate and the sky'],
    #     icons=['check'] * 3
    )


    num_holes = widgets.ToggleButtons(
        options=np.sort(np.array([0,15,30,45, len(stars['x'])//2, len(stars['x'])])),
        value = 45,
        disabled=False,
        button_style='info', # 'success', 'info', 'warning', 'danger' or ''
    #     icons=['check'] * 3
    )
    output2 = widgets.Output()


    #----------------------------------------------------------------------------------------------------------------



    def gerb(b):
        with output2:
            Gal_widget.update()

    # def f(Green_O):
    #     Green_O = num_holes.value
    #     print(Green_O)
    #     return Green_O

    def stars_show_plot(gv, holes, gal_tf, blue ='#307efb', green ='#5acc77', flip = True):
        if flip:
            imflip = 2
            rot = 1
        else:
            imflip = 0
            rot = -1
        
        if gv == 'Plate':
            plate_plot(stars['x'][:holes],stars['y'][:holes], gal_tf = gal_tf, blue = blue, green = green)

        if gv == 'Sky':
            plt.figure(figsize=(15, 15)) 
            plt.imshow(np.rot90(img, imflip))

            plt.plot(rot*stars['x'][:holes]+width/2, width/2-rot*stars['y'][:holes], marker='o',fillstyle='none', ls='None', c=blue, markersize=10)
            if gal_tf:
                plt.plot(rot*gals['x']+width/2, -rot*gals['y']+width/2, marker='o',fillstyle='none', ls='None', c=green, markersize=10)
            plt.plot(width/2, width/2, 'ro', fillstyle='none', markersize=700)
            plt.axis('off')
            plt.show()

        if gv == 'Both':
            fig, (ax1, ax2) = plt.subplots(ncols=2,figsize=(30,15))
            ax1.plot(stars['x'][:holes], stars['y'][:holes], 'o', c=blue)
            if gal_tf:
                ax1.plot(gals['x'], gals['y'], 'o', c=green)
                ax2.plot(rot*gals['x']+width/2, -rot*gals['y']+width/2, marker='o',fillstyle='none', ls='None',c=green, markersize=10)
            ax1.plot(0,0,'ko') 
            ax1.plot(xc, yc, c='k')
            ax1.plot(xc, -1.*yc, c='k')
            ax1.plot(xtab, ytab, c='k')
            ax1.text(xtab[0]+20, ytab[0]+10, plate, fontsize=17)
            ax1.axis('off')  # also removes border

            ax2.imshow(np.rot90(img, imflip))
            ax2.plot(rot*stars['x'][:holes]+width/2, width/2-rot*stars['y'][:holes], marker='o',fillstyle='none', ls='None', c=blue, markersize=10)
            ax2.plot(width/2, width/2, 'ro', fillstyle='none', markersize=700)
            ax2.axis('off')

            plt.show()

    def Display_plots(Galaxies, c = ['#307efb','#5acc77'], f = True): 
    #     print(Green_O)
        blue = c[0]
        green = c[1]
        gv = sky.value
        holes = num_holes.value
        if Galaxies == 'No':
            stars_show_plot(gv, holes, gal_tf=False, blue = blue, green = green, flip = f)
        if Galaxies == 'Yes':
            stars_show_plot(gv, holes, gal_tf=True, blue = blue, green = green, flip = f)



    #--------------------------------------------------------------------------------------------------------
    print('Pick the number of holes to be displayed:')
    display(num_holes, output2)

    print('Pick the background:')
    display(sky, output2)

    sky.observe(gerb, 'value')
    num_holes.observe(gerb, 'value')


    Gal_widget = widgets.interactive(Display_plots, Galaxies=['Yes','No'], c=widgets.fixed([blue, green]), f = widgets.fixed(flip));
    display(Gal_widget)
    
    
    
       
def comp_IMG(stars, gals, opt = '',hexa = False):
#     global prior_star_num, prior_gal_num
    comp_img = widgets.ToggleButtons(
        options=['Blue Circle','Green Circle','Both'],
        disabled=False,
        button_style='info', # 'success', 'info', 'warning', 'danger' or ''
        tooltips=['Displays one of the blue objects', 
                  'Displays one of the green objects',
                 'Displayes the objects side by side'],

    )

    prior_star_num = 1
    prior_gal_num = 1


 
    def comp_star_plot(ax, img, s, opt):
        right = 0
        if opt:
            right = 240
        ax.text(20+right, 25+right/24, 'Star '+stars['apogee_id'][s], fontsize=17, c='white')
        ax.text(90+right, 40+right/24, '(APOGEE ID)', fontsize=10, c='white')
        ax.imshow(img)
        ax.axis('off')
        
    def comp_gal_plot(ax, img, g, opt):
        right = 0
        if opt:
            right = 330
        ax.text(20+right, 25+right/33, 'Galaxy '+gals['id'][g], fontsize=17, c='white')
        ax.text(90+right, 40+right/33, '(MaNGA ID)', fontsize=10, c='white')
        ax.imshow(img)
        ax.axis('off')
        
    def comp_img_plot(green_circle,blue_circle, opt, hexa):

        s = blue_circle-1
        g = green_circle-1
        
        scale = 0.4
        width = 512
        
        if comp_img.value == 'Blue Circle':

            fig, ax = plt.subplots(figsize=(10, 10))

            img_s = SkyServer.getJpegImgCutout(ra=stars['ra'][s],
                                               dec=stars['dec'][s], width=width, height=width, scale=scale, opt=opt)
            comp_star_plot(ax, img_s, s, opt)
            plt.show()

        if comp_img.value == 'Green Circle':
            
            IFU = np.round(gals['ifudsgn'].to_numpy(dtype=int)/100)
            scale = (IFU[g] + 20)/width
            
            fig, ax = plt.subplots(figsize=(10, 10))

            img_g = SkyServer.getJpegImgCutout(ra=gals['ifura'][g],
                                               dec=gals['ifudec'][g], width=width, height=width, scale=scale, opt=opt)
            comp_gal_plot(ax, img_g, g, opt)
            if hexa:
                ax.add_patch(patches.RegularPolygon([width/2,width/2], 6, radius = IFU[g]/scale/2, orientation=np.pi/2, ec = 'white', fill=False))
            plt.show()

        if comp_img.value == 'Both':
            fig, ax = plt.subplots(ncols=2,figsize=(20,10))

            img_s = SkyServer.getJpegImgCutout(ra=stars['ra'][s],
                                               dec=stars['dec'][s], width=width, height=width, scale=scale, opt=opt)
            
            
            IFU = np.round(gals['ifudsgn'].to_numpy(dtype=int)/100)
            scale = (IFU[g] + 20)/width
            img_g = SkyServer.getJpegImgCutout(ra=gals['ifura'][g],
                                               dec=gals['ifudec'][g], width=width, height=width, scale=scale, opt=opt)

            comp_star_plot(ax[0], img_s, s, opt)
            comp_gal_plot(ax[1], img_g, g, opt)
            if hexa:
                ax2.add_patch(patches.RegularPolygon([width/2,width/2], 6, radius = IFU[g]/scale/2, orientation=np.pi/2, ec = 'white', fill=False))
            plt.show()

            





    def gerb2(b):
        with output3:
            interactive_plot_comb.update()

    output3 = widgets.Output()
    interactive_plot_comb = widgets.interactive(comp_img_plot, green_circle=np.arange(len(gals))+1, blue_circle=np.arange(2*len(gals))+1, opt=widgets.fixed(opt), hexa=widgets.fixed(hexa))


    print('Pick Which Images are Displayed:')
    display(comp_img, output3)
    display(interactive_plot_comb)
    # display(interactive_plot_comb)
    comp_img.observe(gerb2,'value')
    
    
    

def dist_analysis(stars, gals, color = ['#307efb','#5acc77']):
    dist_star = stars['dist'][(stars['dist'] > 0)]
    dist_star_xlim = np.linspace(min(dist_star),round(max(dist_star)*1.1), 10).astype(int)[1:]

    dist_star_xlim2 = np.linspace(min(dist_star_xlim),round(max(dist_star)*100), 5).astype(int)[1:-1]

    dis_gal = gals['dist'][(gals['dist'] > 0)]*1e6
    dist_bet_xlim = np.linspace(max(dist_star_xlim2), min(dis_gal),5).astype(int)[1:-1]


    dis_gal_xlim = np.linspace(min(dis_gal), round(max(dis_gal)*1.1), 10).astype(int)

    event_xlim = np.concatenate([dist_star_xlim, dist_star_xlim2, dist_bet_xlim, dis_gal_xlim], axis= 0)

    blue, green = color

    def disevent(x):
        fig, ax = plt.subplots(figsize=(20, 5))
    #     ax.eventplot(stars['dist'], label = 'Blue Circles (stars)', color=blue)
    #     ax.eventplot(gals['dist']*1e6, color=green, label = 'Green Circles (galaxies)')

        ax.plot(stars['dist'], np.ones_like(stars['dist']), label = 'Blue Circles (stars)', color=blue, marker='o', ls='', markeredgecolor='k', markersize=15)
        ax.plot(gals['dist']*1e6, np.ones_like(gals['dist']*1e6), color=green, label = 'Green Circles (galaxies)', marker='o', ls='', markeredgecolor='k', markersize=25)

        ax.set_xlim(0,x)
        ax.legend(fontsize=20)
        ax.set_xlabel('Distance in pc', fontsize=30)
        if x > 1e6:
            ax.set_xlabel('Distance in Millions of pc', fontsize=30)
        if x > 1e7:
            ax.set_xlabel('Distance in 10 Millions of pc', fontsize=30)
        if x > 1e8:
            ax.set_xlabel('Distance in 100 Millions of pc', fontsize=30)
        ax.get_yaxis().set_visible(False)
        ax.xaxis.set_tick_params(labelsize=20)
    #     ax.tick_params(axis='both', which='minor', labelsize=20)
        plt.show()

    max_scroll = widgets.interactive(disevent, x=widgets.widgets.SelectionSlider(
        options= event_xlim,
        value=event_xlim[5],
        description='Max Distance',
        disabled=False,
        continuous_update=False,
        orientation='horizontal',
        readout=True
    ));

    display(max_scroll)

def mass_analysis(stars, gals, color = ['#307efb','#5acc77'], sort = ['dist', 'dist'], comb=True):
    blue, green = color
    def mass_plot(x, comb = True):

        stars_sorted = stars.sort_values(sort[0])
        gals_sorted = gals.sort_values(sort[1])

        fig, ax = plt.subplots(figsize=(15, 10))
        ax.bar(np.arange(len(stars_sorted[(stars_sorted['mass'] > 0)])),
               stars_sorted['mass'][(stars_sorted['mass'] > 0)], color = blue,
              label = 'Blue Circles (stars)')
        if comb:
            ax.bar(len(stars_sorted[(stars_sorted['mass'] > 0)]),
                   sum(stars_sorted['mass'][(stars_sorted['mass'] > 0)]), color = 'yellow',
                  label = 'Every Star on the Plate Combined', width=5)
            ax.bar(np.arange(len(stars_sorted[(stars_sorted['mass'] > 0)])+3, 
                             len(stars_sorted[(stars_sorted['mass'] > 0)])+3+len(gals_sorted['id'])),
                   gals_sorted['nsa_elpetro_mass'], color = green,
                  label = 'Green Circles (galaxies)')
        else:
            ax.bar(np.arange(len(stars_sorted[(stars_sorted['mass'] > 0)]), 
                             len(stars_sorted[(stars_sorted['mass'] > 0)])+len(gals_sorted['id'])),
                   gals_sorted['nsa_elpetro_mass'], color = green,
                  label = 'Green Circles (galaxies)')
        ax.set_ylabel('Mass in Suns', fontsize=30)
        if x > 1e12:
            ax.set_ylabel('Mass in Trillions of Suns', fontsize=30)
        elif x > 1e11:
            ax.set_ylabel('Mass in 100 Billions of Suns', fontsize=30)
        elif x > 1e10:
            ax.set_ylabel('Mass in 10 Billions of Suns', fontsize=30)
        elif x > 1e9:
            ax.set_ylabel('Mass in Billions of Suns', fontsize=30)
        elif x > 1e8:
            ax.set_ylabel('Mass in 100 Millions of Suns', fontsize=30)
        elif x > 1e7:
            ax.set_ylabel('Mass in 10 Millions of Suns', fontsize=30)
        elif x > 1e6:
            ax.set_ylabel('Mass in Millions of Suns', fontsize=30)
        elif x > 1e5:
            ax.set_ylabel('Mass in 100 Thousands of Suns', fontsize=30)
        ax.set_title('Mass of Objects (sorted by '+sort[0]+')', fontsize=30)
        ax.legend(fontsize=20)
        ax.get_xaxis().set_visible(False)
        ax.set_ylim(0,x)
        plt.show()


    mass_scroll = widgets.interactive(mass_plot, x=widgets.widgets.SelectionSlider(
            options= np.round([max(stars['mass']), 10, 100, 1e3, 1e5, min(gals['nsa_elpetro_mass'])*1.5, np.mean(gals['nsa_elpetro_mass']), max(gals['nsa_elpetro_mass'])],2),
    #         value=event_xlim[5],
            description='Max Mass',
            disabled=False,
            continuous_update=False,
            orientation='horizontal',
            readout=True), comb = widgets.fixed(comb));


    display(mass_scroll)

def zero_point():
    zero_points = {
    "U": [1790], # Unused for SDSS
    "B": [4063], # Unused for SDSS
    "V": [3636], # Unused for SDSS
    "R": [3064], # Unused for SDSS
    "I": [2416], # Unused for SDSS
    "J": [1589], # 2MASS Filters (Vega Mag)
    "H": [1021], # 2MASS Filters (Vega Mag)
    "Ks":[640],  # 2MASS Filters (Vega Mag)
    "u": [3631], # sdss filter (AB)
    "g": [3631], # sdss filter (AB)
    "r": [3631], # sdss filter (AB)
    "i": [3631], # sdss filter (AB)
    "z": [3631]  # sdss filter (AB)
    }

    zp = pd.DataFrame(zero_points)
    return(zp)
zp = zero_point()

def Sun_Abs_Mag():
    Abs_mag = { 
    'U': [5.61, 6.33],
    'B': [5.44, 5.31],
    'V': [4.81, 4.80],
    'R': [4.43, 4.60],
    'I': [4.10, 4.51],
    'J': [3.67, 4.54],
    'H': [3.32, 4.66],
    'Ks': [3.27, 5.08],
    'u': [5.49, 6.39],
    'g': [5.23, 5.11],
    'r': [4.53, 4.65],
    'i': [4.19, 4.53],
    'z': [4.01, 4.50]
    }

    AMag = pd.DataFrame(Abs_mag, index = ['Vega','AB'])
    return AMag
AMag = Sun_Abs_Mag()

def mag_flux(M1, Ffilter):
    zp = zero_point()
    F = zp[Ffilter][0]*10**((M1)/(-2.5))
    F.name = 'F-'+Ffilter
    return F

def app_abs(m1, dist):
    #dist is in PC
    M1 = m1 - 5*np.log10(dist) + 5
    return M1

def abs_mag_Lum(M1, Ffilter, Mag = None):
    M_sun = 4.83
    AMag = Sun_Abs_Mag()
    if Mag:
        M2 = AMag[Ffilter].loc[Mag]
    elif Ffilter in ['U','B','V','R','I','J','Ks']:
        M2 = AMag[Ffilter].loc['Vega']
    elif Ffilter in ['u','g','r','i','z']:
        M2 = AMag[Ffilter].loc['AB']

    L2_Lsun = 10**((M_sun - M2)/(2.5))

    L1_L2 = 10**((M1 - M2)/(-2.5))

    return L1_L2*L2_Lsun

def mag_flux(M1, Ffilter):
    zp = zero_point()
    F = zp[Ffilter][0]*10**((M1)/(-2.5))
    F.name = 'F-'+Ffilter
    return F

def app_abs(m1, dist):
    #dist is in PC
    M1 = m1 - 5*np.log10(dist) + 5
    return M1

def abs_mag_Lum(M1, Ffilter, Mag = None):
    M_sun = 4.83
    AMag = Sun_Abs_Mag()
    if Mag:
        M2 = AMag[Ffilter].loc[Mag]
    elif Ffilter in ['U','B','V','R','I','J','Ks']:
        M2 = AMag[Ffilter].loc['Vega']
    elif Ffilter in ['u','g','r','i','z']:
        M2 = AMag[Ffilter].loc['AB']

    L2_Lsun = 10**((M_sun - M2)/(2.5))

    L1_L2 = 10**((M1 - M2)/(-2.5))

    return L1_L2*L2_Lsun

def plot_lum_flux(stars, starsfilter, gals, galsfilter, color = ['#307efb', '#5acc77'], plot='luminosity'):
    blue, green = color
    fig=plt.figure(figsize = (15, 10))
    ax1 = plt.subplot(122)
    ax2 = plt.subplot(121)
    
    if plot == 'flux':
        yg = mag_flux(gals[galsfilter[0]], galsfilter[1])
        ys = mag_flux(stars[starsfilter[0]], starsfilter[1])
        ax2.set_ylabel('Flux (Jy)', fontsize=16)
    elif plot == 'luminosity':
        yg = abs_mag_Lum(gals[galsfilter[0]], galsfilter[1])
        ys = abs_mag_Lum(app_abs(stars[starsfilter[0]],stars['dist']), starsfilter[1])
        ax2.set_ylabel('Luminosity (L_sun)', fontsize=16)
    
    
    
    ax1.scatter(1e6*gals['dist'], yg, c=green, label='Green Circles (galaxies)')
    
    ax1.legend(fontsize=15)
    ax1.set_xlim(0.5e8,3e8)
    ax1.set_xlabel('Distance in 100 Million Pc', fontsize=16)
    ax1.get_shared_y_axes().join(ax1, ax2)
    ax1.set_yticklabels([])

    ax2.scatter(stars['dist'], ys, c=blue, label = 'Blue Circles (stars)')
    
    
    
    ax2.set_xlabel('Distance in Pc', fontsize=16)
    ax2.set_xlim(0,10e3)
    
#     ax2.set_ylim(0,max([max(yg),max(ys)]))
#     print(max([max(yg),max(ys)]))
    ax2.legend(fontsize=15)

    plt.show()
    
def star_gal_analysis(stars, gals, color = ['#307efb', '#5acc77'], sort = ['dist', 'dist'], comb=True):
    def analysis(x):
        if x == 'Distance':
            dist_analysis(stars, gals, color = color)
        if x == 'Mass':
            mass_analysis(stars, gals, color = color, sort = sort, comb=comb)
        if x == 'Brightness':
            fluxlum = widgets.ToggleButtons(
                    options=['Flux','Luminosity'],
                    description=' ',
                    disabled=False,
                    value= 'Flux',
                    button_style='info', # 'success', 'info', 'warning', 'danger' or ''
                    tooltips=['Flux is the amount of light that we can see here on Earth', 
                              'Luminosity is the actual light given off by the object'],

                )
            def plotLF(x):
                if x == 'Flux':
                    plot_lum_flux(stars, ['mag_k','Ks'], gals, ['mag_u', 'u'], plot='flux')
                elif x == 'Luminosity':
                    plot_lum_flux(stars, ['mag_k','Ks'], gals, ['absmag_u', 'u'], plot='luminosity')

            fluxlum_interactive = widgets.interactive(plotLF, x=fluxlum, 
                                                       stars=widgets.fixed(stars),
                                                       gals=widgets.fixed(gals), 
                                                      color=widgets.fixed(color));

            display(fluxlum_interactive)

    analysis_buttons = widgets.ToggleButtons(
        options=['Distance','Mass','Brightness'],
        description='Analysis',
        disabled=False,
        value= None,
        button_style='info', # 'success', 'info', 'warning', 'danger' or ''
        tooltips=['Compare the distances between stars and galaxies', 
                  'Compare the masses of stars and galaxies',
                 'Compare the brightness of stars and galaxies'],

    )

    analysis_interactive = widgets.interactive(analysis, x=analysis_buttons, 
                                               stars=widgets.fixed(stars),
                                               gals=widgets.fixed(gals), 
                                              color=widgets.fixed(color));

    display(analysis_interactive)
