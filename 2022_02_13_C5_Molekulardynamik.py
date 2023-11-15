# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 11:01:20 2022

@author: sonne
"""
#0. Imports

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.patches as patches
import matplotlib.animation as animation
import matplotlib.ticker as ticker
import tkinter as Tk #Interface
import numpy as np #Numerische Mathematik
import itertools #für Iteration
import random #Zufallszahlen
from scipy.constants import k #Boltzmann-Konstante

#1. Dictionaries für die Edelgase

Masse = {"Helium" : 6.642e-27, 
         "Neon" : 3.351e-26, 
         "Argon" : 6.634e-26, 
         "Krypton" : 1.392e-25, 
         "Xenon": 2.180e-25}

Durchmesser = {"Helium" : 36.58,  # 1.4e-10m
         "Neon" : 44.22,       # 1.58e-10m
         "Argon" : 65.96,      # 1.88e-10m
         "Krypton" : 76.15,    # 2.00e-10m
         "Xenon": 87.07}       # 2.18e-10m

Farbe = {"Helium" : "blue", 
         "Neon" : "darkblue", 
         "Argon" : "blueviolet", 
         "Krypton" : "purple", 
         "Xenon": "indigo"}

LennardJones_epsilon = {"Helium" : 14e-23, 
         "Neon" : 50e-23, 
         "Argon" : 167e-23, 
         "Krypton" : 225e-23, 
         "Xenon": 320e-23}

LennardJones_sigma = {"Helium" : 2.56e-10, 
         "Neon" : 2.74e-10, 
         "Argon" : 3.4e-10, 
         "Krypton" : 3.65e-10, 
         "Xenon": 3.98e-10}

#2. Figure für Plot erstellen

#2.1 für die Simulation

fig = Figure(figsize = (6,6))
ax = fig.add_subplot(111)

#2.2 für Beschleunigungspfeil

arrow = Figure(figsize = (1.5,1.5)) #Zweites Plotfenster für den Beschleunigungsvektor
arr = arrow.add_subplot(111)
arr.set_xlim(-1,1)
arr.set_ylim(-1,1)
arr.axes.xaxis.set_visible(False) #Achsen zur Übersichtlichkeit ausgeblendet
arr.axes.yaxis.set_visible(False)
Inaktiv_text = arr.text(0,0,'Aktiviere \n"Teilchen verfolgen"', ha="center", va = "center", fontsize = 8.) #Text in der Figure zu Beginn
Pfeil = patches.Arrow(0, 0, 0, 0) #Pfeilobjekt mit Länge 0 erstellen
patch = arr.add_patch(Pfeil) #Pfeil zu Plot hinzufügen

#3. Interface mit Tkinter

class Interface():
    
    def __init__(self,  Teilchentracker_aus = True):
        self.root = Tk.Toplevel()
        self.Teilchentracker_aus = Teilchentracker_aus
   
    #3.1 Tkinter-Fenster konfigurieren
        
        self.root.title("C5 Molecular Dynamics") #Titel des Fensters
        self.root.geometry("1400x800") #Größe des Fensters in Pixel
        self.root.config(bg = "white") #weißer Hintergrund
        self.root.columnconfigure(0, weight=3)
        self.root.columnconfigure(1, weight=1)
        self.root.columnconfigure(2, weight =1)
        
        self.Ueberschrift = Tk.Label(self.root,text="Thermal motion", font = "Verdana 20 bold", \
                                     bg = "white").grid(column=0, row=0)
    
        #3.2 Canvas für die Simulation und für den Beschleunigungspfeil erstellen
        
        self.canvas = FigureCanvasTkAgg(fig, master=self.root) #für Simulation
        self.canvas.get_tk_widget().grid(column=0, row=1, rowspan = 9, sticky = Tk.N)
        
        self.Label_Beschleunigungspfeil = Tk.Label(self.root, text = "Acceleration",\
                                                   font = "Verdana 10 bold", bg = "white").grid(column = 2, row = 1)
        self.canvas_arrow = FigureCanvasTkAgg(arrow, master=self.root) #für Pfeil
        self.canvas_arrow.get_tk_widget().grid(column=2, row =2, rowspan = 2, sticky = Tk.N, pady = 10)

        #3.2 Schieberegler für Änderung der Temperatur
    
        self.Label_Temperatur = Tk.Label(self.root, text = "Temperature in K", font = "Verdana 10 bold",\
                                         bg = "white").grid(column = 1, row =1) 
        self.Slider_Temperatur = Tk.Scale(self.root, from_=1, to=2000, orient = "horizontal",\
                                          bg = "white") #Schieberegler
        self.Slider_Temperatur.grid(column = 1, row = 2) #Schieberegler platzieren
        self.Slider_Temperatur.set(300) #Startwert
        self.Button_Temperatur = Tk.Button(self.root, text="Change temperature", bg= "lightgreen", \
                                           compound = "left", width = 18, command= \
                                           self.update_Temperatur).grid(column = 1, row = 3) 
                                           #Knopf, ruft Funktion für Temperatur auf
                                           
    
        #3.3 Schieberegler für Änderung der Teilchenzahl
            
        self.Label_Teilchenzahl = Tk.Label(self.root, text = "Number of particles", \
                                           font = "Verdana 10 bold", bg = "white").grid(column = 1, row = 4, sticky= Tk.S)
        self.Slider_Teilchenzahl = Tk.Scale(self.root, from_=1, to=20,\
                                            orient = "horizontal", bg = "white") #Schieberegler
        self.Slider_Teilchenzahl.grid(column = 1, row = 5) #Schieberegler platzieren
        self.Slider_Teilchenzahl.set(5) #Startwert
        self.Button_Teilchenzahl = Tk.Button(self.root, text="Change number of particles",\
                                             bg = "lightgreen", compound = "left", width = 18,\
                                             command=self.update_Teilchenzahl).grid(column = 1, row = 6) 
                                             #Knopf, ruft Funktion für Teilchenzahl auf

        #3.4 Dropdownmenü für Änderung der Teilchenart
        
        self.Label_Teilchenart = Tk.Label(self.root, text = "Gas type", font = "Verdana 10 bold",\
                                          bg = "white").grid(column = 1, row = 7, sticky = Tk.S)
        Edelgase = ["Helium","Neon", "Argon","Krypton","Xenon"] #Liste der Optionen
        Variable = Tk.StringVar() #Definition des Wert des Widgets, Hält eine Zeichenfolge; Standardwert ""
        Variable.set(Edelgase[0]) #gibt an welches Element der Liste im Menü angezeigt wird
        self.dropdown = Tk.OptionMenu(self.root, Variable, *Edelgase, command= \
                                      self.update_Teilchenart).grid(column = 1, row = 8)  #Widget für Dropdown-Menü erstellen
        
        #3.5 Label mit Informationen zur aktuellen Simulation
        
        self.Infos = Tk.Label(self.root, text = "Informationen", bg = "white", font = \
                              "Verdana 10 bold").grid(column = 2, row = 7, sticky = Tk.S)
        self.Label_Infos = Tk.Label(self.root, text = "Infos", justify = "left") #Label erstellen
        self.Label_Infos.grid(column = 2, row = 8) #Label platzieren

        #3.6 Teilchentracker zum An- und Ausschalten

        self.Label_Teilchentracker = Tk.Label(self.root, text = "Track particle", font = \
                                              "Verdana 10 bold", bg = "white").grid(column = 2, row = 4, sticky = Tk.S)
        self.Button_teilchen_verfolgen = Tk.Button(self.root, fg = "darkgreen", text="Teilchen verfolgen", bg = "white",\
                                                   height = 2, command=self.teilchen_verfolgen).grid(column = 2, row = 5, \
                                                   rowspan = 2, sticky = Tk.N, pady = 12)
                                                   #Knopf, aktiviert Teilchentracker  
                        

        #3.7 Stopp-Knopf zum Beenden des Programms
        
        self.Beenden = Tk.Button(self.root, text = "Interrupt", fg= "white", bg="maroon",\
                                 command = self.stopp, width = 65).grid(column = 1, row = 9, columnspan = 2)
                                 #Knopf, ruft Stopp-Funktion auf
    
    #Funktionen für Tkinter-Schaltflächen:
        
    def update_Temperatur(self):        
        box.Temperatur = self.Slider_Temperatur.get() #Wert des Schiebereglers abrufen
        box.start_Animation() #Startbedingungen aktualisieren
        
    def update_Teilchenzahl(self):            
        box.Teilchenzahl = self.Slider_Teilchenzahl.get() #Wert des Schiebereglers abrufen
        box.particles = [Particle(i) for i in range(box.Teilchenzahl)] #neue Teilchen erstellen
        box.start_Animation() #Startbedingungen aktualisieren
    
    def update_Teilchenart(self, Variable): 
        partikel.Teilchenart = str(Variable) #Teilchenart als String speichern (für Infolabel)
        partikel.m = Masse.get(Variable) #Masse im Dictionary nachschlagen und aktualisieren
        partikel.R = Durchmesser.get(Variable) #Teilchenradius im Dictionary nachschlagen und aktualisieren
        partikel.color = Farbe.get(Variable) #Farbe im Dictionary nachschlagen und aktualisieren
        partikel.epsilon = LennardJones_epsilon.get(Variable) #Parameter im Dictionary nachschlagen und aktualisieren
        partikel.sigma = LennardJones_sigma.get(Variable) #Parameter im Dictionary nachschlagen und aktualisieren
        box.start_Animation() #Startbedingungen aktualisieren
        
    def teilchen_verfolgen(self): #Möglichkeit die Farbe eines Teilchens zu ändern und den Beschleunigungsvektor zu verfolgen
        global patch, Inaktiv_text
        if self.Teilchentracker_aus: #wenn noch nicht aktiv
            Inaktiv_text.remove() #Text aus Plot entfernen
            arrow.canvas.draw()
            self.Button_teilchen_verfolgen = Tk.Button(self.root, foreground = "white",\
                                                       text="Teilchen entfolgen", bg = "darkgreen", height = 2,\
                                                       command=self.teilchen_verfolgen).grid(column = 2, row = 5,\
                                                       rowspan = 2, sticky = Tk.N, pady = 12)
                                                       #Knopf ändert sein Aussehen                                     
            self.Teilchentracker_aus = False 
            
        else: #falls schon aktiv
            Inaktiv_text = arr.text(0,0,'Aktiviere \n"Teilchen verfolgen"', ha="center", va = "center", fontsize = 8.)
                #Text in Plot einfügen
            arrow.canvas.draw()
            self.Button_teilchen_verfolgen = Tk.Button(self.root, foreground = "darkgreen",\
                                                       text="Teilchen verfolgen", bg = "white", height = 2,\
                                                       command=self.teilchen_verfolgen).grid(column = 2, row = 5,\
                                                       rowspan = 2, sticky = Tk.N, pady = 12)
                                                       #Knopf ändert sein Aussehen                                     
            patch.remove() #Pfeil entfernen
            Pfeil = patches.Arrow(0, 0, 0, 0) #einen neuen Pfeil der Länge 0 erstellen
            patch = arr.add_patch(Pfeil) #Pfeil hinzufügen
            arrow.canvas.draw() #Pfeil anzeigen
            self.Teilchentracker_aus = True 
            
    def stopp(self):
       self.root.destroy() #Tkinter-Fenster schließen
       self.root.quit() #Programmausführung stoppen
         
interface = Interface() #auf die Klasse Interface() mit interface zugreifen  
    
#4. Beschleunigungspfeil für das erste Teilchen

def pfeil(Beschleunigung):
    global patch
    if interface.Teilchentracker_aus == False: #nur wenn Teilchentracker aktiv ist
        Betrag_Beschleunigung = np.sqrt(Beschleunigung[0]**2 + Beschleunigung[1]**2)
        if Betrag_Beschleunigung != 0: 
            patch.remove() #Pfeil entfernen
            Pfeil_x = Beschleunigung[0]/np.abs(Beschleunigung[0]) \
                * np.log(np.abs(Beschleunigung[0]))/50 #logarithmisch skalierte Beschleunigung, 
                                                        #um die nötigen Größenordnungen abzudecken
            Pfeil_y = Beschleunigung[1]/np.abs(Beschleunigung[1]) \
                * np.log(np.abs(Beschleunigung[1]))/50 #logarithmisch skalierte Beschleunigung
            Pfeil = patches.FancyArrow(0, 0, Pfeil_x, Pfeil_y, color = "maroon", width = 0.05, overhang = 0.2,\
                                       head_width = 0.25, head_length = 0.3) #Pfeil mit den Komponenten der 
                                                                                #Beschleunigung erstellen
            patch = arr.add_patch(Pfeil) #Pfeil hinzufügen
            arrow.canvas.draw() #Pfeil anzeigen
   
    
#5. Kritischer Radius

kritischerRadius = 4e-9 # entspricht 40% der Boxgröße, um Rechenaufwand zu reduzieren

#6. Teilchen als Klasse:

class Particle(): #Ordnet jedem Teilchen einen Radius (R), Masse (m), Farbe (color),
                    #und die Parameter für das Lennard-Jones-Potential (sigma, epsilon) zu
   def __init__(self, R = 36.58, m = 6.642e-27, color = "blue", epsilon = 14e-23, sigma = 2.56e-10, Teilchenart = "Helium"):
       self.R, self.m, self.color, self.epsilon, self.sigma, self.Teilchenart = R, m, color, epsilon, sigma, Teilchenart
    
partikel = Particle() #auf Klasse Particle() als partikel zugreifen

# 7. Funktionen für die Bewegung der Teilchen in der Box

class Box(): #enthält die Funktionen für die Bewegung der Teilchen in der Box

    def __init__(self, Teilchenzahl = 5, dt=4E-15, Temperatur = 300, Boxgroesse = 1e-8,\
                 Anfangsgeschwindigkeit = 1367.8, E_gesamt = 3.1e-20): 
                 #Default-Werte für Anzahl der Teilchen, Zeitintervall, Temperatur, Boxgröße, Anfangsgeschwindigkeit, Gesamtenergie
        self.dt, self.Teilchenzahl, self.Temperatur, self.Boxgroesse, self.Anfangsgeschwindigkeit, \
            self.E_gesamt  = dt, Teilchenzahl, Temperatur, Boxgroesse, Anfangsgeschwindigkeit, E_gesamt
        self.particles = [Particle(i) for i in range(self.Teilchenzahl)] #für jedes Teilchen eine Instanz der Particle-Klasse erstellen
    
    #7.1 Startbedingungen für die Simulation berechnen und festlegen
                  
    def start_Animation(self):
        self.scatter = ax.scatter([],[], s= partikel.R) #Streudiagramm mit Teilchen als Marker, Größe entsprechend des Teilchenradius
        self.Anfangsgeschwindigkeit = self.mittlere_Geschwindigkeit(self.Temperatur) #Anfangsgeschwindigkeit aus Temperatur berechnen
        self.E_gesamt = self.gesamtenergie(self.Teilchenzahl, self.Anfangsgeschwindigkeit) #Gesamtenergie aus kinetischer Energie bestimmen
        Infos = "Edelgas: " + partikel.Teilchenart + \
            "\nMasse: %10.3e kg \nGesamtenergie: %10.3e J \nMittlere Geschwindigkeit: %8.2f m/s" \
            % (partikel.m, box.E_gesamt, box.Anfangsgeschwindigkeit) #Text für das Info-Label
        interface.Label_Infos.configure(text = Infos) #Labelinhalt aktualisieren
        box.startpositionen() #Startpositionen-Funktion aufrufen             
        for particle in self.particles:
            angle = random.uniform(-1, 1) #zufälliger Winkel für Richtung der Geschwindigkeit
            particle.v = np.array([(np.cos(angle * np.pi/2)), (np.sin(angle * np.pi/2))]) \
                * self.Anfangsgeschwindigkeit #Anfangsgeschwindigkeit als Array definieren
            particle.a = np.zeros(2) #Beschleunigung zum Zeitpunkt t=0 ist Null
    
    #7.1.1 Anfangsgeschwindigkeit der Teilchen als mittlere Geschwindigkeit festlegen
    
    def mittlere_Geschwindigkeit(self, T):
        return np.sqrt(3*k*T/partikel.m) #Als mittlere Geschwindigkeit über Temperatur berechnen
     
    #7.1.2 Gesamtenergie aller Teilchen berechnen
    
    def gesamtenergie(self, Teilchenzahl, v):
        return Teilchenzahl * 0.5 * partikel.m * v**2 #Summe der kinetischen Energie
    
    #7.1.3 Startpostitionen der Teilchen zufällig bestimmen
    
    def startpositionen(self):
        for particle in self.particles:
            particle.r = 100*np.random.uniform(0, self.Boxgroesse/100, size=2) #Startposition zufällig 
                                                                                #innerhalb des Kastens festlegen
       
        #Wiederholung der zufälligen Teilchenverteilung bei Überlappung von 2 Teilchen zu Beginn der Animation
        
        for particle, particle2 in itertools.combinations(self.particles, 2): #für jedes Teilchenpaar
              # Abstand berechnen
              x_diff = particle.r[0] - particle2.r[0]
              y_diff = particle.r[1] - particle2.r[1]
              Abstand = np.sqrt(x_diff**2 + y_diff**2)
              
              if Abstand < 1.12*partikel.sigma: #wenn Abstand kleiner abstoßende Wechselwirkunegn
                  box.startpositionen() #neue Startpositionen berechnen
      
    #7.2 Trajektorien der Teilchen über Velocity-Verlet-Algorithmus bestimmen                 

    def zeitliche_Entwicklung(self, particles, Boxgroesse, dt, E_gesamt): 
        box.kollision_Box(particles, Boxgroesse) #elastische Stöße mit Wand berücksichtigen
        for particle in particles:
            particle.r += dt * particle.v + dt**2*particle.a #Ort nach Velocity-Verlet-Algorithmus bestimmen
            particle.a_vorher = particle.a #Wert für die Beschleunigung für nächsten Zeitschritt speichern
            particle.a = np.zeros(2) #Beschleunigung wieder auf Null setzen vor neuer Evaluation des Potentials
            box.beschleunigung(particles) #Beschleunigung aus dem Potential berechnen
            particle.v = (particle.v + dt/2 * (particle.a + particle.a_vorher)) \
                * box.normierung(particles, E_gesamt) #Geschwindigkeit nach Velocity-Verlet-Algorithmus bestimmen und normieren
        pfeil(box.particles[0].a)  #Beschleunigungspfeil updaten
                   
    #7.2.1 Elastische Stöße mit den Wänden der Box             
                    
    def kollision_Box(self, particles, Boxgroesse): 
       for particle in self.particles:
           for i in range(0,2): #für x- und y-Koordinate
               if particle.r[i] >= Boxgroesse: 
                       particle.r[i] = Boxgroesse #Verhindert 'Tunneling', wo die Geschwind. eines sehr schnellen 
                                             # Teilchens vor Rückkehr in den Kasten zweimal gespiegelt wird
                       particle.v[i] *=-1 #Spiegelung der Geschwindigkeit
               if particle.r[i] <= 0:
                       particle.r[i] = 0 #s.o.
                       particle.v[i] *=-1   
                       
    #7.2.2 Abstand und Beschleunigung der Teilchen bestimmen
    
    def beschleunigung(self, particles):
       for particle, particle2 in itertools.combinations(self.particles, 2):  #über alle Paare von Teilchen iterieren
                    
           #Abstand berechnen
           x_diff = particle.r[0] - particle2.r[0]
           y_diff = particle.r[1] - particle2.r[1]
           Abstand = np.sqrt(x_diff**2 + y_diff**2)
           
           #Wechselwirkung aus Potential berechnen:
          
           if Abstand < kritischerRadius: #nur Wechselwirkung bestimmen, wenn innerhalb des kritischen Radius
                Wechselwirkung = self.lennardJones_Kraft(Abstand) #Abstand in Lennard-Jones-Potential einsetzen
                particle.a[0] -= 1/(partikel.m) * Wechselwirkung * x_diff/Abstand
                particle.a[1] -= 1/(partikel.m) * Wechselwirkung * y_diff/Abstand
                particle2.a[0] += 1/(partikel.m) * Wechselwirkung * x_diff/Abstand
                particle2.a[1] += 1/(partikel.m) * Wechselwirkung * y_diff/Abstand
   
    #7.2.3 Lennard-Jones-Potential
    
    def lennardJones_Kraft(self, Distanz): #Kraft als Gradient des LennardJones-Potentials in Abhängigkeit 
                                            #vom Abstand der Teilchen
        return (-24 * partikel.epsilon) * (2 *(partikel.sigma**12 / Distanz**13) -(partikel.sigma**6 / Distanz**7))         
    
    #7.2.4 Geschwindigkeiten für Energieerhaltung normieren
    
    def normierung(self, particles, E_gesamt): 
        Summe_v=0 
        for particle in particles: 
            Summe_v += particle.v**2 #alle Geschwindigkeitsquadrate aufaddieren
        return np.sqrt(E_gesamt /(0.5*Particle().m*Summe_v)) #neue Gesamtenergie bestimmen und Skalierungsfaktor zurückgeben
    
    #7.3 Position der Teilchen zurückgeben
    
    def position(self, particles): #Funktion um den Ort jedes Teilchens an die Animation zu übergeben
        return [particle.r for particle in particles]     

   
box = Box() #auf die Klasse Box() als box zugreifen


#8. Animation starten und Programm ausführen       

        
def particle_Farbe(particles):
    for particle in box.particles:
     particle.color = partikel.color #jedem Teilchen die Farbe des Edelgases aus der Particle-Klasse zuweisen
    if interface.Teilchentracker_aus == False: #wenn Teiclhentracker aktiviert
        box.particles[0].color = "red" #erstes Teilchen wird rot eingefärbt
    return [particle.color for particle in box.particles] #Farben zurückgeben
   
def init(): #Box darstellen
    ax.set_xlim (0, box.Boxgroesse) #Boxgröße einstellen
    ax.set_ylim (0, box.Boxgroesse)
    ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(0,12e-9, 2e-9))) #Ticklabels einstellen
    ax.yaxis.set_major_locator(ticker.FixedLocator(np.arange(0,12e-9, 2e-9)))
    ax.xaxis.set_ticklabels(np.arange(0,11,2)) #Beschriftungen in nm festlegen
    ax.yaxis.set_ticklabels(np.arange(0,11,2))
    ax.set_xlabel("Boxbreite in nm") #Achsenbeschriftung
    return box.scatter,
 
def update(frame): #Funktion für die Animation (FunAn)
    box.zeitliche_Entwicklung(box.particles, box.Boxgroesse, box.dt, box.E_gesamt) #Funktion für zetilche Entwicklung aufrufen
    box.scatter.set_offsets(np.array(box.position(box.particles))) #neuen Ort der Marker übernehmen
    box.scatter.set_color(particle_Farbe(box.particles)) #Farbe der Teilchen ggf. ändern
    return box.scatter,
 
box.start_Animation() #Funktion für Startbedingungen aufrufen

ani = animation.FuncAnimation(fig, update , frames=range(10000), init_func = init, blit=True,\
                              interval = 1/2000, repeat = True) #Animation abrufen
Tk.mainloop() #Tkinter-Fenster aufrufen