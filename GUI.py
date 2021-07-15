import re
from tkinter import *
from tkinter import simpledialog
from tkinter import messagebox
import tkinter.scrolledtext as st
import tkinter as tk
from tkinter import filedialog as fd
import logging
import time
import pathlib

from drawing import *
from Biologic import *

SCRIPT_DIR = pathlib.Path(".")

## initiate logging
# start_time = time.strftime("%Y%m%d_%H%M%S")
start_time = time.strftime("%Y%m%d")
logging.basicConfig(
    filename=f"./log/{start_time}.log",
    encoding="utf-8",
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)s - %(name)s: %(message)s"
)
log = logging.getLogger()
# log.setLevel(logging.DEBUG)

# def _switch_gui_logging_level(menu:tk.Menu):
#     if logging.root.level > logging.DEBUG:
#         log.setLevel(logging.DEBUG)
#         log.info("Debug logging activated")
#         menu.entryconfig(1, label='Deactivate debug logging')
#     else:
#         log.setLevel(logging.INFO)
#         log.info("Debug logging deactivated")
#         menu.entryconfig(1, label='Activate debug logging')


# Dialog für die Bewertung einer Struktur
class eval_Dialog():

    def __init__(self, master, title):
        self.master = master
        self.dialog = Toplevel(master)
        self.dialog.geometry("400x120")
        self.dialog.title(title)
        # Überschrift
        w1 = Label(self.dialog, text = "Struktur (als Klammerausdruck):", font = ("Helvetica", 11))
        w1.grid(row = 0, column = 0)
        # Eingabefeld
        self.e1 = Entry(self.dialog, font = ("Arial", 11))
        self.e1.bind("<Return>", self.inStrukt)
        self.e1.grid(row = 0, column = 1)
        # Ausgabe-Label
        w2 = Label(self.dialog, text = "Anzahl Wasserstoffbruecken: ", font = ("Helvetica", 11))
        w2.grid(row = 1, column = 0, pady =4)
        # Ausgabefeld
        self.w3 = Label(self.dialog, foreground = "blue", font = ("Arial", 11))
        self.w3.grid(row = 1, column = 1)
        # Schaltknopf
        close_button = Button(self.dialog,text = "Close", comman = self.close, font = ("Helvetica", 11))
        close_button.grid(row= 2, pady = 5, padx = 5)

    # Methode zum Schließen des Dialogs
    def close(self):
        self.dialog.destroy()

    # Überprüfe die Struktur auf gültige Zeichen
    def checkStrukt(self,s,n):
        s = s.strip('\n')
        p = re.compile("[^(, ), .]")
        return p.search(s) == None and len(s) == n

    # Auswertung der Struktur
    def inStrukt(self, text):
        s = self.e1.get()
        if self.checkStrukt(s, len(self.master.logic.get_rna())):
            num = self.master.logic.eval_rna(s)
            self.w3.configure(text = str(num))
        else:
            tk.messagebox.showinfo(title = "Syntax error", message = "keine korrekte Klammerstruktur")
    

# Hauptfenster
class App(Tk):
    """ das Hauptfenster """
    def __init__(self):
        super().__init__()
        # Konfiguration
        self.title('HSWT RNA Structure')
        self.geometry('800x1000')


# class ScrollableFrame(Frame):
#     def __init__(self, container, *args, **kwargs):
#         super().__init__(container, *args, **kwargs)
#         canvas = Canvas(self)
#         scrollbar = Scrollbar(self, orient='vertical', command=canvas.yview)
#         self.scrollable_frame = Frame(canvas)

#         self.scrollable_frame.bind(
#             '<Configure>', 
#             lambda e: canvas.configure(
#                 scrollregion=canvas.bbox('all')
#             )
#         )

#         canvas.create_window((0,0), window=self.scrollable_frame, anchor='nw')

#         canvas.configure(yscrollcommand=scrollbar.set)

#         canvas.pack(side='left', fill='both', expand=True)
#         scrollbar.pack(side='right', fill='y')


# Der Rahmen des Hauptfensters
class SeqMainFrame(Frame):

    def __init__(self, container, application):
        # container ist das übergeordnete Fenster
        super().__init__(container)
        # die logic enthält die Algorithmen
        self.logic = application

        # DNA-Anzeige
        dna_label = Label(self, text = "DNA Sequenz:", anchor = 'w', font = ("Helvetica", 11))
        dna_label.grid(row = 0, column = 0, sticky = 'w',pady = 1, padx = 5)
        # DNA-Darstellung
        self.dna_text = st.ScrolledText(self, width = 60, height = 3, font = ("Arial", 11))
        self.dna_text.grid(column = 0, pady = 1, padx = 5)
        self.dna_text.insert(tk.INSERT,"")
        # die Sequenz ist hier nicht veränderbar 
        self.dna_text.configure(state ='disabled')

        # RNA-Anzeige
        rna_label = Label(self, text = "RNA Sequenz:", anchor = 'w', font = ("Helvetica", 11))
        rna_label.grid(row = 3, column = 0, sticky = 'w',pady = 1, padx = 5)
        # RNA-Darstellung
        self.rna_text = st.ScrolledText(self, width = 60, height = 3, font = ("Arial", 11), foreground = "blue")
        self.rna_text.grid(column = 0, pady = 1, padx = 5)
        self.rna_text.insert(tk.INSERT,"")
        # die RNA ist hier nicht veränderbar
        self.rna_text.configure(state ='disabled')
        
        # Unterrahmen für die Darstellung der optimalen Anzahl von Wasserstoffbrücken
        number_frame = tk.Frame(self)
        number_frame.grid(column = 0, pady = 5, padx = 5, sticky = 'w')

        anzahl_label = Label(number_frame, text = "Optimale Anzahl von Wasserstoffbrücken: ", anchor = 'w', font = ("Helvetica", 11))
        anzahl_label.pack(side = 'left')

        self.opt_anzahl = Label(number_frame, text = '0', foreground = "blue", font = ("Arial", 11))
        self.opt_anzahl.pack(side = 'left')

        # Unterrahmen for die Darstellung der optimalen Sekundärstruktur
        struct_frame = tk.Frame(self)
        struct_frame.grid(column = 0, pady = 5, padx = 5, sticky = 'w')

        struct_label = Label(struct_frame, text = "Optimale Struktur: ", anchor = 'w', font = ("Helvetica", 11))
        struct_label.pack(side = 'left')
        
        self.struct_text = st.ScrolledText(self, width = 60, height = 3, font = ("Arial", 11), foreground = "blue")
        self.struct_text.grid(column = 0, pady = 1, padx = 5)

        # Unterrahmen für die Zeichnung der optimalen Struktur
        canvas_frame=Frame(self,width=700,height=400)
        canvas_frame.grid(column = 0, pady = 5, padx = 5, sticky = 'w') #.grid(row=0,column=0)
        self.canvas = Canvas(canvas_frame, relief = RIDGE, bd = 2, bg = "white", width = 700, height = 400, scrollregion=(0,0,1200,1200))
        # Scrollbars rechts und unten
        hbar=Scrollbar(canvas_frame,orient=HORIZONTAL)
        hbar.pack(side=BOTTOM,fill=X)
        hbar.config(command=self.canvas.xview)
        vbar=Scrollbar(canvas_frame,orient=VERTICAL)
        vbar.pack(side=RIGHT,fill=Y)
        vbar.config(command=self.canvas.yview)
        self.canvas.config(width=700,height=400)
        self.canvas.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
        self.canvas.pack(side=LEFT,expand=True,fill=BOTH)

        # Darstellung der translatierten Aimosäurensequenz
        aas_label = Label(self, text = "Aminosäuren-Sequenz:", anchor = 'w', font = ("Helvetica", 11))
        aas_label.grid(column = 0, sticky = 'w',pady = 5, padx = 5)

        self.aas_text = st.ScrolledText(self, width = 40, height = 3, font = ("Arial", 11), foreground = "blue")
        self.aas_text.grid(column = 0, pady = 1, padx = 5)

        # Darstellung der gefundenen Proteine
        prot_label = Label(self, text = "Proteine:", anchor = 'w', font = ("Helvetica", 11))
        prot_label.grid(column = 0, sticky = 'w',pady = 5, padx = 5)

        self.prot_text = st.ScrolledText(self, width = 40, height = 3, font = ("Arial", 11), foreground = "blue")
        self.prot_text.grid(column = 0, pady = 1, padx = 5)

        # Erzeuge die Menuleiste
        menubar = tk.Menu(container)
        container.config(menu=menubar)

        # Erzeuge ein File-Menu mit Eintrag 'Exit'
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label='Exit',command=container.destroy)
        menubar.add_cascade(label="File",menu=file_menu)

        # Erzeuge ein DNA-Menu
        dna_menu = tk.Menu(menubar, tearoff=0)

        # Menu-Einträge
        dna_menu.add_command(label='Enter sequence',command=self.read_dna)
        dna_menu.add_command(label='Read fasta file',command=self.read_fasta)
        dna_menu.add_command(label='Retrieve from genbank',command=self.access_genbank)
        
        # füge das Menu in die Menubar ein
        menubar.add_cascade(label="DNA",menu=dna_menu)

        # Erzeuge ein RNA-Menu
        rna_menu = tk.Menu(menubar, tearoff=0)

        # Menu-Einträge
        rna_menu.add_command(label='Transcribe',command=self.transcribe)
        rna_menu.add_command(label='Evaluate Structure',command=self.eval_rna)
        rna_menu.add_command(label='Optimal Structure (Folding)',command=self.fold_rna)
         
        # füge das Menu in die Menubar ein
        menubar.add_cascade(label="RNA",menu=rna_menu)

        # Erzeuge ein Protein-Menu
        protein_menu = tk.Menu(menubar, tearoff=0)

        # Menu-Einträge
        protein_menu.add_command(label='Translate',command=self.translate)
        protein_menu.add_command(label='Extract proteins',command=self.extract_proteins)
         
        # füge das Menu in die Menubar ein
        menubar.add_cascade(label="Protein",menu=protein_menu)

        # # Erzeuge About-Menu mit debug log toggle
        # about_menu = tk.Menu(menubar, tearoff=0)
        # about_menu.add_command(
        #     label='Activate debug logging',
        #     command=lambda: _switch_gui_logging_level(about_menu)
        # )
        # menubar.add_cascade(label='About',menu=about_menu)

        # den Frame im Container anzeigen
        options = {'padx': 5, 'pady': 5}
        self.pack(**options)

    # Berechnung der Proteine. Es wird der Leserahmen abgefragt, die Proteine berechnet
    # und dann in der GUI angezeigt.
    def extract_proteins(self):
        rframe_no = simpledialog.askinteger("Reading Frames","Welchen Reading Frame wollen Sie benutzen (0-5)?]", parent=self, minvalue=0,maxvalue=5) 
        proteins = self.logic.get_proteins(rframe_no)
        self.prot_text.config(state=tk.NORMAL)
        self.prot_text.delete('1.0','end')
        for p in proteins:
            self.prot_text.insert(END, p+'\n')
        self.prot_text.configure(state ='disabled')

    # Translation
    def translate(self):
        aa_seq = self.logic.translate(0)
        self.aas_text.config(state=tk.NORMAL)
        self.aas_text.delete('1.0','end')
        self.aas_text.insert('1.0', aa_seq)
        self.aas_text.configure(state ='disabled')

    # Faltung der RNA. Es wird die Anzahl der Wasserstoffbrücken angezeigt. Die optimale
    # Struktur als Klammerausdruck und als Zeichnung
    def fold_rna(self):
        opt_anz = self.logic.foldRna()
        self.opt_anzahl.configure(text = str(opt_anz))
        struct = self.logic.trace()
        self.struct_text.config(state=tk.NORMAL)
        self.struct_text.delete('1.0','end')
        self.struct_text.insert('1.0', str(struct))
        self.struct_text.configure(state ='disabled')
        drawRnaStructure(self.logic.get_rna(), struct, self.canvas)

    # Für die Auswertung einer Struktur gibt es einen eigenen Dialog
    def eval_rna(self):
        eval = eval_Dialog(self, "RNA Struktur")

    # Transkription der DNA
    def transcribe(self):
        rna = self.logic.transcribe()
        if not rna is None:
            self.rna_text.config(state=tk.NORMAL)
            self.rna_text.delete('1.0','end')
            self.rna_text.insert('1.0', rna)
            self.rna_text.configure(state ='disabled')

    # Für das Einlesen einer FASTA-Datei wird ein File-Dialog verwendet.
    def read_fasta(self):
        self.dna_text.config(state=tk.NORMAL)
        self.dna_text.delete('1.0','end')
        # öffne einen Dialog zur Dateiauswahl
        filename = fd.askopenfilename()
        # Verwende die Applikation
        dna = self.logic.read_fasta_file(filename)
        self.dna_text.insert('1.0', dna)
        self.dna_text.configure(state ='disabled')
        self.clear_rna()

    # Ein einfacher Dialog für die Eingabe einer DNA
    def read_dna(self):
        self.ask_dialog = Toplevel(self)
        self.ask_dialog.geometry("400x150")
        self.ask_dialog.title("DNA Eingabe")
        ask_label = Label(self.ask_dialog, text = "Bitte DNA eingeben:", anchor = 'w', font = ("Helvetica", 11))
        ask_label.grid(row = 0, column = 0, sticky = 'w',pady = 1, padx = 5)
        self.ask_text = st.ScrolledText(self.ask_dialog, width = 40, height = 3, font = ("Arial", 11))
        self.ask_text.grid(column = 0, pady = 1, padx = 5)
        self.ask_text.bind("<Return>", self.dna_entered)
  
        self.ask_text.insert(tk.INSERT,"""""")

        bottom_frame = tk.Frame(self.ask_dialog)

        ok_button = Button(bottom_frame,text = "Ok", command = self.dna_input_ok)
        cancel_button = Button(bottom_frame,text = "Cancel")
        ok_button.pack(side = 'left')
        cancel_button.pack(side = 'left')

        bottom_frame.grid(row = 3, column = 0, pady = 5, padx = 5)

    # Return-Taste gedrückt
    def dna_entered(self, text):
        self.dna_input_ok()

    # Auf ok gedrückt. Die DNA wird angezeigt.
    def dna_input_ok(self):
        dna = self.ask_text.get("1.0", END)
        dna = dna.strip('\n')
        dna = dna.upper()
        self.logic.set_dna(dna)
        if self.logic.validate_dna():
            self.dna_text.config(state=tk.NORMAL)
            self.dna_text.delete('1.0','end')
            self.dna_text.insert('1.0', dna)
            self.dna_text.configure(state ='disabled')
            self.ask_dialog.destroy()
            self.clear_rna()
        else: 
            tk.messagebox.showinfo(title = "Syntax error", message = "Keine korrekte DNA-Sequenz")

    # Für den Zugriff auf NCBI wird zuerst die email-Adresse und die Accession-Number abgefragt.
    def access_genbank(self):
        ask_email = simpledialog.askstring("Email Adresse", "Bitte email eingeben:") 
        ask_id = simpledialog.askstring("Accession Number", "Bitte accession number eingeben:") 

        dna = self.logic.read_from_ncbi(ask_id, ask_email)
 
        self.dna_text.config(state=tk.NORMAL)
        self.dna_text.delete('1.0','end')
        self.dna_text.insert('1.0', dna)
        self.dna_text.configure(state ='disabled')
        self.clear_rna()
        
    # Nach Eingabe einer neuen DNA sollen alle anderen Anzeigen gelöscht/bereinigt werden.
    def clear_rna(self):
        self.rna_text.config(state=tk.NORMAL)
        self.rna_text.delete('1.0','end')
        self.rna_text.configure(state ='disabled')
        self.opt_anzahl.config(text = "")
        #self.opt_struct.config(text = "")
        self.struct_text.config(state=tk.NORMAL)
        self.struct_text.delete('1.0','end')
        self.struct_text.configure(state ='disabled')
        self.canvas.delete(ALL)
        self.aas_text.config(state=tk.NORMAL)
        self.aas_text.delete('1.0','end')
        self.aas_text.configure(state ='disabled')
        self.prot_text.config(state=tk.NORMAL)
        self.prot_text.delete('1.0','end')
        self.prot_text.configure(state ='disabled')

    
# Hauptfunktion
if __name__ == "__main__":
    log.info("Application started.")
    app = App()
    # erzeuge ein biologic-Objekt und verwende dies in der App
    logic = biologic()
    frame = SeqMainFrame(app, logic)
    app.mainloop()
