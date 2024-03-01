#----------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------Les données et modules---------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------

import sys
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QLineEdit, QPushButton, QTextEdit, QVBoxLayout, QFileDialog, \
    QMessageBox, QInputDialog
import random
import numpy as np

#table de correspondance des codons ARN et leurs acides aminés
rna_codon_table = {
    "AAA": "Lysine", "AAC": "Asparagine", "AAG": "Lysine", "AAU": "Asparagine",
    "ACA": "Threonine", "ACC": "Threonine", "ACG": "Threonine", "ACU":"Threonine",
    "AGA": "Arginine", "AGC": "Serine", "AGG": "Arginine", "AGU": "Serine",
    "AUA": "Isoleucine", "AUC": "Isoleucine", "AUG": "Methionine", "AUU": "Isoleucine",
    "CAU": "Histidine", "CAC": "Histidine", "CAG": "Glutamine",
    "CCA": "Proline", "CCC": "Proline", "CCG": "Proline","CCU":"Proline",
    "CGA": "Arginine", "CGC": "Arginine", "CGG": "Arginine", "CGU": "Arginine",
    "CUA": "Leucine", "CUC": "Leucine", "CUG": "Leucine", "CUU": "Leucine",
    "GAA": "Glutamic acid", "GAC": "Glutamic acid", "GAG": "Glutamic acid","GAU":"Aspatic acid",
    "GCA": "Alanine", "GCC": "Alanine", "GCG": "Alanine","GCU":"Alanine",
    "GGA": "Glycine", "GGC": "Glycine", "GGG": "Glycine", "GGU": "Glycine",
    "GUA": "Valine", "GUC": "Valine", "GUG": "Valine", "GUU": "Valine",
    "UAA": "STOP", "UAC": "Tyrosine", "UAG": "STOP", "UAU": "Tyrosine",
    "UCA": "Serine", "UCC": "Serine", "UCG": "Serine", "UCU": "Serine",
    "UGA": "STOP", "UGC": "Cysteine", "UGG": "Tryptophan", "UGU": "Cysteine",
    "UUA": "Leucine", "UUC": "Leucine", "UUG": "Leucine", "UUU": "Leucine","CAA":"Glutamine"
}


#----------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------Les fonctionalités-------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------


#fonction qui génére une chaine d'ADN aléatoire
def generate_dna(length):
    dna=""
    for _ in range(length):
        dna+=random.choice("ATCG")
    return dna

#fonction qui choisi des chaines d'ADN a partir d'un fichier Fasta
def read_fasta_file(file_path):
    sequences = {}
    current_sequence_name = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()

            if line.startswith('>'):  
                current_sequence_name = line[1:] 
                sequences[current_sequence_name] = "" 
            else:
                sequences[current_sequence_name] += line 

    return sequences


#fonction qui verifie la validité d'une chaine d'ADN
def is_valid(dna):
    valid=set("ATCG")
    for nuc in dna:
        if nuc not in valid:
            return False
    return True


#fonction qui calcule la fréquence des base nucléiques
def freq_nuc(dna):
    freq={}
    for nuc in "ATCG":
        freq[nuc]=dna.count(nuc)
    return freq


#fonction qui transcrie une chaine d'ADN en ARN
def from_dna_to_rna(dna):
    rna=dna.replace("T","U")
    return rna


#Fonction qui transcrie une chaine d'ARN en protéines
def from_rna_to_protein(rna):
    codons=[rna[i:i+3] for i in range(0,len(rna),3)]
    prot=""
    for codon in codons :
        prot=prot+"-"+rna_codon_table[codon]
    return prot


#fonction qui calcule le complément inverse d'une chaine d'ADN
def complement_inverse_dna(dna):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement[base] for base in dna[::-1])


#fonction qui calcule le taux de GC dans une chaine d'ADN
def gc_content(dna):
    return((dna.count("C")+dna.count("G"))/len(dna))*100


#fonction qui calcule la fréquence des codons dans une chaine d'ADN
def codon_frequencies(dna):
    codons = [dna[i:i+3] for i in range(0, len(dna), 3)]
    frequencies = {}
    for codon in codons:
            if len(codon)==3:
                if codon not in frequencies:
                    frequencies[codon] = 0
                frequencies[codon]+= 1
    return frequencies


#fonction qui fait des mutations aléatoires sur une chaine d'ADN
def mutate_dna(dna, num_mutations):
    mutated_dna = list(dna)
    for _ in range(num_mutations):
        pos = random.randint(0, len(dna) - 1)
        new_nucleotide = None
        while new_nucleotide is None or new_nucleotide == mutated_dna[pos]:
            new_nucleotide = random.choice("ATCG")
        mutated_dna[pos] = new_nucleotide
    mutated_dna = "".join(mutated_dna)
    return mutated_dna


#fonction qui chercheun motif dans une chaine d'ADN et affiche les positions du motif
def motif_search(dna,motif):
    if len(motif)>len(dna):
        return "la taille du motif ne peut pas dépasser la taille de la séquance d'ADN"
    else:
        pos=""
        for i in range(0,len(dna),1):
            if dna[i:i+len(motif)]==motif.upper():
                pos=pos+" "+str(i+1)
        return pos
    

#fonction qui génère la chaine d'ADN consensus et la matrice profil a partir d'un fichier fasta
def profile_mat(fasta_file):
    # Lecture du fichier FASTA
    with open(fasta_file, "r") as file:
        lines = file.readlines()

    # Initialisation des variables
    dna_sequences = []
    l = 0

    # Extraction des séquences d'ADN du fichier FASTA
    for line in lines:
        if line.startswith('>')and is_valid(lines[l+1].strip()):
            # Ajout de la séquence actuelle à la liste
            dna_sequences.append(lines[l + 1].strip())
        l += 1

    # Construction de la matrice profil
    mat = np.array([list(seq) for seq in dna_sequences], dtype=str)
    profile_mat = np.zeros((4, len(dna_sequences[0])))

    for i, base in enumerate("ACGT"):
        profile_mat[i, :] = np.sum(mat == base, axis=0)

    # Construction de la chaîne d'ADN consensus
    bases = "ACGT"
    index = np.argmax(profile_mat, axis=0)
    consensus = ''.join(bases[i] for i in index)

    return consensus, profile_mat

#fonction qui verifie si un fichier est de type Fasta
def is_fasta_file(file_path):
    try:
        with open(file_path, 'r') as file:
            first_line = file.readline().strip()

            # Vérifier si la première ligne commence par ">"
            if not first_line.startswith('>'):
                return False

            # Vérifier si les lignes suivantes contiennent une séquence ADN valide
            for line in file:
                line = line.strip()
                if not line:
                    continue 
                elif line.startswith('>'):
                    return True  
                elif not line[0].upper() in {'A', 'C', 'G', 'T'}:
                    return False  

            # Si toutes les conditions sont remplies, le fichier est probablement au format FASTA avec séquence ADN
            return True
    except FileNotFoundError:
        return False  
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier : {e}")
        return False

#----------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------La classe principale-----------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------

class DNAAnalyzerApp(QWidget):
    #constructeur de la class qui herite de la class QWidget
    def __init__(self):
        super().__init__()
        
        #Les variables qui vont contenir le path du fichier fasta , la sequence/les séquences d'ADN , et si ils ont étaient importés d'un fichier
        
        self.file_path = ""
        self.dna_sequence = ""
        self.sequences={}
        self.from_file=False

        self.init_ui()

    #fonction qui s'occupe de l'interface 
    def init_ui(self):
        layout = QVBoxLayout()
        
        # Les boutons pour générer ou importer une/des chaines d'ADN
        
        self.start_label = QLabel("----------------------        Choisisez entre:", self)
        layout.addWidget(self.start_label)

        choose_file_button = QPushButton("1- Importer la chaine d'ADN d'un fichier FASTA", self)
        choose_file_button.clicked.connect(self.choose_file)
        layout.addWidget(choose_file_button)

        random_dna_button = QPushButton("2- Générer une chaine d'ADN aléatoire", self)
        random_dna_button.clicked.connect(self.generate_random_dna)
        layout.addWidget(random_dna_button)

        # Le champ de text qui va afficher la/les chaines d'ADN
        self.dna_label = QLabel("----------------------        Votre chaine d'ADN:", self)
        self.dna_label.hide()
        layout.addWidget(self.dna_label)

        self.dna_entry = QTextEdit(self)
        self.dna_entry.setPlaceholderText("")
        self.dna_entry.hide() 
        layout.addWidget(self.dna_entry)


        # Les boutons pour les 10 fonctionalités
        self.functions_label = QLabel("----------------------        Choisissez entre ces fonctions:", self)
        self.functions_label.hide()
        layout.addWidget(self.functions_label)

        self.validate_button = QPushButton("1- Vérifier la validité de la chaîne ADN", self)
        self.validate_button.clicked.connect(self.validate_sequence)
        self.validate_button.hide()
        layout.addWidget(self.validate_button)

        self.freq_nuc_button = QPushButton("2- Calculer les fréquences des bases nucléiques dans la chaîne ADN", self)
        self.freq_nuc_button.clicked.connect(self.calculate_nucleotide_frequencies)
        self.freq_nuc_button.hide()
        layout.addWidget(self.freq_nuc_button)

        self.transcribe_button = QPushButton("3- Transcrire la chaîne ADN en une chaîne ARN", self)
        self.transcribe_button.clicked.connect(self.transcribe_to_rna)
        self.transcribe_button.hide()
        layout.addWidget(self.transcribe_button)

        self.translate_button = QPushButton("4- Transcrire la chaîne ARN résultante en protéines", self)
        self.translate_button.clicked.connect(self.translate_to_protein)
        self.translate_button.hide()
        layout.addWidget(self.translate_button)

        self.reverse_complement_button = QPushButton("5- Calculer le complément inverse de la chaîne ADN", self)
        self.reverse_complement_button.clicked.connect(self.calculate_reverse_complement)
        self.reverse_complement_button.hide()
        layout.addWidget(self.reverse_complement_button)

        self.gc_content_button = QPushButton("6- Calculer le taux de GC de la séquence ADN", self)
        self.gc_content_button.clicked.connect(self.calculate_gc_content)
        self.gc_content_button.hide()
        layout.addWidget(self.gc_content_button)

        self.codon_freq_button = QPushButton("7- Calculer les fréquences de codons dans la chaîne ADN", self)
        self.codon_freq_button.clicked.connect(self.calculate_codon_frequencies)
        self.codon_freq_button.hide()
        layout.addWidget(self.codon_freq_button)

        self.mutate_button = QPushButton("8- Réaliser des mutations aléatoires sur la chaîne ADN", self)
        self.mutate_button.clicked.connect(self.mutate_dna)
        self.mutate_button.hide()
        layout.addWidget(self.mutate_button)

        self.motif_search_button = QPushButton("9- Chercher un motif dans la chaîne ADN", self)
        self.motif_search_button.clicked.connect(self.search_motif)
        self.motif_search_button.hide()
        layout.addWidget(self.motif_search_button)

        self.profile_mat_button = QPushButton("10- Générer la chaîne ADN consensus et la matrice profil", self)
        self.profile_mat_button.clicked.connect(self.generate_profile_mat)
        self.profile_mat_button.hide()
        layout.addWidget(self.profile_mat_button)


        # Le champ de text qui va afficher les résultats des fonctionalités
        self.result_text = QTextEdit(self)
        self.result_text.setReadOnly(True)
        self.result_text.setPlaceholderText("Les resultats vont s'afficher ici!")
        layout.addWidget(self.result_text)


        # Le bouton qui permettra a l'utilisateur de sauvgarder les resultats dans un fichier
        save_button = QPushButton("Sauvegarder les résultats dans un fichier", self)
        save_button.clicked.connect(self.save_results)
        layout.addWidget(save_button)

        self.setLayout(layout)

        self.setGeometry(100, 100, 800, 600)  # Set window size here
        self.setWindowTitle('DNA Analyzer')
        self.show()

    # Fonction qui permet d'afficher les fonctionalités aprés avoir générer/importer une/des chaines d'ADN
    def showall(self):
        self.functions_label.show()
        self.gc_content_button.show()
        self.mutate_button.show()
        self.freq_nuc_button.show()
        self.validate_button.show()
        self.translate_button.show()
        self.codon_freq_button.show()
        self.transcribe_button.show()
        self.reverse_complement_button.show()
        self.motif_search_button.show()
        self.profile_mat_button.show()
        self.dna_entry.show()
        self.dna_label.show()
    

    # Fonctions qui permet de lire les séquences d'ADN d'un fichier
    def choose_file(self):
        self.file_path, _ = QFileDialog.getOpenFileName(self, 'Choose DNA File', '', 'Text Files (*.txt);;All Files (*)')
        if is_fasta_file(self.file_path):
            if self.file_path:
                self.sequences = read_fasta_file(self.file_path)

                # Afficher les noms des séquences dans un message box
                sequence_names = "\n".join(self.sequences.keys())
                QMessageBox.information(self, 'Séquences chargées', f'Les séquences suivantes ont été chargées :\n{sequence_names}')
                result=""
                for name,sequence in self.sequences.items():
                    if is_valid(sequence):
                        result=result+f"{name}:\n{sequence}\n"
                self.dna_entry.setPlainText(result)
                self.validate_sequence(validate_all=True)
                self.showall()
                self.from_file=True 
        else:
            self.show_error("il faut que le fichier soit de type FASTA!!")
                
    # Fonction qui génère une séquence d'ADN aléatoire
    def generate_random_dna(self):
        length, ok = QInputDialog.getInt(self, "Générer chaine aléatoire d'ADN ", 'Veuillez entrez la taille de la chaine: ')
        if ok:
            self.dna_sequence = generate_dna(length)
            self.dna_entry.setPlainText(self.dna_sequence)
            self.validate_sequence()
            self.showall()
            self.from_file=False

    """
    Chaqu'une des fonctions ci-dessous :
        1-verifie d'abord si la/les sequences ont étaient importés d'un fichier.
        2-puis fait ce qu'elle a faire
    de cette maniére les fonctions peuvent traiter chaque séquence d'un fichier fasta séparément(dans le cas ou il ya plusieurs!)

    """
    def validate_sequence(self, validate_all=False):
        if validate_all or self.from_file:
            if not self.sequences:
                self.show_error("Veuillez charger des séquences à partir d'un fichier FASTA d'abord!")
                return

            validation_results = []
            for sequence_name, sequence in self.sequences.items():
                if is_valid(sequence):
                    validation_results.append(f"La séquence {sequence_name} est valide.")
                else:
                    validation_results.append(f"La séquence {sequence_name} est invalide!!!!!!!!!")
            result_text = "\n".join(validation_results)
            QMessageBox.information(self,"Validation",f"{result_text}")
        else:
            if self.dna_sequence:
                if is_valid(self.dna_sequence):
                    QMessageBox.information(self, 'Validation', 'La séquence d\'ADN est valide!')
                else:
                    QMessageBox.warning(self, 'Validation', "Chaine d'ADN invalide!!!!")
            else:
                self.show_error("Veuillez entrer une chaine MERDE d'ADN d'abord!")

    def calculate_nucleotide_frequencies(self):
        if self.from_file:
            results=[]
            for name,sequence in self.sequences.items():
                if is_valid(sequence):
                    nucleotide_frequencies = freq_nuc(sequence)
                    text=""
                    for nuc,freq in nucleotide_frequencies.items():
                        text+=f"{nuc} : {freq}\n"
                    results.append(f"-la fréquence des nuclétides de {name}:\n{text}")
            final_result="\n".join(results)
            self.display_result(final_result)
        elif self.dna_sequence:
            nucleotide_frequencies = freq_nuc(self.dna_sequence)
            print(nucleotide_frequencies)
            result="Fréquences des nucléotides :\n"
            for nuc,freq in nucleotide_frequencies.items():
                result+=f"{nuc} : {freq}\n"
            self.display_result(result)
        else:
            self.show_error("Veuillez entrez une chaine d'ADN d'abord!")

    def calculate_reverse_complement(self):
        if self.from_file:
            results=[]
            for name,sequence in self.sequences.items():
                if is_valid(sequence):
                    results.append(f"-le complément inverse de la chaine d'ADN '{name}' est: {complement_inverse_dna(sequence)}")
            self.display_result("\n".join(results))
        elif self.dna_sequence:
            reverse_complement_sequence = complement_inverse_dna(self.dna_sequence)
            self.display_result(f"Le complément inverse de la chaine d'ADN: {reverse_complement_sequence}")
        else:
            self.show_error("Veuillez entrez une chaine d'ADN d'abord!")

    def calculate_codon_frequencies(self):
        if self.from_file:
            results=[]
            for name,sequence in self.sequences.items():
                if is_valid(sequence):
                    codon_frequencie = codon_frequencies(sequence)
                    text=""
                    for codon,freq in codon_frequencie.items():
                        text+=f"{codon} : {freq}\n"
                    results.append(f"-la fréquence des codons de {name}:\n{text}") 
            final_result='\n'.join(results)
            self.display_result(final_result)

        elif self.dna_sequence:
            if is_valid(self.dna_sequence):
                codon_frequencies_result = codon_frequencies(self.dna_sequence)
                result="La fréquence des codons:\n"
                for codon,freq in codon_frequencies_result.items():
                    result=result+codon+":"+str(freq)+"\n"
                self.display_result(result)
            else:
                self.show_error("Chaine d'ADN invalide!!!!")
        else:
            self.show_error("Veuillez entrez une chaine d'ADN d'abord!")

    def transcribe_to_rna(self):
        if self.from_file:
            results=[]
            for name,sequence in self.sequences.items():
                if is_valid(sequence):
                    results.append(f"-La chaine d'ARN aprés transcription de la chaine {name} :{from_dna_to_rna(sequence)}")
            self.display_result("\n".join(results))
        elif self.dna_sequence:
            rna_sequence = from_dna_to_rna(self.dna_sequence)
            self.display_result(f"La chaine d'ARN aprés transcription: {rna_sequence}")
        else:
            self.show_error("Veuillez entrez une chaine d'ADN d'abord!")

    def calculate_gc_content(self):
        if self.from_file:
            results=[]
            for name,sequence in self.sequences.items():
                if is_valid(sequence):
                    results.append(f"-Le taux de GC de la chaine d'ADN {name}: {gc_content(sequence):.2f}%")
            self.display_result("\n".join(results))
        elif self.dna_sequence:
            if is_valid(self.dna_sequence):
                gc_content_result = gc_content(self.dna_sequence)
                self.display_result(f"Le taux de GC de la chaine d'ADN: {gc_content_result:.2f}%")
            else: 
                self.show_error("Chaine d'ADN invalide!!!!")
        else:
            self.show_error("Veuillez entrez une chaine d'ADN d'abord!")
    
    def mutate_dna(self):
        num_mutations, ok = QInputDialog.getInt(self, "Mutation d'ADN", 'Veuillez entez le nombre de mutations: ')
        if ok:
            if self.from_file:
                results=[]
                for n,s in self.sequences.items():
                    if is_valid(s):
                        results.append(f"-La chaine d'ADN {n} aprés {num_mutations} mutations :{mutate_dna(s,num_mutations)} ")
                self.display_result("\n".join(results))
            elif self.dna_sequence:
                if is_valid(self.dna_sequence):
                    mutated_dna_sequence = mutate_dna(self.dna_sequence, num_mutations)
                    self.display_result(f"La chaine d'AN aprés {num_mutations} mutations : {mutated_dna_sequence}")
                else:
                    self.show_error("Chaine d'ADN invalide!!!!")

    def search_motif(self):
        motif, ok = QInputDialog.getText(self, 'Chercher un motif', 'Veuillez entrez le motif:')
        if ok:
            if is_valid(motif.upper()):
                if self.from_file:
                    results=[]
                    errors=[]
                    for n,s in self.sequences.items():
                        if is_valid(s):
                            if len(motif)<len(s):
                                results.append(f"-Les positions du motif {motif} dans la chaine d'ADN {n}:{motif_search(s,motif)} ")
                            else:
                                errors.append(f"- Le motif est plus long que la chaine d'ADN {n}!!!!")                  
                    self.display_result("\n".join(results))
                    self.show_error("\n".join(errors))
                elif self.dna_sequence:
                    if is_valid(self.dna_sequence):
                        if len(motif)<len(self.dna_sequence):
                            motif_positions = motif_search(self.dna_sequence, motif)
                            self.display_result(f"Les Postions du motif {motif}:  {motif_positions}")
                        else:
                            self.show_error("le motif ne peut pas etre plus long que la chaine d'ADN !!!!")
                    else:
                        self.show_error("Chaine d'ADN invalide!!!!")
            else:
                self.show_error("Le motif est invalid!!")

    def generate_profile_mat(self):
        if self.file_path:
            first_length=len(next(iter(self.sequences.values()),""))
            is_equal=all(len(value)==first_length for value in self.sequences.values())
            if is_equal:
                consensus, profile_matrix = profile_mat(self.file_path)
                self.display_result(f"La chaîne ADN consensus: {consensus}\n\nLa matrice profil:\n{profile_matrix}")
            else:
                self.show_error("Les séquences doivent etre de la meme taille!!")
        else:
            self.show_error("Veuilez d'abord choisir un fichier FASTA!!!")
    
    def translate_to_protein(self):
        if self.from_file:
            results=[]
            for n,s in self.sequences.items():
                if is_valid(s):
                    results.append(f"-La séquence de protéines de la chaine {n} :{from_rna_to_protein(from_dna_to_rna(s))}")
            self.display_result("\n".join(results))
        elif self.dna_sequence:
            if is_valid(self.dna_sequence):
                protein_sequence = from_rna_to_protein(from_dna_to_rna(self.dna_sequence))
                self.display_result(f"La séquence de Protéines:  {protein_sequence}")
            else:
                self.show_error("Chaine d'ADN invalide!!!!")
        else:
            self.show_error("Veuillez entrez une chaine d'ADN d'abord!")

    # Fonction qui affiche les resultats dans le champ de resultats
    def display_result(self, result):
        current_content = self.result_text.toPlainText()
        new_content = current_content + "\n\n" + result
        self.result_text.setPlainText(new_content)


    #fonction qui enregistre les resultats dans un fichier. 
    def save_results(self):
        result_text = self.result_text.toPlainText()
        file_path, _ = QFileDialog.getSaveFileName(self, 'Save Results', '', 'Text Files (*.txt);;All Files (*)')
        if file_path:
            with open(file_path, 'w') as file:
                file.write(result_text)
            QMessageBox.information(self, 'Success', 'Results saved successfully!')

    def show_error(self, message):
        QMessageBox.critical(self, 'Error', message)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = DNAAnalyzerApp()
    sys.exit(app.exec_())