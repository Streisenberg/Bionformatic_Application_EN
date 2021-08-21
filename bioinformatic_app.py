import streamlit as st
import sys
import os
from Bio.Seq import Seq 
from Bio import SeqIO
from Bio.SeqUtils import GC
import neatbio.sequtils as utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import Counter
from Bio.Data import CodonTable
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Emboss.Applications import NeedleCommandline
from Bio import pairwise2, Align
from Bio.Align import substitution_matrices
#from PIL import Image
from chembl_webresource_client.new_client import new_client
import base64
from SessionState import SessionState
#Aşağıdaki libraryler herokuya yüklenmeden önce requirements dosyasına eklenecek
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import seaborn as sns
from numpy.random import seed
from numpy.random import randn
from scipy.stats import mannwhitneyu
import json
import pickle
import uuid
import re
from functions import download_button
#import PyPDF2
from PIL import Image



st.set_option('deprecation.showPyplotGlobalUse', False)

sns.set(style='ticks')

sys.path.append('/usr/local/lib/python3.7/site-packages/')

st.set_option('deprecation.showfileUploaderEncoding', False)

session_state = SessionState.get(name="", button_sent=False)


def mannwhitney(descriptor, verbose=False):

    # seed the random number generator
    seed(1)

    # actives and inactives
    selection = [descriptor, 'bioactivity_class']
    df = df_class[selection]
    active = df[df.bioactivity_class == 'active']
    active = active[descriptor]

    selection = [descriptor, 'bioactivity_class']
    df = df_class[selection]
    inactive = df[df.bioactivity_class == 'inactive']
    inactive = inactive[descriptor]

    # compare samples
    stat, p = mannwhitneyu(active, inactive)
    #print('Statistics=%.3f, p=%.3f' % (stat, p))

    # interpret
    alpha = 0.05
    if p > alpha:
        interpretation = 'Same distribution (fail to reject H0)'
    else:
        interpretation = 'Different distribution (reject H0)'
    
    results = pd.DataFrame({'Descriptor':descriptor,
                            'Statistics':stat,
                            'p':p,
                            'alpha':alpha,
                            'Interpretation':interpretation}, index=[0])
    filename = 'mannwhitneyu_' + descriptor + '.csv'
    results.to_csv(filename)

    return results

def pIC50(input):
    pIC50 = []

    for i in input['standard_value_norm']:
        molar = i*(10**-9) # Converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)
        
    return x

def norm_value(input):
    norm = []
    for i in input['standard_value']:
        i = float(i)
        i = int(i)
        if i > 100000000:
          i = 100000000
        norm.append(i)

    input['standard_value_norm'] = norm
    x = input.drop('standard_value', 1)
        
    return x

# Inspired by: https://codeocean.com/explore/capsules?query=tag:data-curation
def lipinski(smiles, verbose=False):


    moldata= []
    for elem in smiles:

        mol=Chem.MolFromSmiles(elem) 
        moldata.append(mol)

    baseData= np.arange(1,1)
    i=0  
    for mol in moldata:        
       
      desc_MolWt = Descriptors.MolWt(mol)
      desc_MolLogP = Descriptors.MolLogP(mol)
      desc_NumHDonors = Lipinski.NumHDonors(mol)
      desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
          
      row = np.array([desc_MolWt,
                      desc_MolLogP,
                      desc_NumHDonors,
                      desc_NumHAcceptors])   
  
      if(i==0):
          baseData=row
      else:
          baseData=np.vstack([baseData, row])
      i=i+1      
  
    columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]   
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)
    
    return descriptors

    #You can check it out rdkit library to further reading for these functions

def get_table_download_link(df):
    """Generates a link allowing the data in a given panda dataframe to be downloaded
    in:  dataframe
    out: href string
    """
    pdf = PyPDF2.PdfFileReader(df)
    b64 = base64.b64encode(pdf.encode()).decode()  # some strings <-> bytes conversions necessary here
    href = f'<a href="data:file/pdf;base64,{b64}" download="CSV_Dosyası.csv">CSV Dosyasını Bilgisayarına Yükle</a>'

    return href

def get_binary_file_downloader_html(bin_file, file_label='File'):

    with open(bin_file, 'rb') as f:
        data = f.read()
    bin_str = base64.b64encode(data).decode()
    href = f'<a href="data:application/octet-stream;base64,{bin_str}" download="{os.path.basename(bin_file)}">Download {file_label}</a>'
    return href

def gc_content(seq):
    result = GC(seq)
    return result

def at_content(seq):
    result = float(str(seq).count("A") + str(seq).count("T")) / len(seq) * 100
    return result

def plot(seq, adenin_color, timin_color, guanin_color, sitozin_color):
    nucleotides = ["Adenin", "Timin", "Guanin", "Sitozin"]
    df = pd.DataFrame(seq, index=nucleotides)
    barlist = plt.bar(nucleotides, seq, label="Nucleotides", width=.2)
    barlist[0].set_color(adenin_color)
    barlist[1].set_color(timin_color)
    barlist[2].set_color(guanin_color)
    barlist[3].set_color(sitozin_color)
    
    st.pyplot()

def main():

    st.title("Bioinformatic App")

    st.sidebar.header("Please Choose")
    secenekler = ["Home","DNA Sequencing", "Multiple Sequence Alignment", "Sequence Scoring", "Target Protein Analysis", "Find SNP and INDEL"]
    select_box = st.sidebar.selectbox("The action you want to do: ", secenekler)

    if select_box == "Home":
        st.header("Welcome to the Home Page")
        st.text("")
        st.write("There will be information section in here. Please Patient (:")
        st.text("")
        st.text("")
        st.write("Made by **_Enes Daşdemir_**")
    elif select_box == "DNA Sequencing":
        st.subheader("DNA Sequence Analysis")

        st.warning("Please import your DNA Sequence from the file upload section in the left bar.")

        seq_dosya = st.sidebar.file_uploader("Please upload your FASTA file", type=["FASTA","fa"])
        
        if seq_dosya is not None:
            dna = SeqIO.read(seq_dosya, "fasta")
            st.write(dna)

            dna_sekansi = dna.seq

            details = st.radio("Details", ("Explanation", "Show the sequence"))
            if details == "Explanation":
                st.text("")
                st.write(dna.description)
            elif details == "Show the sequence":
                st.text("")
                st.write(dna.seq)

            st.text("")
            st.text("")

            st.subheader("Nucleotids information")

            st.text("")


            if ("M" and "L") in str(dna_sekansi):
                st.write("Please enter a **DNA sequence** so that the nucleotide information can be calculated!")
                

            else:
                
                adenin = int(str(dna_sekansi).count("A"))
                guanin = int(str(dna_sekansi).count("G"))
                citosin = int(str(dna_sekansi).count("C"))
                timin = int(str(dna_sekansi).count("T"))
                st.write("Number of **Adenine** = {0} ".format(adenin))
                st.write("Number of **Thymine** = {0} ".format(timin))
                st.write("Number of **Guanine** = {0} ".format(guanin))
                st.write("Number of **Cytosine** = {0} ".format(citosin))

                st.text("")
                st.text("")

                if st.checkbox("Show on Graphics"):
                    adenin_color = st.beta_color_picker('Choose a color for Adenine', "#F50000")
                    timin_color = st.beta_color_picker('Choose a color for Thymine', "#00DE2D")
                    guanin_color = st.beta_color_picker('Choose a color for Guanine', "#1A00FF")
                    sitozin_color = st.beta_color_picker('Choose a color for Cytosine', "#000000")
                    
                    numbers = [adenin, timin, guanin, citosin]
                    plot(numbers, adenin_color, timin_color, guanin_color, sitozin_color)
                
                st.text("")

                st.subheader("Ingredient rates")
                st.text("")

                gc_orani = round(gc_content(dna_sekansi), 2)
                at_orani = round(at_content(dna_sekansi),2)

                st.write("**GC** ratio = % {0}".format(gc_orani))
                st.write("**AT** ratio = % {0}".format(at_orani))

                st.text("")

                st.subheader("Protein Synthesis")
                aa = dna_sekansi.translate()
                aa_frekansi = Counter(str(aa))
                st.text("")

                if st.button("Click to see Translation table"):
                    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
                    st.text(standard_table)
                
                st.text("")

                if st.checkbox("Transcription"):
                    transkribe = dna_sekansi.transcribe()
                    st.write(transkribe)
                elif st.checkbox("Translation"):
                    transle = dna_sekansi.translate()
                    st.write(transle)
                elif st.checkbox("Complement"):
                    st.write(dna_sekansi.complement())
                elif st.checkbox("Frequency of Amino Acid"):
                    st.write(aa_frekansi)
                elif st.checkbox("Graph of Amino Acid"):
                    st.text("")
                    aa_color = st.beta_color_picker("You can change color from here")
                    plt.bar(aa_frekansi.keys(), aa_frekansi.values(), color = aa_color)
                    st.pyplot()
                elif st.checkbox("Full Amino Acid Names"):
                    st.text("")
                    aa_ismi = str(aa).replace("*", "")
                    aa3 = utils.convert_1to3(aa_ismi)
                    st.write("**Letter Display**")
                    st.text(aa_ismi)
                    st.write("**************************************")
                    st.write("**Abbreviation Display**")
                    st.text(aa3)
                    st.write("**************************************")
                    st.write("**Open Name Display**")
                    st.text(utils.get_acid_name(aa3))              

    elif select_box == "Multiple Sequence Alignment":

        st.warning("You can load the files you want to compare from the bar on the left.")
        seq1 = st.sidebar.file_uploader("1.FASTA File", type=["fasta", "fa"])

        seq2 = st.sidebar.file_uploader("2.FASTA File", type=["fasta", "fa"])

        if seq1 and seq2 is not None:
            sekans1 = SeqIO.read(seq1, "fasta")
            sekans2 = SeqIO.read(seq2, "fasta")

            st.text("")

            st.write("**Sequence 1** = {0}".format(sekans1.description))
            st.text("")
            st.write("**Sequence 2** = {0}".format(sekans2.description))
            st.text("")
            st.write("**Aligned state:**")

            alignments = pairwise2.align.globalxx(sekans1.seq, sekans2.seq)
            st.text(pairwise2.format_alignment(*alignments[0]))
            st.text("")

            if st.checkbox("Click to see BLOSUM62 Rating"):

                secenek = st.radio("Seçenekler",("Aligned Sequence","All Sequence"))

                if secenek == "Aligned Sequence":

                    blosum62 = substitution_matrices.load("BLOSUM62")
                    alignment_blossum = pairwise2.align.localds(sekans1.seq, sekans2.seq, blosum62, -10, -1)
                    st.text(pairwise2.format_alignment(*alignment_blossum[0]))
                    
                elif secenek == "All Sequence":
                    blosum62_2 = substitution_matrices.load("BLOSUM62")
                    full_alignment_blossum = pairwise2.align.localds(sekans1.seq, sekans2.seq, blosum62_2, -10, -1)
                    st.text(pairwise2.format_alignment(*full_alignment_blossum[0], full_sequences=True))

                  #Bu kısımda geliştirme olarak gap-penalty ve gap-extension gibi değerleri kullanıcının değiştirebileceği gibi ayarlayabiliriz. 
                  #Geliştirme olarak sekansları taşımak yerine yazılabilir hale de getirebilirim!!!!!


            st.write("Made by **_Enes Daşdemir_**")
    elif select_box == "Sequence Scoring":

        
        st.text("")
        st.text("")
        st.text("")
        st.subheader("Do your own calculation: ")
        
        seq1_puan = st.text_area("Sequence 1")
        seq2_puan = st.text_area("Sequence 2")
        substitution_matrices.load()
        option = st.selectbox('Matrix you want to use?',("BENNER22", 'BENNER6', 'BENNER74', 'BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'DAYHOFF', 'FENG', 'GENETIC', 'GONNET1992', 'JOHNSON', 'JONES', 'LEVIN', 'MCLACHLAN', 'MDM78', 'NUC.4.4', 'PAM250', 'PAM30', 'PAM70', 'RAO', 'RISLER', 'SCHNEIDER', 'STR'))
        st.write('Matrix you choose:', option)
        try:        
            aligner = Align.PairwiseAligner()
            if option == "BENNER22":
                matrix = substitution_matrices.load("BENNER22")
                st.text(matrix)
                aligner.substitution_matrix = matrix
                score = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score)) 
            elif option == "BLOSUM62":
                matrix2 = substitution_matrices.load("BLOSUM62")
                st.text(matrix2)
                aligner.substitution_matrix = matrix2
                score2 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score2))
            elif option == "BENNER6":
                matrix3 = substitution_matrices.load("BENNER6")
                st.text(matrix3)
                aligner.substitution_matrix = matrix3
                score3 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score3))
            elif option == "BENNER74":
                matrix4 = substitution_matrices.load("BENNER74")
                st.text(matrix4)
                aligner.substitution_matrix = matrix4
                score4 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score4))
            elif option == "BLOSUM45":
                matrix5 = substitution_matrices.load("BLOSUM45")
                st.text(matrix5)
                aligner.substitution_matrix = matrix5
                score5 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score5))
            elif option == "BLOSUM50":
                matrix6 = substitution_matrices.load("BLOSUM50")
                st.text(matrix6)
                aligner.substitution_matrix = matrix6
                score6 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score6))
            elif option == "BLOSUM80":
                matrix7 = substitution_matrices.load("BLOSUM80")
                st.text(matrix7)
                aligner.substitution_matrix = matrix7
                score7 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score7))
            elif option == "BLOSUM90":
                matrix8 = substitution_matrices.load("BLOSUM90")
                st.text(matrix8)
                aligner.substitution_matrix = matrix8
                score8 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score8))
            elif option == "DAYHOFF":
                matrix9 = substitution_matrices.load("DAYHOFF")
                st.text(matrix9)
                aligner.substitution_matrix = matrix9
                score9 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score9))
            elif option == "FENG":
                matrix10 = substitution_matrices.load("FENG")
                st.text(matrix10)
                aligner.substitution_matrix = matrix10
                score10 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score10))
            elif option == "GENETIC":
                matrix11 = substitution_matrices.load("GENETIC")
                st.text(matrix11)
                aligner.substitution_matrix = matrix11
                score11 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score11))
            elif option == "GONNET1992":
                matrix12 = substitution_matrices.load("GONNET1992")
                st.text(matrix12)
                aligner.substitution_matrix = matrix12
                score12 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score12))
            elif option == "JOHNSON":
                matrix13 = substitution_matrices.load("JOHNSON")
                st.text(matrix13)
                aligner.substitution_matrix = matrix13
                score13 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score13))
            elif option == "JONES":
                matrix14 = substitution_matrices.load("JONES")
                st.text(matrix14)
                aligner.substitution_matrix = matrix14
                score14 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score14))
            elif option == "LEVIN":
                matrix15 = substitution_matrices.load("LEVIN")
                st.text(matrix15)
                aligner.substitution_matrix = matrix15
                score15 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score15))
            elif option == "MCLACHLAN":
                matrix16 = substitution_matrices.load("MCLACHLAN")
                st.text(matrix16)
                aligner.substitution_matrix = matrix16
                score16 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score16))
            elif option == "MDM78":
                matrix17 = substitution_matrices.load("MDM78")
                st.text(matrix17)
                aligner.substitution_matrix = matrix17
                score17 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score17))
            elif option == "NUC.4.4":
                matrix18 = substitution_matrices.load("NUC.4.4")
                st.text(matrix18)
                aligner.substitution_matrix = matrix18
                score18 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score18))
            elif option == "PAM250":
                matrix19 = substitution_matrices.load("PAM250")
                st.text(matrix19)
                aligner.substitution_matrix = matrix19
                score19 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score19))
            elif option == "PAM30":
                matrix20 = substitution_matrices.load("PAM30")
                st.text(matrix20)
                aligner.substitution_matrix = matrix20
                score20 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score20))
            elif option == "PAM70":
                matrix21 = substitution_matrices.load("PAM70")
                st.text(matrix21)
                aligner.substitution_matrix = matrix21
                score21 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score21))
            elif option == "RAO":
                matrix22 = substitution_matrices.load("RAO")
                st.text(matrix22)
                aligner.substitution_matrix = matrix22
                score22 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score22))
            elif option == "RISLER":
                matrix23 = substitution_matrices.load("RISLER")
                st.text(matrix23)
                aligner.substitution_matrix = matrix23
                score23 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score23))
            elif option == "SCHNEIDER":
                matrix24 = substitution_matrices.load("SCHNEIDER")
                st.text(matrix24)
                aligner.substitution_matrix = matrix24
                score24 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score24))
            elif option == "STR":
                matrix25 = substitution_matrices.load("STR")
                st.text(matrix25)
                aligner.substitution_matrix = matrix25
                score25 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Score of Sequence** = {0} ".format(score25))
            
            #aligner.substitution_matrix = matrix
            #score = aligner.score("ACDQ", "ACDQ")
            #st.write(score)
        except:
            st.text("")

    elif select_box == "Target Protein Analysis":
        try:
            st.text("")
            arama = st.sidebar.text_input("Write the protein that you want to research", "coronavirus")
            target = new_client.target
            target_query = target.search(arama) #Your target protein name that you want to search
            targets = pd.DataFrame.from_dict(target_query)
            st.write("**Chembl Data**:")
            st.write(targets)
            st.text("")
            
            hedef_protein = st.number_input("Enter the number to the left of the Single Protein you wish to further search for.", min_value=0 ,value=4, format="%d")
            selected_target = targets.target_chembl_id[hedef_protein]
            st.text("")
            st.write("**ChEMBL ID** of the protein that you choose: {0}".format(selected_target))
            activity = new_client.activity
            res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
            df = pd.DataFrame.from_dict(res)
            
            st.text("")
        except:
            st.warning("Invalid Value")
        

        if hedef_protein is not None:
            

            if st.checkbox("Click to see the IC50 value of the protein you selected."):
                if df.empty:
                    st.warning("Please choose a **Single Protein**!")
                else:
                        
                    st.write(df)
                    st.text("")
                
                    st.markdown(download_button(df, 'IC50_Data.csv', 'Download the CSV File'), unsafe_allow_html=True)
            
            
            st.text("")
            st.text("")
            st.markdown("<h3 style='text-align: center; color: red;'>Click the button for calculate the Molecular Activities</h3>", unsafe_allow_html=True)
            st.text("")
            
            df2 = df[df.standard_value.notna()]
            bioactivity_class = []

            mol_cid = []
            canonical_smiles = []
            standard_value = []
            
            for unit in df2.standard_value:
                if float(unit) >= 10000:
                    bioactivity_class.append("inactive")
                elif float(unit) <= 1000:
                    bioactivity_class.append("active")
                else:
                    bioactivity_class.append("intermediate")
            

            for i in df2.molecule_chembl_id:

                mol_cid.append(i)
            
            for i in df2.canonical_smiles:

                canonical_smiles.append(i)
                
            for i in df2.standard_value:

                standard_value.append(i)
                    
                    

            data_tuples = list(zip(mol_cid, canonical_smiles, standard_value, bioactivity_class))
            df3 = pd.DataFrame( data_tuples,  columns=['molecule_chembl_id', 'canonical_smiles', 'standard_value','bioactivity_class' ])
            st.text("")
            if df.empty:
                st.warning("Please choose a **Single Protein**!")
            else:
                if st.checkbox("Calculate the Molecular Activity"):
                    st.text("")
                    st.text("")
                    st.write(df3)
                    st.text("")
                    st.markdown(download_button(df3, 'General_Data.csv', 'Download the CSV File'), unsafe_allow_html=True)

                    st.text("")
                    if st.selectbox("Just show the active ones",("Active","")):
                        active_data = (df3.loc[df3['bioactivity_class'] == "active"])
                        st.write(active_data)
                        st.text("")
                        st.markdown(download_button(active_data, 'Active_Data.csv', 'Download the CSV File'), unsafe_allow_html=True)

            
                st.text("")
                st.text("")
                st.markdown("<h3 style='text-align: center; color: red;'>Click the button for calculate Lipinski Descriptors</h3>", unsafe_allow_html=True)
                st.text("")
                
                button_sent = st.checkbox("Lipinski Descriptors")
                
                if button_sent:
                    session_state.button_sent = True

                if session_state.button_sent:
                    st.subheader("Lipinski Data:")
                    st.write("**MW** = Molecular Weight")
                    st.write("**LogP** = Molecule Solubility")
                    st.write("**NumHDonors** = Hydrogen Bounds Donors")
                    st.write("**NumHAcceptors** = Hydrogen Bounds Acceptors")
        
                    exploratory_data = df3
                    df_lipinski = lipinski(exploratory_data.canonical_smiles)
                    #st.write(df_lipinski)
                    df_combined = pd.concat([exploratory_data,df_lipinski], axis=1)
                    st.subheader("Combined Data:")
                    st.write(df_combined)
                    st.markdown(download_button(df_combined, 'Combined_Data.csv', 'Download the CSV File'), unsafe_allow_html=True)
                    df_norm = norm_value(df_combined)
                    #st.write(df_norm)
                    df_final = pIC50(df_norm)
                    st.subheader("Dataset with IC50 converted to pIC50:")
                    st.write(df_final)
                    st.markdown(download_button(df_final, 'pIC50_Data.csv', 'Download the CSV File'), unsafe_allow_html=True)
                    df_class = df_final
                    #İntermediate ekleyen kod
                    #df_class = df_final[df_final.bioactivity_class != "intermediate"]

                

                    def mannwhitney(descriptor, verbose=False):

                        # seed the random number generator
                        seed(1)

                        # actives and inactives
                        selection = [descriptor, 'bioactivity_class']
                        df = df_class[selection]
                        active = df[df.bioactivity_class == 'active']
                        active = active[descriptor]

                        selection = [descriptor, 'bioactivity_class']
                        df = df_class[selection]
                        inactive = df[df.bioactivity_class == 'inactive']
                        inactive = inactive[descriptor]

                        selection = [descriptor, "bioactivity_class"]
                        df = df_class[selection]
                        intermediate = df[df.bioactivity_class == "intermediate"]
                        intermediate = intermediate[descriptor]

                        # compare samples
                        stat, p = mannwhitneyu(active, inactive)
                        #print('Statistics=%.3f, p=%.3f' % (stat, p))

                        # interpret
                        alpha = 0.05
                        if p > alpha:
                            interpretation = 'Same distribution (fail to reject H0)'
                        else:
                            interpretation = 'Different distribution (reject H0)'
                        
                        results = pd.DataFrame({'Descriptor':descriptor,
                                                'Statistics':stat,
                                                'p':p,
                                                'alpha':alpha,
                                                'Interpretation':interpretation}, index=[0])
                        filename = 'mannwhitneyu_' + descriptor + '.csv'
                        results.to_csv(filename)

                        return results

                    st.text("")
                    st.text("")
                    session_state.grafik = st.checkbox("Active/Inactive Molecules Graph")
                    session_state.mw = st.checkbox("Molecular Weight/Solubility Graph")
                    session_state.pic50 = st.checkbox("pIC50/Molecular Activity Graph")
                    session_state.logp = st.checkbox("Solubility/Molecular Activity Graph")
                    session_state.donors = st.checkbox("Hydrogen Bounds Donors/Molecular Activity Graph")
                    session_state.acceptors = st.checkbox("Hydrogen Bounds Acceptors/Molecular Activity Graph")
                    session_state.molecul = st.checkbox("Molecular Weight Graph")

                    

                    if session_state.grafik:
                        st.write("**********************************")
                        st.text("")
                        st.subheader("**Active/Inactive Molecules Graph**")

                        plt.figure(figsize=(5.5, 5.5))

                        sns.countplot(x='bioactivity_class', data=df_class, edgecolor='black')

                        plt.xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
                        plt.ylabel('Frequency', fontsize=14, fontweight='bold')
                        
                        st.pyplot()
                        #st.markdown(get_table_download_link(veri), unsafe_allow_html=True)
                        
                        #Buralara PDF indirici eklenecek

                    if session_state.mw:
                        st.write("**********************************")
                        st.text("")
                        st.subheader("**Molecular Weight/Solubility Graph**")

                        plt.figure(figsize=(5.5, 5.5))
                        sns.scatterplot(x='MW', y='LogP', data=df_class, hue='bioactivity_class', size='pIC50', edgecolor='black', alpha=0.7)

                        plt.xlabel('Molecular Weight', fontsize=14, fontweight='bold')
                        plt.ylabel('LogP', fontsize=14, fontweight='bold')
                        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
                        st.pyplot()
                        
                        #Buralara PDF indirici eklenecek
                        st.write("**Mann-Whitney U Test Data**:")
                        st.write(mannwhitney("MW"))

                    if session_state.pic50:
                        st.write("**********************************")
                        st.text("")
                        st.subheader("**pIC50/Molecular Activity Graph**")

                        plt.figure(figsize=(5.5, 5.5))

                        sns.boxplot(x = 'bioactivity_class', y = 'pIC50', data = df_class)

                        plt.xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
                        plt.ylabel('pIC50', fontsize=14, fontweight='bold')
                        st.pyplot()
                        #Buralara PDF indirici eklenecek

                        st.write("**Mann-Whitney U Test Data**:")
                        st.write(mannwhitney("pIC50"))
                    
                    if session_state.logp:
                        st.write("**********************************")
                        st.text("")
                        st.subheader("**Solubility/Molecular Activity Graph**")

                        plt.figure(figsize=(5.5, 5.5))

                        sns.boxplot(x = 'bioactivity_class', y = 'LogP', data = df_class)

                        plt.xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
                        plt.ylabel('LogP', fontsize=14, fontweight='bold')
                        st.pyplot()
                        #Buralara PDF indirici eklenecek

                        st.write("**Mann-Whitney U Test Data**:")
                        st.write(mannwhitney("LogP"))
                    
                    if session_state.donors:
                        st.write("**********************************")
                        st.text("")
                        st.subheader("**Hydrogen Bounds Donors/Molecular Activity Graph**")

                        plt.figure(figsize=(5.5, 5.5))

                        sns.boxplot(x = 'bioactivity_class', y = 'NumHDonors', data = df_class)

                        plt.xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
                        plt.ylabel('NumHDonors', fontsize=14, fontweight='bold')
                        st.pyplot()
                        #Buralara PDF indirici eklenecek

                        st.write("**Mann-Whitney U Test Data**:")
                        st.write(mannwhitney("NumHDonors"))

                    if session_state.acceptors:
                        st.write("**********************************")
                        st.text("")
                        st.subheader("**Hydrogen Bounds Acceptors/Molecular Activity Graph**")

                        plt.figure(figsize=(5.5, 5.5))

                        sns.boxplot(x = 'bioactivity_class', y = 'NumHAcceptors', data = df_class)

                        plt.xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
                        plt.ylabel('NumHAcceptors', fontsize=14, fontweight='bold')
                        st.pyplot()
                        #Buralara PDF indirici eklenecek

                        st.write("**Mann-Whitney U Test Data**:")
                        st.write(mannwhitney("NumHAcceptors"))

                    if session_state.molecul:
                        st.write("**********************************")
                        st.text("")
                        st.subheader("**Molecular Weight Graph**")
                        plt.figure(figsize=(5.5, 5.5))

                        sns.boxplot(x = 'bioactivity_class', y = 'MW', data = df_class)

                        plt.xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
                        plt.ylabel('Molecular Weight', fontsize=14, fontweight='bold')
                        st.pyplot()


                
            
            

    elif select_box == "Find SNP and INDEL":

        align_dosya = st.sidebar.file_uploader("Upload your FASTA file", type=["FASTA","fa"])
        

        y=0
        alignment = SeqIO.read(align_dosya, "fasta")
        for r in range(0,len(alignment[1].seq)):
            if alignment[0,r] != alignment[1,r]:
                if alignment[0,r] != "-" and alignment[1,r] != "-":
                    y=y+1
                    st.write(r, alignment[0,r], alignment[1,r], y)
                else:
                    y=0



    

        



                    



if __name__ == "__main__":
    main()