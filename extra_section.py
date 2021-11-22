"""
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
"""
