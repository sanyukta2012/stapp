import streamlit as st

# --- Page config for better appearance ---
st.set_page_config(
    page_title="EvoMetrics | Phylogenetic Tree Constructor",
    page_icon="🧬",
    layout="wide",
)

# --- Custom CSS for background and style ---
st.markdown("""
    <style>
        body {
            background: linear-gradient(120deg, #e0f2fe 0%, #f0fdfa 100%);
            background-image: url('https://images.unsplash.com/photo-1465101046530-73398c7f28ca?auto=format&fit=crop&w=1200&q=80');
            background-size: cover;
            background-attachment: fixed;
        }
        .stApp {
            background: rgba(255,255,255,0.92);
            border-radius: 18px;
            padding: 2.5rem;
            box-shadow: 0 6px 36px rgba(0,0,0,0.08);
        }
        h1, h2, h3, h4 {
            color: #0e7490 !important;
            letter-spacing: 0.5px;
        }
        .stSidebar {
            background: rgba(236,254,255,0.85);
        }
        .css-1v0mbdj p, .css-1v0mbdj ul, .css-1v0mbdj li {
            font-size: 1.08rem;
        }
        .stAlert {
            border-radius: 12px;
            background: #e0f2fe;
            color: #0369a1;
        }
        .stExpander {
            border-radius: 12px;
            background: #f0fdfa;
        }
    </style>
""", unsafe_allow_html=True)
from Bio import SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from io import StringIO
import matplotlib.pyplot as plt
import pandas as pd

st.title("Phylogenetic Tree Constructor")

st.write("""
Upload one or more FASTA files **or** paste your sequences in FASTA format below.  
Select the method (UPGMA or Neighbor Joining) to construct your phylogenetic tree.
""")

method = st.radio(
    "Select Tree Construction Method:",
    ("UPGMA (Molecular Clock)", "Neighbor Joining (No Molecular Clock)")
)

tab1, tab2 = st.tabs(["Upload FASTA Files", "Paste Sequences"])

records = []

with tab1:
    uploaded_files = st.file_uploader("Upload FASTA files", type=["fasta", "fa", "txt"], accept_multiple_files=True)
    if uploaded_files:
        for uploaded_file in uploaded_files:
            file_content = uploaded_file.read().decode("utf-8")
            file_io = StringIO(file_content)
            file_records = list(SeqIO.parse(file_io, "fasta"))
            records.extend(file_records)

with tab2:
    example_fasta = """>seq1
ATGCTAGCTAGCTACGATCG
>seq2
ATGCTAGCTAGCTACGATGG
>seq3
ATGCTAGCTAGATACGATCG
"""
    user_input = st.text_area("Paste FASTA sequences here:", value=example_fasta, height=200)
    if st.button("Build Tree from Pasted Sequences"):
        fasta_io = StringIO(user_input)
        try:
            pasted_records = list(SeqIO.parse(fasta_io, "fasta"))
            records.extend(pasted_records)
        except Exception as e:
            st.error(f"Error parsing pasted sequences: {str(e)}")

# Remove duplicate IDs if any
unique_records = []
seen_ids = set()
for rec in records:
    if rec.id not in seen_ids:
        unique_records.append(rec)
        seen_ids.add(rec.id)
records = unique_records

if len(records) >= 2:
    # Align sequences by padding with gaps to the max length
    max_len = max(len(rec.seq) for rec in records)
    for rec in records:
        seq = str(rec.seq).ljust(max_len, "-")
        rec.seq = rec.seq.__class__(seq)
    alignment = MultipleSeqAlignment(records)

    # Distance calculation
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()

    # Tree construction based on user selection
    if method.startswith("UPGMA"):
        tree = constructor.upgma(dm)
        tree_type = "UPGMA Phylogenetic Tree"
    else:
        tree = constructor.nj(dm)
        tree_type = "Neighbor Joining Phylogenetic Tree"

    # Display distance matrix as a DataFrame
    st.subheader("Distance Matrix")
    df = pd.DataFrame(dm.matrix, index=dm.names, columns=dm.names)
    st.dataframe(df)

    # Draw tree with branch distances
    st.subheader(tree_type)
    fig, ax = plt.subplots(figsize=(8, 5))
    Phylo.draw(
        tree,
        axes=ax,
        branch_labels=lambda c: f"{c.branch_length:.5f}" if c.branch_length is not None else None,
        do_show=False
    )
    st.pyplot(fig)

    # Show Newick string
    st.subheader("Newick Format")
    newick_io = StringIO()
    Phylo.write(tree, newick_io, "newick")
    st.code(newick_io.getvalue())

elif len(records) == 1:
    st.warning("Only one sequence found. Please provide at least two sequences.")

else:
    st.info("Please upload FASTA files or paste at least two sequences.")

st.markdown("---")
st.markdown("**Developed with [Biopython](https://biopython.org/) and [Streamlit](https://streamlit.io/)**")

