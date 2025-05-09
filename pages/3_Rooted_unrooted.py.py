import streamlit as st

# --- Page config for better appearance ---
st.set_page_config(
    page_title="EvoMetrics | Tree Enumeration",
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
            background: rgba(255,255,255,0.93);
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
        .stMetric {
            background: #e0f2fe;
            border-radius: 12px;
            padding: 0.5rem 1rem;
            margin-bottom: 1rem;
        }
    </style>
""", unsafe_allow_html=True)
import matplotlib.pyplot as plt
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
from Bio import Phylo
import math
import random

# --- Helper functions ---

def count_rooted_trees(n):
    # Number of rooted binary trees: (2n-3)! / [2^(n-2) * (n-2)!]
    if n < 2:
        return 0
    return math.factorial(2 * n - 3) // (2 ** (n - 2) * math.factorial(n - 2))

def count_unrooted_trees(n):
    # Number of unrooted binary trees: (2n-5)! / [2^(n-3) * (n-3)!]
    if n < 3:
        return 0
    return math.factorial(2 * n - 5) // (2 ** (n - 3) * math.factorial(n - 3))

def random_distance_matrix(names):
    # Generate a lower triangle random distance matrix for demonstration
    n = len(names)
    matrix = []
    for i in range(n):
        row = []
        for j in range(i):
            row.append(round(random.uniform(0.1, 1.0), 3))  # random distances
        matrix.append(row)
    return _DistanceMatrix(names, matrix)

def plot_tree(tree, rooted=True):
    fig = plt.figure(figsize=(6, 4))
    Phylo.draw(tree, do_show=False, rooted=rooted)
    st.pyplot(fig)
    plt.close(fig)

# --- Streamlit App ---

st.title("Phylogenetic Tree Constructor")
st.write("Visualize rooted and unrooted phylogenetic trees and calculate the number of possible trees for a given number of nodes.")

n_nodes = st.number_input("Number of taxa (nodes/leaves)", min_value=3, max_value=12, value=5, step=1)
taxa = [f"Taxon_{i+1}" for i in range(n_nodes)]

st.subheader("Number of Possible Trees")
col1, col2 = st.columns(2)
with col1:
    st.metric("Rooted Trees", f"{count_rooted_trees(n_nodes):,}")
with col2:
    st.metric("Unrooted Trees", f"{count_unrooted_trees(n_nodes):,}")

st.caption("Powered by Streamlit & Biopython")
