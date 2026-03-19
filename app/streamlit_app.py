"""Streamlit app for MMHC Bayesian network structure learning.

Run with:
    streamlit run app/streamlit_app.py
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st
from pyvis.network import Network

from mmhc import MMHC, MMHCConfig, MMHCResult, make_rainy, make_student
from mmhc.graph_utils import adjacency_to_networkx, is_dag

# ---------------------------------------------------------------------------
# Ground truth definitions for built-in datasets
# ---------------------------------------------------------------------------
GROUND_TRUTHS: dict[str, dict[str, object]] = {
    "Student network (built-in)": {
        "columns": ["difficulty", "intelligence", "SAT", "grade", "letter"],
        "edges": [
            ("difficulty", "grade"),
            ("intelligence", "grade"),
            ("intelligence", "SAT"),
            ("grade", "letter"),
        ],
        "description": (
            "**Ground truth:** difficulty -> grade, intelligence -> grade, "
            "intelligence -> SAT, grade -> letter"
        ),
    },
    "Rainy network (built-in)": {
        "columns": ["sprinkler", "rain", "grassWet"],
        "edges": [
            ("rain", "sprinkler"),
            ("rain", "grassWet"),
            ("sprinkler", "grassWet"),
        ],
        "description": (
            "**Ground truth:** rain -> sprinkler, rain -> grassWet, "
            "sprinkler -> grassWet"
        ),
    },
}


def _build_gt_adjacency(gt: dict[str, object]) -> np.ndarray:
    """Build ground truth adjacency matrix from edge list."""
    columns: list[str] = gt["columns"]  # type: ignore[assignment]
    edges: list[tuple[str, str]] = gt["edges"]  # type: ignore[assignment]
    n = len(columns)
    adj = np.zeros((n, n), dtype=np.int64)
    for src, dst in edges:
        adj[columns.index(src), columns.index(dst)] = 1
    return adj


def _render_pyvis(adj: np.ndarray, labels: list[str], height: str = "450px") -> str:
    """Render a pyvis graph and return the HTML content."""
    g = adjacency_to_networkx(adj, labels)
    net = Network(
        height=height,
        width="100%",
        directed=True,
        notebook=False,
        cdn_resources="remote",
    )
    net.from_nx(g)

    for node in net.nodes:
        node["size"] = 35
        node["color"] = {
            "background": "#4A90D9",
            "border": "#2C6FB5",
            "highlight": {"background": "#6DB3F8", "border": "#2C6FB5"},
        }
        node["font"] = {"size": 18, "color": "#FFFFFF", "bold": {"color": "#FFFFFF"}}
        node["shape"] = "box"
        node["margin"] = 12
        node["borderWidth"] = 2

    for edge in net.edges:
        edge["color"] = {"color": "#888888", "highlight": "#4A90D9"}
        edge["width"] = 2.5
        edge["arrows"] = {"to": {"enabled": True, "scaleFactor": 1.2}}
        edge["smooth"] = {"type": "curvedCW", "roundness": 0.15}

    net.set_options("""{
        "physics": {
            "solver": "forceAtlas2Based",
            "forceAtlas2Based": {
                "gravitationalConstant": -120,
                "centralGravity": 0.01,
                "springLength": 250,
                "springConstant": 0.04
            },
            "stabilization": {"iterations": 200}
        },
        "interaction": {
            "hover": true,
            "tooltipDelay": 100
        }
    }""")

    with tempfile.NamedTemporaryFile(suffix=".html", delete=False, mode="w") as f:
        net.save_graph(f.name)
        return Path(f.name).read_text()


def _compare_edges(
    learned_adj: np.ndarray,
    gt_adj: np.ndarray,
    labels: list[str],
) -> pd.DataFrame:
    """Compare learned edges against ground truth."""
    n = len(labels)
    rows = []

    # Check ground truth edges
    for i in range(n):
        for j in range(n):
            if gt_adj[i, j] == 1:
                src, dst = labels[i], labels[j]
                if learned_adj[i, j] == 1:
                    rows.append({"Edge": f"{src} -> {dst}", "Status": "Correct"})
                elif learned_adj[j, i] == 1:
                    rows.append({"Edge": f"{src} -> {dst}", "Status": "Reversed"})
                else:
                    rows.append({"Edge": f"{src} -> {dst}", "Status": "Missing"})

    # Check for extra edges not in ground truth
    for i in range(n):
        for j in range(n):
            if learned_adj[i, j] == 1 and gt_adj[i, j] == 0 and gt_adj[j, i] == 0:
                rows.append(
                    {"Edge": f"{labels[i]} -> {labels[j]}", "Status": "Extra"}
                )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="MMHC — Bayesian Network Learner",
    page_icon="🔗",
    layout="wide",
)

st.title("MMHC — Bayesian Network Structure Learning")
st.markdown(
    "Upload a dataset (Excel or CSV) with **discrete/categorical columns**, "
    "configure the algorithm, and discover the causal structure."
)

# ---------------------------------------------------------------------------
# Sidebar — data source
# ---------------------------------------------------------------------------
st.sidebar.header("Data Source")
data_source = st.sidebar.radio(
    "Choose data source",
    ["Upload file", "Student network (built-in)", "Rainy network (built-in)"],
)

data: pd.DataFrame | None = None

if data_source == "Upload file":
    uploaded = st.sidebar.file_uploader(
        "Upload Excel (.xlsx) or CSV (.csv)",
        type=["xlsx", "xls", "csv"],
    )
    if uploaded is not None:
        data = (
            pd.read_csv(uploaded)
            if uploaded.name.endswith(".csv")
            else pd.read_excel(uploaded)
        )
elif data_source == "Student network (built-in)":
    n_samples = st.sidebar.slider("Number of samples", 500, 20000, 5000, step=500)
    data = make_student(n_samples, random_state=42)
else:
    n_samples = st.sidebar.slider("Number of samples", 500, 20000, 5000, step=500)
    data = make_rainy(n_samples, random_state=42)

# ---------------------------------------------------------------------------
# Sidebar — algorithm config
# ---------------------------------------------------------------------------
st.sidebar.header("Algorithm Parameters")
alpha = st.sidebar.slider(
    "alpha (significance level)",
    min_value=0.001,
    max_value=0.20,
    value=0.05,
    step=0.005,
    help="Lower = stricter independence test = fewer edges",
)
eta = st.sidebar.slider(
    "eta (BDeu equivalent sample size)",
    min_value=0.1,
    max_value=20.0,
    value=1.0,
    step=0.1,
    help="Higher = stronger prior favouring simpler graphs",
)
max_iterations = st.sidebar.slider(
    "Max iterations",
    min_value=10,
    max_value=500,
    value=100,
    step=10,
)
early_stop = st.sidebar.slider(
    "Early stop rounds",
    min_value=1,
    max_value=20,
    value=5,
)
random_seed = st.sidebar.number_input(
    "Random seed",
    min_value=0,
    max_value=999999,
    value=42,
    help="For reproducible results",
)

# ---------------------------------------------------------------------------
# Main area — data preview
# ---------------------------------------------------------------------------
if data is not None:
    st.subheader("Data Preview")
    col1, col2 = st.columns([3, 1])
    with col1:
        st.dataframe(data.head(50), use_container_width=True)
    with col2:
        st.metric("Rows", data.shape[0])
        st.metric("Columns", data.shape[1])
        st.markdown("**Column types:**")
        for col_name in data.columns:
            n_unique = data[col_name].nunique()
            st.text(f"  {col_name}: {n_unique} values")

    # Show ground truth for built-in datasets
    if data_source in GROUND_TRUTHS:
        gt = GROUND_TRUTHS[data_source]
        st.info(gt["description"])

    # -----------------------------------------------------------------------
    # Run algorithm
    # -----------------------------------------------------------------------
    run_button = st.button("Run MMHC", type="primary", use_container_width=True)

    if run_button:
        config = MMHCConfig(
            alpha=alpha,
            eta=eta,
            max_iterations=max_iterations,
            early_stop_rounds=early_stop,
            random_seed=int(random_seed),
        )

        with st.spinner("Running MMHC algorithm..."):
            model = MMHC(config=config)
            result: MMHCResult = model.fit(data)

        st.session_state["result"] = result
        st.session_state["data_source"] = data_source

    # -----------------------------------------------------------------------
    # Display results
    # -----------------------------------------------------------------------
    if "result" in st.session_state:
        result = st.session_state["result"]
        st.divider()
        st.subheader("Results")

        # Metrics row
        m1, m2, m3, m4 = st.columns(4)
        n_edges = int(result.adjacency_matrix.sum())
        m1.metric("Edges found", n_edges)
        m2.metric("BDeu Score", f"{result.score:.2f}")
        m3.metric("Iterations", result.n_iterations)
        m4.metric("Converged", "Yes" if result.converged else "No")

        # DAG validation
        if is_dag(result.adjacency_matrix):
            st.success("Valid DAG (no cycles)")
        else:
            st.error("Warning: Result contains cycles")

        # -------------------------------------------------------------------
        # Ground truth comparison (built-in datasets only)
        # -------------------------------------------------------------------
        active_source = st.session_state.get("data_source", data_source)
        has_gt = active_source in GROUND_TRUTHS

        if has_gt:
            gt = GROUND_TRUTHS[active_source]
            gt_adj = _build_gt_adjacency(gt)

            st.subheader("Ground Truth vs Learned")

            gt_col, learned_col = st.columns(2)
            with gt_col:
                st.markdown("**Ground Truth**")
                gt_html = _render_pyvis(gt_adj, gt["columns"], height="400px")  # type: ignore[arg-type]
                st.components.v1.html(gt_html, height=420, scrolling=False)

            with learned_col:
                st.markdown("**Learned by MMHC**")
                learned_html = _render_pyvis(
                    result.adjacency_matrix, result.column_names, height="400px"
                )
                st.components.v1.html(learned_html, height=420, scrolling=False)

            # Edge comparison table
            comparison = _compare_edges(
                result.adjacency_matrix, gt_adj, result.column_names
            )
            if not comparison.empty:
                st.markdown("**Edge Comparison**")

                def _color_status(val: str) -> str:
                    colors = {
                        "Correct": "background-color: #2ECC71; color: white",
                        "Reversed": "background-color: #F39C12; color: white",
                        "Missing": "background-color: #E74C3C; color: white",
                        "Extra": "background-color: #9B59B6; color: white",
                    }
                    return colors.get(val, "")

                st.dataframe(
                    comparison.style.map(_color_status, subset=["Status"]),
                    use_container_width=True,
                    hide_index=True,
                )

                n_correct = int((comparison["Status"] == "Correct").sum())
                n_total_gt = int(gt_adj.sum())
                n_reversed = int((comparison["Status"] == "Reversed").sum())

                if n_reversed > 0:
                    st.caption(
                        f"**{n_correct}/{n_total_gt}** edges direction-correct, "
                        f"**{n_reversed}** reversed. "
                        "Reversed edges are expected — BDeu scoring cannot distinguish "
                        "between Markov-equivalent DAGs from observational data alone."
                    )

        else:
            # -------------------------------------------------------------------
            # Full-width visualization for uploaded data
            # -------------------------------------------------------------------
            st.subheader("Network Visualization")
            html_content = _render_pyvis(
                result.adjacency_matrix, result.column_names
            )
            st.components.v1.html(html_content, height=520, scrolling=False)

        # -------------------------------------------------------------------
        # Edge list
        # -------------------------------------------------------------------
        st.subheader("Edge List")
        edges = []
        for i, src in enumerate(result.column_names):
            for j, dst in enumerate(result.column_names):
                if result.adjacency_matrix[i, j] == 1:
                    edges.append({"Source": src, "Target": dst})

        if edges:
            edge_df = pd.DataFrame(edges)
            st.dataframe(edge_df, use_container_width=True, hide_index=True)
        else:
            st.info("No edges found. Try increasing alpha or using more data.")

        # -------------------------------------------------------------------
        # Adjacency matrix
        # -------------------------------------------------------------------
        st.subheader("Adjacency Matrix")
        adj_df = pd.DataFrame(
            result.adjacency_matrix,
            index=result.column_names,
            columns=result.column_names,
        )
        st.dataframe(
            adj_df.style.map(
                lambda v: "background-color: #4A90D9; color: white" if v == 1 else ""
            ),
            use_container_width=True,
        )

        # -------------------------------------------------------------------
        # Per-node scores
        # -------------------------------------------------------------------
        st.subheader("Per-Node BDeu Scores")
        score_df = pd.DataFrame({
            "Node": result.column_names,
            "BDeu Score": result.node_scores,
            "Parents": [
                ", ".join(
                    result.column_names[p]
                    for p in range(len(result.column_names))
                    if result.adjacency_matrix[p, i] == 1
                ) or "(none)"
                for i in range(len(result.column_names))
            ],
        })
        st.dataframe(score_df, use_container_width=True, hide_index=True)

        # -------------------------------------------------------------------
        # Downloads
        # -------------------------------------------------------------------
        st.subheader("Download Results")
        dl1, dl2 = st.columns(2)
        with dl1:
            csv_adj = adj_df.to_csv()
            st.download_button(
                "Download adjacency matrix (CSV)",
                csv_adj,
                file_name="adjacency_matrix.csv",
                mime="text/csv",
            )
        with dl2:
            if edges:
                csv_edges = edge_df.to_csv(index=False)
                st.download_button(
                    "Download edge list (CSV)",
                    csv_edges,
                    file_name="edge_list.csv",
                    mime="text/csv",
                )

else:
    st.info("Upload a file or select a built-in dataset to get started.")
