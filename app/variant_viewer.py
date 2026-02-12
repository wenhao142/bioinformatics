"""
Simple Postgres-backed variant browser for local use.
- Query by chromosome/start/end and optional sample_id.
- Shows results table and lets user download CSV.
"""

import os
import io
 streamlit as st
import pandas as pd
import psycopg2
import psycopg2.extras


def get_conn():
    return psycopg2.connect(
        host=os.getenv("POSTGRES_HOST", "postgres"),
        port=os.getenv("POSTGRES_PORT", "5432"),
        dbname=os.getenv("POSTGRES_DB", "omics"),
        user=os.getenv("POSTGRES_USER", "postgres"),
        password=os.getenv("POSTGRES_PASSWORD", "postgres"),
    )


def run_query(chr_, start, end, sample_id):
    sql = """
        SELECT chr, pos, ref, alt, qual, filter_status, sample_id, af, ac, dp
        FROM variants
        WHERE chr = %s AND pos BETWEEN %s AND %s
        {sample_clause}
        ORDER BY pos
        LIMIT 1000;
    """
    clause = ""
    params = [chr_, start, end]
    if sample_id:
        clause = "AND sample_id = %s"
        params.append(sample_id)
    sql = sql.format(sample_clause=clause)
    with get_conn() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
            cur.execute(sql, params)
            rows = cur.fetchall()
    return pd.DataFrame(rows, columns=["chr", "pos", "ref", "alt", "qual", "filter", "sample_id", "af", "ac", "dp"])


def main():
    st.title("Local Variant Viewer")
    st.caption("Query variants stored in Postgres (loaded by auto_loader/auto_caller).")

    chr_ = st.text_input("Chromosome (e.g., chr1, 1)", value="chr1")
    col1, col2 = st.columns(2)
    with col1:
        start = st.number_input("Start", min_value=1, value=1_000_000)
    with col2:
        end = st.number_input("End", min_value=1, value=1_010_000)
    sample_id = st.text_input("Sample ID (optional)", value="")

    if st.button("Run query", type="primary"):
        try:
            df = run_query(chr_, start, end, sample_id.strip() or None)
        except Exception as exc:
            st.error(f"Query failed: {exc}")
        else:
            st.success(f"Returned {len(df)} rows (max 1000).")
            st.dataframe(df, use_container_width=True)
            if not df.empty:
                csv = df.to_csv(index=False).encode("utf-8")
                st.download_button("Download CSV", data=csv, file_name="variants.csv", mime="text/csv")


if __name__ == "__main__":
    main()
