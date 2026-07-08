"""
TransE.py -- EMBEDDING EXPORT ONLY (not an evaluation script).

⚠️  The previous version of this file evaluated TransE with `training=tf,
    testing=tf` -- i.e. it tested on the exact triples it trained on, which makes
    any reported metric meaningless (PROJECT_REPORT.md §6.2). That evaluation has
    been removed.

    For HONEST link-prediction evaluation (held-out split, filtered ranking,
    baselines, multiple KGE models) use:  scripts/kge_benchmark.py

This script trains TransE on the FULL graph purely to export node embeddings for
downstream/unsupervised use (visualization, clustering). Training on all edges is
fine for that purpose; it is NOT a performance claim.
"""
import argparse
import pandas as pd
import torch
from pykeen.triples import TriplesFactory
from pykeen.models import TransE
from pykeen.training import SLCWATrainingLoop


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", default="data/processed/edges_clean_integrated.csv")
    ap.add_argument("--out", default="data/processed/transe_embeddings_full.csv")
    ap.add_argument("--dim", type=int, default=128)
    ap.add_argument("--epochs", type=int, default=100)
    ap.add_argument("--batch", type=int, default=16384)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    df = pd.read_csv(args.edges)[["source_id", "relation", "target_id"]]
    df = df[df["source_id"] != df["target_id"]]
    tf = TriplesFactory.from_labeled_triples(df.to_numpy())
    print(f"Nodes: {tf.num_entities}  Relations: {tf.num_relations}  Triples: {tf.num_triples}")

    device = "cuda" if torch.cuda.is_available() else "cpu"
    torch.manual_seed(args.seed)
    model = TransE(triples_factory=tf, embedding_dim=args.dim, random_seed=args.seed).to(device)
    loop = SLCWATrainingLoop(model=model, triples_factory=tf)
    loop.train(triples_factory=tf, num_epochs=args.epochs, batch_size=args.batch, use_tqdm=True)

    emb = model.entity_representations[0](indices=None).cpu().detach().numpy()
    out = pd.DataFrame(emb)
    out.insert(0, "node", list(tf.entity_to_id.keys()))
    out.to_csv(args.out, index=False)
    print(f"✅ Saved TransE embeddings for {tf.num_entities} nodes -> {args.out}")
    print("ℹ️  For evaluation/metrics use scripts/kge_benchmark.py (held-out split).")


if __name__ == "__main__":
    main()
