import pandas as pd
from pykeen.triples import TriplesFactory
from pykeen.pipeline import pipeline
from pykeen.predict import predict_target

# =========================
# STEP 1: LOAD CSV CLEANLY
# =========================

# Read CSV normally (keep header)
df = pd.read_csv("data/processed/transetableuseforpykeen.csv")

# Ensure correct columns
df = df[['head', 'relation', 'tail']]

# Convert to numpy array
triples = df.values

print("Triples shape:", triples.shape)

# =========================
# STEP 2: CREATE TRIPLES FACTORY
# =========================

tf = TriplesFactory.from_labeled_triples(triples)

# =========================
# STEP 3: TRAIN TransE MODEL
# =========================

result = pipeline(
    training=tf,
    testing=tf,  # REQUIRED (fixes your error)
    model='TransE',
    epochs=5,    # keep small for speed
    random_seed=42,
)

model = result.model

# =========================
# STEP 4: RUN LINK PREDICTION
# =========================

relation = "CONSENSUS_SIMILAR"

all_predictions = []

# Loop through UNIQUE diseases (tails)
unique_tails = df['tail'].unique()

print("Running predictions for", len(unique_tails), "targets...")

for tail in unique_tails:
    try:
        preds = predict_target(
            model=model,
            relation=relation,
            tail=tail,
            triples_factory=tf,
        )

        pred_df = preds.df

        # Add which tail we predicted for
        pred_df['target_tail'] = tail

        all_predictions.append(pred_df)

    except Exception as e:
        print("Error with tail:", tail, "|", e)

# =========================
# STEP 5: COMBINE RESULTS
# =========================

if len(all_predictions) > 0:
    final_df = pd.concat(all_predictions, ignore_index=True)

    print("Total predictions:", final_df.shape)

    # =========================
    # STEP 6: REMOVE EXISTING EDGES
    # =========================

    existing = set(zip(df['head'], df['relation'], df['tail']))

    final_df = final_df[
        ~final_df.apply(
            lambda row: (row['head_label'], relation, row['target_tail']) in existing,
            axis=1
        )
    ]

    print("After removing known edges:", final_df.shape)

    # =========================
    # STEP 7: SORT + SAVE
    # =========================

    final_df = final_df.sort_values(by="score", ascending=False)

    print(final_df.head(20))

    final_df.to_csv("data/processed/transe_new_predictions.csv", index=False)

else:
    print("No predictions generated.")