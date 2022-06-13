import msprime
import numpy as np
from IPython.core.display import HTML
from jupyterquiz import display_quiz

def _setup_data():
    
    
    return ts
    
class Workshop:
    css = "<style>h6 {color: white; background-color: green; padding: 4px; display: block; }</style>"

    def __init__(self):
        self.ts = self.simulate_ts()
        assert self.ts.num_mutations == self.ts.num_sites  # check it is like an infinite sites model
        self.ts.dump("simulated.trees")

    def simulate_ts(self):
        pop_size=1000
        seq_length=1_000_000
        
        sweep_model = msprime.SweepGenicSelection(
            position=seq_length/2, start_frequency=0.001, end_frequency=0.999, s=0.25, dt=1e-6)
        
        return msprime.sim_mutations(
            msprime.sim_ancestry(
                10,
                model=[sweep_model, msprime.StandardCoalescent()],
                population_size=pop_size,
                sequence_length=seq_length,
                recombination_rate=1e-8,
                random_seed=654321,  # only needed for repeatabilty
                ),
            # add finite-site mutations to the ts using the Jukes & Cantor model, creating SNPs
            rate=2e-8,
            random_seed=123456
        )

        
    def Q1(self):
        display_quiz([
            {
                "question": "How many edges in the tree sequence:",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": self.ts.num_edges,
                        "correct": True,
                        "feedback": "Correct."
                    },
                    {
                        "type": "range",
                        "range": [ -100000000, 100000], 
                        "correct": False,
                        "feedback": "Try again (hint: look at `ts.num_edges()`)"
                    },
                ]
            },
            {
                "question": "How many sites in the tree sequence:",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": self.ts.num_sites,
                        "correct": True,
                        "feedback": "Correct."
                    },
                    {
                        "type": "range",
                        "range": [ -100000000, 100000], 
                        "correct": False,
                        "feedback": "Try again (hint: look at `ts.num_sites()`)"
                    },
                ]
            },
            {
                "question": "How many mutations in the tree sequence:",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": self.ts.num_mutations,
                        "correct": True,
                        "feedback": 
                            "Correct: there are the same number of mutations "
                            "as sites, because in this tree sequence we have only "
                            "created VARIABLE sites, and have only one mutation per "
                            "site (like an infinite sites model)"
                    },
                    {
                        "type": "range",
                        "range": [ -100000000, 100000], 
                        "correct": False,
                        "feedback": "Try again (hint: look at `ts.num_sites()`)"
                    },
                ]
            },
        ])

    def Q2(self):
        display_quiz([{
            "question":
                "What is the age of the root in the first tree (to 1 d.p.)",
            "type": "numeric",
            "precision": 1,
            "answers": [
                {
                    "type": "value",
                    "value": round(self.ts.node(self.ts.first().root).time, 1),
                    "correct": True,
                    "feedback": "Correct."
                },
                {
                    "type": "range",
                    "range": [ -100000000, 1000000], 
                    "correct": False,
                    "feedback":
                        "Try again (hint: the root node in the first tree "
                        f"has ID {self.ts.first().root};"
                        " feed that to the `ts.node()` method)"
                },
            ]
        }])

    def Q3(self):
        tree = self.ts.at(400_000)
        mut_counts = np.bincount([tree.num_samples(m.node) for m in tree.mutations()])
        display_quiz([
            {
                "question":
                    "How many mutations in total are there in the tree at position 400Kb",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": tree.num_mutations,
                        "correct": True,
                        "feedback": "Correct."
                    },
                    {
                        "type": "range",
                        "range": [ -1000000, 1000000], 
                        "correct": False,
                        "feedback":
                            "Try again"
                    },
                ]
            },
            {
                "question":
                    "How many singleton mutations (above one sample, i.e. on a terminal branch)",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": int(mut_counts[1]),
                        "correct": True,
                        "feedback": "Correct."
                    },
                    {
                        "type": "range",
                        "range": [ -1000000, 1000000], 
                        "correct": False,
                        "feedback":
                            "Try again"
                    },
                ]
            },
            {
                "question":
                    "How many doubleton mutations (above 2 samples)",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": int(mut_counts[2]),
                        "correct": True,
                        "feedback": "Correct."
                    },
                    {
                        "type": "range",
                        "range": [ -1000000, 1000000], 
                        "correct": False,
                        "feedback":
                            "Try again"
                    },
                ]
            },
        ])

def setup():
    return Workshop()
        

