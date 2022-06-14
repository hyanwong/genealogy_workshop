import msprime
import tskit
import numpy as np
from IPython.core.display import HTML
from jupyterquiz import display_quiz

def _setup_data():
    
    
    return ts
    
class Workshop:
    css = """<style>
        dl {border: green 1px solid; margin-top: 1em}
        dt {color: white; background-color: green; padding: 4px; display: block; }
        dd {padding: 4px;}
    </style>"""

    # Used for making SVG formatting smaller
    small_class = "x-lab-sml"
    small_style = (
        ".x-lab-sml .sym {transform:scale(0.6)} .x-lab-sml .lab {font-size:7pt;}" # All labels small
        ".x-lab-sml .x-axis .tick .lab {"
        "font-weight:normal;transform:rotate(90deg);text-anchor:start;dominant-baseline:central;}"
    )


    def __init__(self):
        self.ts = self.simulate_ts()
        assert len(self.ts.site(12).mutations) == 2  # check there are sites with multiple mutations
        self.ts.dump("simulated.trees")

    def simulate_ts(self):
        pop_size=1000
        seq_length=1_000_000
        
        sweep_model = msprime.SweepGenicSelection(
            position=seq_length/2, start_frequency=0.001, end_frequency=0.999, s=0.25, dt=1e-6)
        
        ts = msprime.sim_ancestry(
            10,
            model=[sweep_model, msprime.StandardCoalescent()],
            population_size=pop_size,
            sequence_length=seq_length,
            recombination_rate=1e-8,
            random_seed=654321,  # only needed for repeatabilty
        )
        ts = msprime.sim_mutations(ts, rate=2e-8, random_seed=203)
        tables = ts.dump_tables()
        tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
        tables.individuals.packset_metadata([
            tables.individuals.metadata_schema.validate_and_encode_row({"name": n})
            for n in ["Ada", "Bob", "Cat", "Dee", "Eli", "Fi", "Guy", "Hal", "Ida", "Jo"]
        ])
        return tables.tree_sequence()
        
    def Q1(self):
        display_quiz([
            {
                "question": "How many edges in this tree sequence?",
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
                "question": "How many sites in this tree sequence?",
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
                "question": "How many mutations in this tree sequence?",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": self.ts.num_mutations,
                        "correct": True,
                        "feedback": 
                            "Correct: there may be more mutations than sites because "
                            "there could be multiple mutations at a single site"
                    },
                    {
                        "type": "range",
                        "range": [ -100000000, 100000], 
                        "correct": False,
                        "feedback": "Try again (hint: look at `ts.num_mutations()`)"
                    },
                ]
            },
        ])

    def Q2(self):
        display_quiz([
            {
                "question": "How many mutations at site 11?",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": len(self.ts.site(11).mutations),
                        "correct": True,
                        "feedback": "Correct."
                    },
                    {
                        "type": "range",
                        "range": [ -100000000, 100000], 
                        "correct": False,
                        "feedback": "Try again (hint: look at `ts.site(11)`)"
                    },
                ]
            },
            {
                "question": "How many mutations at site 12:",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": len(self.ts.site(12).mutations),
                        "correct": True,
                        "feedback": "Correct."
                    },
                    {
                        "type": "range",
                        "range": [ -100000000, 100000], 
                        "correct": False,
                        "feedback": "Try again (hint: look at `ts.site(12)`)"
                    },
                ]
            },
        ])

    def Q3(self):
        tree = self.ts.at(400_000)
        mut_counts = np.bincount([tree.num_samples(m.node) for m in tree.mutations()])
        display_quiz([
            {
                "question":
                    "How many mutations in total are there in the tree at position "
                    "400 Kb?",
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
                    "How many singleton mutations (above one sample, i.e. on a terminal "
                    "branch) are there in the tree at position 400 Kb?",
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
                    "How many doubleton mutations (above 2 samples) are there in the "
                    "tree at position 400 Kb?",
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

    def Q4(self):
        display_quiz([{
            "question":
                "What is the age (to the nearest generation) of the root in the first "
                "tree?",
            "type": "numeric",
            "precision": 0,
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
                        "Try again (hint: the root node in the tree "
                        f"has ID {self.ts.first().root};"
                        " feed that to the `ts.node()` method)"
                },
            ]
        }])

    def Q5(self):
        correct_name = self.ts.individual(self.ts.node(14).individual).metadata["name"]
        display_quiz([{
            "question":
                "What is the name of the individual associated with sample node 14 "
                "in the original tree sequence?",
            "type": "multiple_choice",
            "answers": [
                {
                    "answer": i.metadata["name"],
                    "correct": i.metadata["name"] == correct_name,
                    "feedback": (
                        "Correct" if i.metadata["name"] == correct_name
                        else "Sorry, that's not right."
                        )
                }
                for i in self.ts.individuals()
            ]
        }])

    def Q6(self):
        for i, v in enumerate(self.ts.variants()):
            if i == 0:
                a1 = int(v.genotypes[0])
            if i == self.ts.num_sites-1:
                a2 = int(v.genotypes[1])

        display_quiz([
            {
                "question":
                    "What is the genotypic state of the first sample at the first site?",
                "type": "multiple_choice",
                "answers": [
                    {
                        "answer": f"{s}", "correct": (a1 == s),
                        "feedback": ("Correct" if a1 == s else "Sorry, that's not right.")
                    }
                    for s in [0, 1]
                ]
            },
            {
                "question":
                    "What is the genotypic state of the second sample at the last site?",
                "type": "multiple_choice",
                "answers": [
                    {
                        "answer": f"{s}", "correct": (a2 == s),
                        "feedback": ("Correct" if a2 == s else "Sorry, that's not right.")
                    }
                    for s in [0, 1]
                ]
            },
        ])



def setup():
    return Workshop()
        

